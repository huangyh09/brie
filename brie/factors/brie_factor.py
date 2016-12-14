# This file is to obtain the sequences for introns.
# Features: 
# 1. frameshift A: whether A is a multiple of 3
# 2. Translatable.C1: existance of stop codon UAA, UAG, UGA
# 3. Translatable.C1-C2
# 4. Translatable.C1-A
# 5. Translatable.C1-A-C2
# 6. AltAGpos: Distance to the first alternative AG in I1_3p from its 5'ss
# 7. AltGTpos: Distance to the first alternative GT in I2_5p from its 3'ss

# Yuanhua Huang, 2016/10/04


import sys
import time
import h5py
import itertools
import subprocess
import numpy as np
import multiprocessing
from optparse import OptionParser
from diceseq import load_annotation, FastaFile

PROCESSED = 0
TOTAL_GENE = 0
START_TIME = time.time()

def show_progress(RV=None):
    global PROCESSED, TOTAL_GENE, START_TIME
    
    PROCESSED += 1
    bar_len = 30
    run_time = time.time() - START_TIME 
    percents = 100.0 * PROCESSED / TOTAL_GENE
    filled_len = int(round(bar_len * percents / 100))
    bar = '=' * filled_len + '-' * (bar_len - filled_len)
    
    sys.stdout.write('\r[%s] %.2f%% processed in %.1f sec.' 
        % (bar, percents, run_time))
    sys.stdout.flush()
    return RV


def rev_seq(seq):
    _tmp = []
    _tmp[:] = seq
    for j in range(len(_tmp)):
        if _tmp[j] == "A": _tmp[j] = "T"
        elif _tmp[j] == "T": _tmp[j] = "A"
        elif _tmp[j] == "G": _tmp[j] = "C"
        elif _tmp[j] == "C": _tmp[j] = "G"
    RV = "".join(_tmp[::-1])
    return RV

def DNA2RNA(seq):
    _tmp = []
    _tmp[:] = seq
    for j in range(len(_tmp)):
        if _tmp[j] == "A": _tmp[j] = "U"
        elif _tmp[j] == "T": _tmp[j] = "A"
        elif _tmp[j] == "G": _tmp[j] = "C"
        elif _tmp[j] == "C": _tmp[j] = "G"
    RV = "".join(_tmp)
    return RV

def cDNA2RNA(seq):
    _tmp = []
    _tmp[:] = seq
    for j in range(len(_tmp)):
        if _tmp[j] == "T": _tmp[j] = "U"
    RV = "".join(_tmp)
    return RV

def get_motif(seq_full, motif, mode="counts"):
    """get the counts of motif in a sequence"""
    cnt = 0
    for i in range(len(seq_full)-len(motif)+1):
        if seq_full[i:i+len(motif)] == motif:
            cnt += 1
    if mode == "counts":
        return cnt
    elif mode == "frequency":
        return cnt / (len(seq_full)-len(motif)+1.0)
    elif mode == "normalized":
        return cnt / (len(seq_full)-len(motif)+1.0) / (0.25**len(motif))
    else:
        return None


def get_kmer_all(kmax=5, kmin=1, seqs="ATGC"):
    RV = []
    for i in range(kmin, kmax+1):
        for _seq in itertools.product(seqs, repeat=i): 
            RV.append("".join(_seq))
    return RV


def get_factor(tran, ref_file, phast_file):
    """
    tran: a triplets of skippint exon, i.e., exon1-AS_exon-exon3;
    ref_file: the file name of the genome reference in fasta format;
    factors: N*5-like array, with number, motif, region, methods, and
    comments. 
        -region: 7 regions as Barash et al 2010, Nature.
        -methods: logLen, 2ndStr, counts, frequency.
    """
    RV = np.zeros(695)
    if tran.exonNum != 3:
        print("This is not a triplet of exons. Please check.")
        return RV
    else:
        exons = tran.exons
        chrom = tran.chrom
    fastaFile = FastaFile(ref_file)

    # logLength #8
    logLen = [np.log(exons[0,1]-exons[0,0]+1),  #C1
              np.log(exons[1,0]-exons[0,1]-1),  #I1
              np.log(exons[1,1]-exons[1,0]+1),  #A
              np.log(exons[2,0]-exons[1,1]-1),  #I2
              np.log(exons[2,1]-exons[2,0]+1)]  #C2
    if tran.strand != "+" and tran.strand != "1":
        logLen = logLen[::-1]
    logLen += [logLen[2]/logLen[1], logLen[2]/logLen[3], logLen[1]/logLen[3]]
        
    # Splice site score #4
    SS_end = [5, 3, 5, 3]
    SS_seq = [fastaFile.get_seq(chrom, exons[0,1]-3,  exons[0,1]+8), #I1'5
              fastaFile.get_seq(chrom, exons[1,0]-17, exons[1,0]+3), #I1'3
              fastaFile.get_seq(chrom, exons[1,1]-3,  exons[1,1]+8), #I2'5
              fastaFile.get_seq(chrom, exons[2,0]-17, exons[2,0]+3)] #I2'3
    if tran.strand != "+" and tran.strand != "1":
        SS_end = SS_end[::-1]
        SS_seq = SS_seq[::-1]
        SS_seq = [rev_seq(x) for x in SS_seq]
    SS_val = splice_site_scores(SS_seq, mode=SS_end)

    # PhastCons score #7
    regID = ["C1", "I1_5p", "I1_3p", "A", "I2_5p", "I2_3p", "C2"]
    regions = [[exons[0,0], exons[0,1]],        #C1
               [exons[0,1]+1, exons[0,1]+300],  #I1_5p
               [exons[1,0]-300, exons[1,0]-1],  #I1_3p
               [exons[1,0], exons[1,1]],        #A
               [exons[1,1]+1, exons[1,1]+300],  #I2_5p
               [exons[2,0]-300, exons[2,0]-1],  #I2_3p
               [exons[2,0], exons[2,1]]]        #C2
    if tran.strand != "+" and tran.strand != "1":
        regions = regions[::-1]

    cons_val = []
    for i in range(len(regions)):
        bashCommand = "bigWigSummary %s %s %d %d 1" %(phast_file, chrom, 
            regions[i][0], regions[i][1]) 
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        output = process.communicate()[0]
        try:
            cons_val.append(float(output))
        except ValueError:
            cons_val.append(0.0) #np.nan
            print "No PhastCons data for %s. Treated as Zero." %tran.tranID
            # print output

    # 1 to 4mers #676
    regK = [3, 2, 3, 4, 3, 2, 3]
    kmer_frq = []
    for i in range(len(regions)):
        kmers = get_kmer_all(kmax=regK[i], kmin=1, seqs="ATGC")
        for i in range(len(regions)):
            seq = fastaFile.get_seq(chrom, regions[i][0], regions[i][1])
            if tran.strand != "+" and tran.strand != "1":
                seq = rev_seq(seq)
            for j in range(len(kmers)):
                kmer_frq.append(get_motif(seq, kmers[j], mode="frequency"))

    RV = np.array(logLen + SS_val + cons_val + kmer_frq)
    return RV


def main():
    print("Welcome to Genomic feature extactor!")

    #part 0. parse command line options
    parser = OptionParser()
    parser.add_option("--anno_file", "-a", dest="anno_file", default=None,
        help="Annotation file of genes")
    parser.add_option("--anno_type", dest="anno_type", default="GTF",
        help="Type of annotation file [default: %default]." )
    parser.add_option("--ref_seq", "-r", dest="ref_seq", default=None,
        help="Genome sequence reference in fasta file.")
    parser.add_option("--phastCons", "-c", dest="phast_file", default=None,
        help="PhastCons conservation scores in bigWig file.")
    parser.add_option("--out_file", "-o", dest="out_file",  
        default="splicing_factor.h5", help="Output in hdf5 file")
    parser.add_option("--nproc", "-p", dest="nproc", type="int", default="1",
        help="The number of subprocesses [default: %default].")

    (options, args) = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        print("use -h or --help for help on argument.")
        sys.exit(1)
    if options.anno_file is None:
        print("Error: need --anno_file for annotation.")
        sys.exit(1)
    else:
        sys.stdout.write("\rloading annotation file...")
        sys.stdout.flush()    
        anno = load_annotation(options.anno_file, options.anno_type)
        sys.stdout.write("\rloading annotation file... Done.\n")
        sys.stdout.flush()
        genes = anno["genes"]
    if options.ref_seq is None:
        print("Error: need --ref_seq for genome reference.")
        sys.exit(1)
    else:
        ref_file = options.ref_seq
    if options.phast_file is None:
        print("No bigWig file for phastCons, ignore these features.")
    else:
        phast_file = options.phast_file
    if options.out_file is None:
        out_file = os.path.dirname(os.path.abspath(ref_file)) + "/brieFactor.h5"
    else:
        out_file = options.out_file
    nproc = options.nproc

    gene_ids = []
    for g in genes:
        gene_ids.append(g.geneID)
    feature_all = np.zeros((len(gene_ids), 9*256))
    
    global TOTAL_GENE
    TOTAL_GENE = len(genes)

    print("extracting features for %d skipping exon triplets with "
          "%d cores..." %(TOTAL_GENE,  nproc))

    if nproc <= 1:
        for g in range(len(genes)):
            feature_all[g,:] = get_factor(genes[g].trans[0], ref_file, 
                phast_file)
            show_progress()
    else:
        pool = multiprocessing.Pool(processes=nproc)
        result = []
        for g in genes:
            result.append(pool.apply_async(get_factor, (g.trans[0], ref_file, 
                phast_file), callback=show_progress))
        pool.close()
        pool.join()
        for g in range(len(result)):
            feature_all[g,:] = result[g].get()
    print("")

    regID = ["C1", "I1_5p", "I1", "I1_3p", "A", "I2_5p", "I2", "I2_3p", "C2"]
    factors = ["Cons.%s" %x for x in regID]

    f = h5py.File(out_file, "w")
    f["factors"] = factors
    f["gene_ids"] = gene_ids
    f["features"] = feature_all
    f.close()


if __name__ == '__main__':
    main()
