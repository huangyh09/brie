# This file is to extact 735 sequences features for Brie:

# 1. Splice-site score #4: I1_5'ss, I1_3'ss, I2_5'ss, I2_3'ss
# 2. logLength #8: C1, I1, A, I2, C2, A/I1, A/I2, I1/I2
# 3. Conservation score #7: C1, I1_5p, I1, I1_3p, A, I2_5p, I2, I2_3p, C2
# 4. 2-4 mers for above 7 regions #716: 3mer, 2mer, 3mer, 4mer, 3mer, 2mer, 3mer

# Yuanhua Huang, 2016/11/10


import os
import sys
import time
import h5py
import numpy as np
import multiprocessing
from optparse import OptionParser, OptionGroup

# import pyximport; pyximport.install()
from utils.gtf_utils import loadgene
from utils.fasta_utils import get_factor, get_factorID, motif_score

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


def main():
    #part 0. parse command line options
    parser = OptionParser()
    parser.add_option("--anno_file", "-a", dest="anno_file", default=None,
        help="Annotation file for genes and transcripts in GTF or GFF3")
    parser.add_option("--ref_seq", "-r", dest="ref_seq", default=None,
        help="Genome sequence reference in fasta file.")
    parser.add_option("--phastCons", "-c", dest="phast_file", default=None,
        help="PhastCons conservation scores in bigWig file.")
    parser.add_option("--out_file", "-o", dest="out_file",  
        default="splicing_factor.h5", help="Output in hdf5 file")

    group = OptionGroup(parser, "Optional arguments")
    group.add_option("--nproc", "-p", type="int", dest="nproc", default=4,
        help="Number of subprocesses [default: %default]")
    group.add_option("--MSA5ss", dest="msa_5ss", default=None,
        help=("Mutiple sequence alignment file for 5'splice-site. It is from "
              "-4 to 7. As default, MSA is based on input 5 splice sites."))
    group.add_option("--MSA3ss", dest="msa_3ss", default=None,
        help=("Mutiple sequence alignment file for 3'splice-site. It is from "
              "-16 to 4. As default, MSA is based on input 3 splice sites."))
    parser.add_option_group(group)

    (options, args) = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        print("Welcome to Brie-factor extactor!\n")
        print("use -h or --help for help on argument.")
        sys.exit(1)
    if options.anno_file is None:
        print("Error: need --anno_file for annotation.")
        sys.exit(1)
    else:
        sys.stdout.write("\rloading annotation file...")
        sys.stdout.flush()    
        genes = loadgene(options.anno_file)
        sys.stdout.write("\rloading annotation file... Done.\n")
        sys.stdout.flush()
    if options.ref_seq is None:
        print("Error: need --ref_seq for genome reference.")
        sys.exit(1)
    else:
        ref_file = options.ref_seq
        if os.path.isfile(ref_file) == False:
            print("Cann't find ref_seq file: %s" %ref_file)
            sys.exit(1)
    
    if options.out_file is None:
        out_file = os.path.dirname(os.path.abspath(ref_file)) + "/brieFactor.h5"
    else:
        out_file = options.out_file
    nproc = options.nproc
    
    if options.phast_file is None:
        print("No bigWig file for phastCons, ignore these features.")
        phast_file = None
    elif os.path.isfile(options.phast_file) == False:
        print("Cann't find phast file: %s\nTurn it off." %options.phast_file)
        phast_file = None
    else:
        phast_file = options.phast_file

    if options.msa_5ss is None:
        msa_5ss = None
    elif os.path.isfile(options.msa_5ss) == False:
        print("Cann't find msa_5ss file: %s\nTurn it off." %msa_5ss)
        msa_5ss = None
    else:
        msa_5ss = np.loadtxt(options.msa_5ss, dtype="str", comments='>')

    if options.msa_3ss is None:
        msa_3ss = None
    elif os.path.isfile(options.msa_3ss) == False:
        print("Cann't find msa_3ss file: %s\nTurn it off." %msa_3ss)
        msa_3ss = None
    else:
        msa_3ss = np.loadtxt(options.msa_3ss, dtype="str", comments='>')

    gene_ids = []
    for g in genes:
        gene_ids.append(g.geneID)
    factors = get_factorID(phast_file!=None)
    # print len(factors), factors[:30], factors[-30:]
    features = np.zeros((len(gene_ids), len(factors)))
    
    global TOTAL_GENE
    TOTAL_GENE = len(genes)

    print("extracting features for %d skipping exon triplets with "
          "%d cores..." %(TOTAL_GENE,  nproc))

    SS_seq = []
    if nproc <= 1:
        for g in range(len(genes)):
            RV = get_factor(genes[g].trans[0], ref_file, 
                phast_file)
            SS_seq.append(RV["SS_seq"])
            features[g, 4:] = RV["factor_val"]
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
            RV = result[g].get()
            SS_seq.append(RV["SS_seq"])
            features[g, 4:] = RV["factor_val"]

    SS_seq = np.array(SS_seq)
    features[:, 0] = motif_score(SS_seq[:,0], msa_5ss)
    features[:, 1] = motif_score(SS_seq[:,1], msa_3ss)
    features[:, 2] = motif_score(SS_seq[:,2], msa_5ss)
    features[:, 3] = motif_score(SS_seq[:,3], msa_3ss)
    print("")

    f = h5py.File(out_file, "w")
    f.create_dataset("factors", data=factors, compression="gzip")
    f.create_dataset("gene_ids", data=gene_ids, compression="gzip")
    f.create_dataset("features", data=features, compression="gzip")
    f.close()


if __name__ == '__main__':
    main()
