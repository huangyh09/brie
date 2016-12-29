# This file is to preprocess the annotation file to keep the alternative
# exons with high quality, by using the following criteria:
# 1. surronding introns are no shorter than 100bp
# 2. not overlapped by any SE triplets or AS-exon
# 3. the length of the alternative exon within a given region, e.g. 50~450bp
# 4. with a minimum distance from TSS or TTS, e.g., 500
# 5. surrounded by AG-GT, i.e., AG-as_exon-GT
# 6. located on specific chromosomes, i.e., 1-22 and X

# To do so, we first choose the chromosomes, and then filter genes overlaping 
# with others, then only keep the genes with given biotype. Finally, based on
# those genes, we will include the exons with high quality.

import numpy as np
from diceseq import FastaFile
from optparse import OptionParser

def get_gene_idx(anno_in):
    """To get the index of all genes in the annotation lines.
    The index number starts from 0."""
    g_idx = []
    now_g = -1
    g_chr = []
    g_start = []
    g_stop = []
    for i in range(len(anno_in)):
        if anno_in[i][0] == "#" or anno_in[i][0] == ">" : continue
        a_line = anno_in[i].split("\t")
        if len(a_line) >= 9: last_g = i
        if len(a_line) >= 9 and a_line[2] == 'gene':
            pre_g = now_g + 0
            now_g = i
            if pre_g != -1: 
                g_idx.append([pre_g, now_g-1])
                g_chr.append(anno_in[pre_g].split("\t")[0])
                g_start.append(int(anno_in[pre_g].split("\t")[3]))
                g_stop.append(int(anno_in[pre_g].split("\t")[4]))
    g_idx.append([now_g, last_g])
    g_chr.append(anno_in[now_g].split("\t")[0])
    g_start.append(int(anno_in[now_g].split("\t")[3]))
    g_stop.append(int(anno_in[now_g].split("\t")[4]))
    return np.array(g_idx), np.array(g_chr), np.array(g_start), np.array(g_stop)


def gene_overlap_check(anno_in, g_idx, g_chr, mode="exon"):
    """check wheter overlap between genes or AS-exons."""
    start_loc, stop_loc = [], []
    if mode == "gene":
        for i in g_idx[:,0]:
            start_loc.append(int(anno_in[i].split("\t")[3]))
            stop_loc.append(int(anno_in[i].split("\t")[4]))
    elif mode == "exon":
        for i in g_idx[:,0]:
            start_loc.append(int(anno_in[i+3].split("\t")[3]))
            stop_loc.append(int(anno_in[i+3].split("\t")[4]))
    start_loc, stop_loc = np.array(start_loc), np.array(stop_loc)

    g_keep = np.ones(len(g_idx), "bool")
    chr_unique = np.unique(g_chr)
    for c in chr_unique:
        idx = np.where(g_chr == c)[0]
        sort_idx = idx[np.argsort(start_loc[idx])]
        for i in range(len(sort_idx)-1):
            if start_loc[sort_idx[i+1]] <= stop_loc[sort_idx[i]]:
                # print "test1"
                # for temp in anno_in[g_idx[sort_idx[i+1],0]:g_idx[sort_idx[i+1],1]+1]:
                #     print temp.split()[:-1]
                # print "test2"
                # for temp in anno_in[g_idx[sort_idx[i],0]:g_idx[sort_idx[i],1]+1]:
                #     print temp.split()[:-1]
                g_keep[sort_idx[i]] = False
                g_keep[sort_idx[i+1]] = False
    anno_out = []
    for i in range(len(g_idx)):
        if g_keep[i]:
            anno_out += anno_in[g_idx[i,0]: g_idx[i,1]+1]
    return anno_out


def as_exon_check(fastaFile, anno_in, g_idx, as_exon_min, as_exon_max, 
    as_exon_tss, as_exon_tts, chroms):
    """check the quality of alternative exon."""
    
    anno_out = []
    for i in range(g_idx.shape[0]):
        vals_g = anno_in[g_idx[i,0]].split("\t")
        _exon_loc = np.array(anno_in[g_idx[i,0]+3].split("\t")[3:5], "int")
        _exon1_loc = np.array(anno_in[g_idx[i,0]+2].split("\t")[3:5], "int")
        _exon3_loc = np.array(anno_in[g_idx[i,0]+4].split("\t")[3:5], "int")

        # in case, reverse in minus strand
        if _exon1_loc[0] > _exon3_loc[0]:
            _exon1_loc, _exon3_loc = _exon3_loc, _exon1_loc

        _exon_len = _exon_loc[1] - _exon_loc[0] + 1

        # check 1: not too short or too long
        if _exon_len < as_exon_min or _exon_len > as_exon_max: 
            continue
        
        # check 2: chromsome
        chrom = vals_g[0]
        if chroms.count(chrom) == 0:
            continue

        # check 3: surrounding splice sites AG--exon--GT 
        if vals_g[6] == "+":
            up_ss3 = fastaFile.get_seq(chrom, _exon_loc[0]-2, _exon_loc[0]-1)
            dn_ss5 = fastaFile.get_seq(chrom, _exon_loc[1]+1, _exon_loc[1]+2)
            if up_ss3 != "AG" or dn_ss5 != "GT":
                continue
            # print up_ss3, dn_ss5
        else:
            up_ss3 = fastaFile.get_seq(chrom, _exon_loc[1]+1, _exon_loc[1]+2)
            dn_ss5 = fastaFile.get_seq(chrom, _exon_loc[0]-2, _exon_loc[0]-1)
            if up_ss3 != "CT" or dn_ss5 != "AC": 
                continue
            # print up_ss3, dn_ss5

        # check 4: not too close to TSS or TTS
        # and not too short introns.
        if vals_g[6] == "+":
            tss_dis = _exon_loc[0] - _exon1_loc[0]
            tts_dis = _exon3_loc[1] - _exon_loc[1]
            up_dis = _exon_loc[0] - _exon1_loc[1]
            dn_dis = _exon3_loc[0] - _exon_loc[1]
        else:
            tts_dis = _exon_loc[0] - _exon1_loc[0]
            tss_dis = _exon3_loc[1] - _exon_loc[1]
            dn_dis = _exon_loc[0] - _exon1_loc[1]
            up_dis = _exon3_loc[0] - _exon_loc[1]
        if tts_dis <= 0 and tss_dis <= 0:
            print "TTS or TSS distance warning: %s" %vals_g[6]
            for temp in anno_in[g_idx[i,0]:g_idx[i,1]+1]:
                print temp.split()[:-1]

        if tss_dis < as_exon_tss or tts_dis < as_exon_tts:
            # print tss_dis, tts_dis
            continue

        # intron distance
        if up_dis < 100 or dn_dis < 100:
            continue

        # high quality
        anno_out += anno_in[g_idx[i,0]:g_idx[i,1]+1]

    g_num, as_exon_num = len(anno_out) / 8, len(anno_out) / 8
    return anno_out, g_num, as_exon_num

def map_ids(id1, id2):
    """this function is to map id2 to id1"""
    j = 0
    idx2 = []
    sidx1 = np.argsort(id1)
    sidx2 = np.argsort(id2)
    sidx0 = np.argsort(sidx1)

    for i in range(len(sidx1)):
        # if i >= 1 and id1[sidx1[i-1]] == id1[sidx1[i]]:
        #     idx2.append(idx2[-1])
        #     continue
        if j >= len(id2):
            idx2.append(None)
            continue
        while id1[sidx1[i]] > id2[sidx2[j]]:
            j += 1
        if id1[sidx1[i]] == id2[sidx2[j]]:
            idx2.append(sidx2[j])
            j += 1
            continue
        else: 
            print "id mapping warning:"
            print id1[sidx1[i-1]], id2[sidx2[j-1]]
            print id1[sidx1[i]], id2[sidx2[j-1]]
            print id1[sidx1[i]], id2[sidx2[j]]
            
            idx2.append(None)
            continue
    return np.array(idx2)[sidx0]


def save_out(anno_in, anno_ref, out_file):
    exon_str_SE = []
    for i in range(0, len(anno_in), 8):
        vals = anno_in[i+3].strip().split("\t")
        exon_str_SE.append(".".join([vals[0],vals[3],vals[4],vals[6]]))

    exon_str_ref = []
    gene_info = []
    for i in range(len(anno_ref)):
        if anno_ref[i].startswith("#"):continue
        vals = anno_ref[i].strip().split("\t")
        if vals[2] != "exon":
            continue
        exon_str_ref.append(".".join([vals[0],vals[3],vals[4],vals[6]]))
        gene_info.append(vals[8])

    ginfo_idx = map_ids(exon_str_SE, exon_str_ref)

    # print len(exon_str_SE), len(exon_str_ref)
    # print len(ginfo_idx), len(gene_info), len(anno_in)
    #exit()

    temp_gene = []
    fid = open(out_file, "w")
    fid.writelines("#annotation file with high-qulity alternative exons.\n")
    for i in range(0, len(anno_in), 8):
    	_gene_id = "#"
        if i is not None:
            # idx  = gene_info[i].find("gene_id")
            idx  = gene_info[ginfo_idx[i/8]].find("gene_id")
            if idx > -1:
                # _gene_id = gene_info[i][idx:].split('"')[1].split("\n")[0]
                _gene_id = gene_info[ginfo_idx[i/8]][idx:].split('"')[1].split("\n")[0]
            _gene_id = _gene_id.split(".")[0]
        _num = temp_gene.count(_gene_id)
        temp_gene.append(_gene_id)
        if _num > 0:
        	_gene_id += ".AS%d" %(_num+1)

        if _gene_id.startswith("#"): 
            if anno_in[i].split()[6] != "-": 
                print "test"
                for temp in anno_in[i:i+8]:
                    print temp.split()[:-1] 
                continue

        ## In gtf format
        vals = anno_in[i+0].strip().split("\t")
        vals[8] = "gene_id \"%s\"" %(_gene_id)
        fid.writelines("\t".join(vals) + "\n")

        vals = anno_in[i+1].strip().split("\t")
        vals[2] = "transcript"
        vals[8] = "gene_id \"%s\"; transcript_id \"%s.in\"" %(_gene_id, _gene_id)
        fid.writelines("\t".join(vals) + "\n")

        vals = anno_in[i+2].strip().split("\t")
        vals[8] = "gene_id \"%s\"; transcript_id \"%s.in\"" %(_gene_id, _gene_id)
        fid.writelines("\t".join(vals) + "\n")

        vals = anno_in[i+3].strip().split("\t")
        vals[8] = "gene_id \"%s\"; transcript_id \"%s.in\"" %(_gene_id, _gene_id)
        fid.writelines("\t".join(vals) + "\n")

        vals = anno_in[i+4].strip().split("\t")
        vals[8] = "gene_id \"%s\"; transcript_id \"%s.in\"" %(_gene_id, _gene_id)
        fid.writelines("\t".join(vals) + "\n")


        vals = anno_in[i+5].strip().split("\t")
        vals[2] = "transcript"
        vals[8] = "gene_id \"%s\"; transcript_id \"%s.out\"" %(_gene_id, _gene_id)
        fid.writelines("\t".join(vals) + "\n")

        vals = anno_in[i+6].strip().split("\t")
        vals[8] = "gene_id \"%s\"; transcript_id \"%s.out\"" %(_gene_id, _gene_id)
        fid.writelines("\t".join(vals) + "\n")

        vals = anno_in[i+7].strip().split("\t")
        vals[8] = "gene_id \"%s\"; transcript_id \"%s.out\"" %(_gene_id, _gene_id)
        fid.writelines("\t".join(vals) + "\n")

    fid.close()


def main():
    #part 0. parse command line options
    parser = OptionParser()
    parser.add_option("--anno_file", "-a", dest="anno_file", default=None,
        help="The annotation file of SE events in gff3 format from rnaseqlib.")
    parser.add_option("--anno_ref", dest="anno_ref", default=None,
        help="The reference annotation file in gtf format.")
    parser.add_option("--reference", "-r", dest="reference", default=None,
        help="The genome reference sequence file in fasta format.")
    parser.add_option("--out_file", "-o", dest="out_file", default=None,
        help="The prefix of out files.")

    parser.add_option("--as_exon_min", dest="as_exon_min", default="50",
        help="the minimum length for the alternative splicing exon.")
    parser.add_option("--as_exon_max", dest="as_exon_max", default="450",
        help="the maximum length for the alternative splicing exon.")
    parser.add_option("--as_exon_tss", dest="as_exon_tss", default="500",
        help="the minimum length for the alternative exon to TSS.")
    parser.add_option("--as_exon_tts", dest="as_exon_tts", default="500",
        help="the minimum length for the alternative exon to TTS.")

    parser.add_option("--add_chrom", dest="add_chrom", default="chrX",
        help="the extra chromosomes besides autosome, e.g., chrX,chrY,chrM")

    (options, args) = parser.parse_args()
    if options.anno_file is None:
        print("Error: need --anno_file for annotation.")
        sys.exit(1)
    else:
        fid = open(options.anno_file, "r")
        anno_in = fid.readlines()
        fid.close()
    if options.out_file is None:
        out_file = ".".join(options.anno_file.split(".")[:-1])
    else:
        out_file = options.out_file
    if options.reference is None:
        print("Error: need --reference for genome sequecne.")
        sys.exit(1)
    else:
        fastaFile = FastaFile(options.reference)

    if options.anno_ref is None:
        anno_ref = None
    else:
        fid = open(options.anno_ref, "r")
        anno_ref = fid.readlines()
        fid.close()
    
    add_chrom = options.add_chrom
    as_exon_min = int(options.as_exon_min)
    as_exon_max = int(options.as_exon_max)
    as_exon_tss = int(options.as_exon_tss)
    as_exon_tts = int(options.as_exon_tts)

    chroms = []
    for i in range(1,23):
        chroms.append("chr%d" %i)
    chroms += add_chrom.split(",")

    # remove overlap splicing events
    # g_idx, g_chr, g_start, g_stop = get_gene_idx(anno_in)
    # anno_in = gene_overlap_check(anno_in, g_idx, g_chr, g_start, g_stop)

    # get gene index
    g_idx, g_chr, g_start, g_stop = get_gene_idx(anno_in)
    print("%d Skipped Exon events are input for quality check." %(len(g_idx)))

    # alternative exon quality check
    anno_out, g_num, ex_num = as_exon_check(fastaFile, anno_in, g_idx, 
        as_exon_min, as_exon_max, as_exon_tss, as_exon_tts, chroms)
    print("%d Skipped Exon events pass the qulity control." %(g_num))
    
    # remove overlapped skipping exons
    g_idx, g_chr, g_start, g_stop = get_gene_idx(anno_out)
    anno_out = gene_overlap_check(anno_out, g_idx, g_chr, "exon")

    g_idx, g_chr, g_start, g_stop = get_gene_idx(anno_out)
    print("%d Skipped Exon events pass the overlapping check." %(len(g_idx)))

    # saving out
    save_out(anno_out, anno_ref, out_file+".gold.gtf")
    

if __name__ == "__main__":
    main()
    