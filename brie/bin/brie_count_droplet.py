# This function is to count reads supporting each read category

import os
import sys
import pysam
import numpy as np
from optparse import OptionParser, OptionGroup

from ..version import __version__
from ..utils.count import SE_effLen
from ..utils.count_droplet import get_droplet_matrix
from ..utils.io_utils import read_brieMM, read_gff, convert_to_annData


def droplet_count(gff_file, sam_file, barcode_file, out_dir=None, nproc=1, 
                  event_type='SE', CB_tag='CB', UMI_tag='UR'):
    """CLI for counting reads supporting isoforms
    """
    ## TODO: Parameter check
    
    ## Load sam file list
    cell_list = np.loadtxt(barcode_file, delimiter = None, dtype=str, ndmin = 2)
    cell_list = cell_list[:, 0]
    print(cell_list[:min(3, cell_list.shape[0])], "...")

    ## Check out_dir
    if out_dir is None:
        sam_dir = os.path.abspath(sam_file)
        out_dir = os.path.dirname(sam_dir) + "/brieCOUNT"
    else:
        out_dir = out_dir
        
    try:
        os.stat(os.path.abspath(out_dir))
    except:
        os.mkdir(os.path.abspath(out_dir))
        
    ## Load GFF file
    sys.stdout.write("\r[BRIE2] loading gene annotations ... ")
    sys.stdout.flush()
    genes = read_gff(gff_file)
    sys.stdout.write("\r[BRIE2] loading gene annotations ... Done.\n")
    sys.stdout.flush()
    
    ## Output gene info
    gene_table = [["GeneID", "GeneName", "TranLens", "TranIDs"]]
    for g in genes:
        tran_ids, tran_lens = [], []
        for t in g.trans:
            tran_ids.append(t.tranID)
            tran_lens.append(str(t.tranL))
        out_list = [g.geneID, g.geneName, ",".join(tran_lens), 
                    ",".join(tran_ids)]
        gene_table.append(out_list)
        
    fid = open(out_dir + "/gene_note.tsv", "w")
    for out_list in gene_table:
        fid.writelines("\t".join(out_list) + "\n")
    fid.close()
        
    ## Output sam total counts
    total_reads = 0
    if not os.path.isfile(sam_file):
        print("Error: No such file\n    -- %s" %sam_file)
        sys.exit(1)
    pysam_stats = pysam.idxstats(sam_file)
    if type(pysam_stats) is not list:
        pysam_stats = pysam_stats.split("\n")
    for tp in pysam_stats: 
        tmp = tp.strip().split("\t")
        if len(tmp) >= 3:
            total_reads += float(tmp[2])
                
    ### Output cell info
    fid = open(out_dir + "/cell_note.tsv", "w")
    for i in range(len(cell_list)):
        fid.writelines("%s\n" %(cell_list[i]))
    fid.close()
    
    ## Generate isoform effective length matrix
    if event_type == 'SE':
        effLen_tensor = np.zeros((len(genes), 2, 3), dtype=np.float32)
        for ig, _gene in enumerate(genes):
            effLen_tensor[ig, :, :] = SE_effLen(_gene, rlen=76)
    else:
        # placeholder, not support yet
        effLen_tensor = np.zeros((len(genes), 1), dtype=np.float32)
    
    ## Load read counts
    print("[BRIE2] counting reads for %d genes in %d cells with %d cores..." 
          %(len(genes), cell_list.shape[0], nproc))

    res = get_droplet_matrix(genes, sam_file, cell_list, out_dir, 
                             event_type, 10, 2, CB_tag, UMI_tag)
    
    ## Don't save into h5ad if not SE
    if event_type != 'SE':
        sys.exit()
    
    ## Save data into h5ad
    print("[BRIE2] save matrix into h5ad ...")
    Rmat_dict = read_brieMM(out_dir + "/read_count.mtx")
    adata = convert_to_annData(Rmat_dict=Rmat_dict, 
                               effLen_tensor=effLen_tensor, 
                               cell_note=np.array(cell_list, dtype='str'), 
                               gene_note=np.array(gene_table, dtype='str'))
    adata.uns['event_type'] = event_type
    adata.uns['total_reads'] = total_reads
    adata.write_h5ad(out_dir + "/brie_count.h5ad")

    
def main():
    # import warnings
    # warnings.filterwarnings('error')

    # parse command line options
    parser = OptionParser()
    parser.add_option("--gffFile", "-a", dest="gff_file", default=None,
        help="GTF/GFF3 file for gene and transcript annotation")
    parser.add_option("--samFile", "-s", dest="sam_file", default=None,
        help=("An indexed bam/sam/cram files"))
    parser.add_option("--barcodes", "-b", dest="barcodes_file", default=None,
        help=("A file containing cell barcodes without header"))
    parser.add_option("--out_dir", "-o", dest="out_dir", default=None, 
        help="Full path of output directory [default: $samFile/brieCOUNT]")

    group = OptionGroup(parser, "Optional arguments")
    group.add_option("--nproc", "-p", type="int", dest="nproc", default="4",
        help="Number of subprocesses [default: %default]")
    group.add_option("--eventType", "-t", dest="event_type", default="SE",
        help="Type of splicing event for check. SE: skipping-exon; "
             "Any: no-checking [default: %default]")
    group.add_option("--cellTAG", dest="cell_tag", default="CB",
        help="Tag for cell barocdes [default: %default]")
    group.add_option("--UMItag", dest="UMI_tag", default="UR",
        help="Tag for UMI barocdes [default: %default]")
    
    # group.add_option("--add_premRNA", action="store_true", dest="add_premRNA", 
    #     default=False, help="Add the pre-mRNA as a transcript")

    parser.add_option_group(group)
    
    (options, args) = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        print("Welcome to brie-count in BRIE v%s!\n" %(__version__))
        print("use -h or --help for help on argument.")
        sys.exit(1)
    if options.sam_file == None:
        print("[BRIE2] Error: need --samFile for indexed & aliged sam/bam/cram file.")
        sys.exit(1)
        
    if options.gff_file == None:
        print("[BRIE2] Error: need --gffFile for gene annotation.")
        sys.exit(1)
    
    # bias mode is not supported yet
    add_premRNA = False
    
    droplet_count(options.gff_file, options.sam_file, options.barcodes_file, 
        options.out_dir, options.nproc, options.event_type, options.cell_tag,
        options.UMI_tag)


if __name__ == "__main__":
    main()
    
