# This file to simulate single cell RNA-seq reads based on real bulk RNA-seq 
# expression profile and input dropout rate and number of reads.

import os
import sys
import time
import subprocess
import numpy as np
from optparse import OptionParser, OptionGroup

# import pyximport; pyximport.install()
from utils import id_mapping

START_TIME = time.time()

def logistic(x):
    """
    Logistic function, mapping (-inf, inf) to (0,1)
    Parameters
    ----------
    x: float, int, array, list
        input variable
    Returns
    -------
    val: float, int, array
        logistic(x)
    """
    return np.exp(x)/(1+np.exp(x))
    
def logit(x, minval=0.001):
    """
    Logit function, mapping (0,1) to (-inf, inf)
    Parameters
    ----------
    x: float, int, array, list
        input variable
    minval: float (optional, default=0.001)
        minimum value of input x, and maximum of 1-x
    Returns
    -------
    val: float, int, array
        logit(x)
    """
    if isinstance(x,  (list, tuple, np.ndarray)):
        x[1-x<minval] = 1-minval
        x[x<minval] = minval
    else:
        x = max(minval, x)
        x = min(1-minval, x)
    val = np.log(x/(1-x))
    return val

def adjust_drop_prob(drop_prob, rate_new=0.3):
    """
    Adjust the drop-out rate based on the input 
    drop-out probability profile.
    
    Parameters:
    -----------
    drop_prob: array like
        the drop-out probability distribution 
    rate_new: float
        the new drop-out rate for output
        
    Returns
    -------
    drop_prob_new: array like
        the updated drop-out probability with the 
        average drop-out rate as rate_new
    
    """
    gaps_all = np.arange(-10, 10, 0.05)
    rate_all = np.zeros(len(gaps_all))
    
    drop_logit = logit(drop_prob)
    for i in range(len(gaps_all)):
        drop_prob_tmp = logistic(drop_logit + gaps_all[i])
        rate_all[i] = np.mean(drop_prob_tmp)
    
    idx = np.argmin(np.abs(rate_all-rate_new))
    drop_prob_new = logistic(drop_logit + gaps_all[idx])
    
    return drop_prob_new


def main():
    #part 0. parse command line options
    parser = OptionParser()
    parser.add_option("--anno_file", "-a", dest="anno_file", default=None,
        help="Annotation file for genes and transcripts.")
    parser.add_option("--ref_file", "-f", dest="ref_file", default=None,
        help="Reference genome in fasta formate.")
    parser.add_option("--out_dir", "-o", dest="out_dir", default=None, 
        help="Directory of the output files.")
    parser.add_option("--dice_file", "-d", dest="dice_file", default=None, 
        help="diceseq output file from bulk RNA-seq.")

    group = OptionGroup(parser, "Optional arguments")
    group.add_option("--tranLevel", action="store_true", dest="tran_level",
        default=False, help="Dropout at transcript level; otherwise gene level")
    group.add_option("--dropoutRate", "-r", dest="dropout_rate", type="float", 
        default=None, help="Dropout rate on average.")
    group.add_option("--dropoutProb", dest="dropout_prob", default=None,  
        help="Dropout probability of transcript. This will ignore the "
        "dropoutRate argument. File formate (tsv with header): gene,tran,prob.")
    group.add_option("--num-reads", "-N", dest="num_reads", type="int", 
        default=1000000, help="Number of reads in total. [default: %default]")
    parser.add_option_group(group)

    group = OptionGroup(parser, "Spanki arguments")
    group.add_option("-m", dest="mismatch_mode", default="random", 
        help="Error mode: random, errorfree, NIST, dm3, flyheads, or custom."
        " [default: %default]")
    group.add_option("--bp", dest="read_len", type="int", default=76, 
        help="Length of each read. [default: %default]")
    group.add_option("--frag", dest="frag_len", type="int", default=200, 
        help="Length of fragments. [default: %default]")
    group.add_option("--ends", dest="ends_num", type="int", default=2, 
        help="Number of reads ends: 1 or 2. [default: %default]")
    parser.add_option_group(group)
    
    ### under development
    # group.add_option("--corr-FPKM", "-c", dest="corr_FPKM", type="float", 
    #     default=0.7, help="Pearson's correlation coefficient between log2 FPKM "
    #     "and dropout probablity. [default: %default]")

    
    (options, args) = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        print("Welcome to dice-simulator for single-cell RNA-seq!\n")
        print("use -h or --help for help on argument.")
        sys.exit(1)
    if options.anno_file == None:
        print("[dice-simu] Error: need --anno_file for annotation.")
        sys.exit(1)
    else: 
        anno_file = options.anno_file
    if options.ref_file == None:
        print("[dice-simu] Error: need --ref_file for reference genome seq.")
        sys.exit(1)
    else:
        ref_file = options.ref_file
    if options.dice_file == None:
        print("[dice-simu] Error: need --dice_file for DICEseq output file.")
        sys.exit(1)
    else:
        dice_data = np.genfromtxt(options.dice_file, skip_header=1, dtype="str")
        tran_ids = dice_data[:,0]
        gene_ids = dice_data[:,1]
        tran_len = dice_data[:,3].astype(float)
        FPKM_all = dice_data[:,4].astype(float)
    if options.tran_level:
        flag_ids = tran_ids
    else:
        flag_ids = gene_ids

    num_reads = options.num_reads
    if options.dropout_prob is None:
        dropout_prob = np.ones(len(dice_data)) * 0.001
    else:
        temp = np.genfromtxt(options.dropout_prob, skip_header=1, dtype="str")
        idx = id_mapping(tran_ids, temp[:, 0])
        dropout_prob = temp[idx,2].astype(float)
        dropout_prob[dropout_prob<0.001] = 0.001
        dropout_prob[dropout_prob>0.999] = 0.999
    if options.dropout_rate is not None:
        idx_drop = FPKM_all > 0
        dropout_prob[idx_drop] = adjust_drop_prob(dropout_prob[idx_drop], 
            options.dropout_rate)

    if options.out_dir is None:
        out_dir = os.path.join(os.path.dirname(options.dice_file), "simuRNA")
    else:
        out_dir = options.out_dir
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    rpk_file = os.path.join(out_dir, "tran_rpk.txt")
    
    np.random.seed(0)
    flag = flag_ids[0]
    keep = np.random.binomial(1, 1-dropout_prob[0])
    FPKM = np.zeros(len(FPKM_all))
    for i in range(len(FPKM_all)):
        if flag != flag_ids[i]:
            flag = flag_ids[i]
            keep = np.random.binomial(1, 1-dropout_prob[i])
        FPKM[i] = keep * FPKM_all[i]
    rpk = FPKM * num_reads * 1000.0 / (np.sum(FPKM*tran_len))
    print("Drop-out rate: %.3f" %np.mean(rpk[idx_drop]==0))

    fid = open(rpk_file, "w")
    fid.writelines("txid\trpk\n")
    for i in range(len(tran_ids)):
        aLine = "%s\t%.4f\n" %(tran_ids[i], rpk[i])
        fid.writelines(aLine)
    fid.close()

    bashCommand = "spankisim_transcripts -o %s -g %s -f %s -t %s " %(out_dir, 
        anno_file, ref_file, rpk_file)
    bashCommand += "-bp %d -frag %d -ends %d -m %s" %(options.read_len, 
        options.frag_len, options.ends_num, options.mismatch_mode)
    print(bashCommand)
    pro = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output = pro.communicate()[0]

    bashCommand = "gzip %s/sim_1.fastq %s/sim_2.fastq" %(out_dir, out_dir) 
    pro = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output = pro.communicate()[0]

    bashCommand = "rm -rf %s/tmp %s/log %s/sim.*" %(out_dir, out_dir, out_dir) 
    pro = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output = pro.communicate()[0]


if __name__ == "__main__":
    main()

