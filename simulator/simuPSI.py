# This file to simulate single cell RNA-seq reads based on real bulk RNA-seq 
# expression profile and input dropout rate and number of reads.

import os
import sys
import time
import subprocess
import numpy as np
import scipy.stats as st
from optparse import OptionParser, OptionGroup

# import pyximport; pyximport.install()

from brie import loadgene
from utils import id_mapping, loadSpanki

START_TIME = time.time()

def logistic(x):
    return np.exp(x)/(1+np.exp(x))
    
def logit(x, minval=0.001):
    if isinstance(x,  (list, tuple, np.ndarray)):
        x[1-x<minval] = 1-minval
        x[x<minval] = minval
    else:
        x = max(minval, x)
        x = min(1-minval, x)
    val = np.log(x/(1-x))
    return val

def generate_prior(psi, corr=0.8, 
        min_sigma=0.1, max_sigma=5, steps=2000):
    """
    Generating prior with correlated to original psi.
    
    Parameters
    ----------
    psi: array like
        PSI values, each element ranges in [0,1]
    corr: float
        Pearson's correlation between psi and prior
        
    Reterns
    -------
    prior: array like
        Generated prior
    """
    sigma_all = np.linspace(min_sigma, max_sigma, steps)
    corr_all  = np.zeros(len(sigma_all))
    psi_logit = logit(psi, minval=0.0001)
    for i in range(len(sigma_all)):
        _prior_logit = psi_logit + np.random.normal(0, sigma_all[i], size=len(psi))
        corr_all[i]  = st.pearsonr(logistic(_prior_logit), psi)[0]
    idx = np.argmin(np.abs(corr_all-corr))
    prior_logit = psi_logit + np.random.normal(0, sigma_all[idx], size=len(psi))
    return logistic(prior_logit)


def main():
    #part 0. parse command line options
    parser = OptionParser()
    parser.add_option("--anno_file", "-a", dest="anno_file", default=None,
        help="Annotation file for genes and transcripts.")
    parser.add_option("--ref_file", "-f", dest="ref_file", default=None,
        help="Reference genome in fasta formate.")
    parser.add_option("--out_dir", "-o", dest="out_dir", default=None, 
        help="Directory of the output files.")
    parser.add_option("--rpk", dest="rpk", type="float", default=100, 
        help="RPK value. [default: %default]")
    parser.add_option("--mode", dest="mode", default="LogitNormal", 
        help="Distribution mode: LogitNormal, UniDiff1, UniDiff2. "
        " [default: %default]")
    parser.add_option("--theta", dest="theta", type="float", default=3.0, 
        help="std for logit normal distribution. [default: %default]")
    parser.add_option("--priorR", dest="priorR", type="float", default=0.8, 
        help="Peason's R for generating priori. [default: %default]")
    parser.add_option("--seed", dest="seed", type="float", default=None, 
        help="Seed for random number. [default: non-seed]")

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
        genes = loadgene(anno_file)
        gene_ids, tran_ids = [], []
        for g in genes:
            for t in g.trans:
                tran_ids.append(t.tranID)
                gene_ids.append(g.geneID)

    if options.ref_file == None:
        print("[dice-simu] Error: need --ref_file for reference genome seq.")
        sys.exit(1)
    else:
        ref_file = options.ref_file

    if options.out_dir is None:
        out_dir = os.path.join(os.path.dirname(options.dice_file), "simuRNA")
    else:
        out_dir = options.out_dir
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    rpk_file = os.path.join(out_dir, "tran_rpk.txt")

    rpk = options.rpk
    if options.seed is not None:
        np.random.seed(int(options.seed))
    if options.mode == "LogitNormal":
        psi = logistic(np.random.normal(0, options.theta, size=len(tran_ids)/2))
    elif options.mode == "UniDiff1":
        x = [0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9]
        psi = np.tile(x, np.ceil(len(tran_ids)/2.0/len(x)))[:len(tran_ids)/2]
    elif options.mode == "UniDiff2":
        x = [0.9, 0.8, 0.7, 0.6, 0.4, 0.3, 0.2, 0.1]
        psi = np.tile(x, np.ceil(len(tran_ids)/2.0/len(x)))[:len(tran_ids)/2]

    #Under study for these modes
    elif options.mode == "Uniform":
        psi = np.ones(len(tran_ids)/2.0) * 0.5
    elif options.mode == "Diff1":
        psi = logistic(np.random.normal(0, options.theta, size=len(tran_ids)/2))
        ### 20% differential splicing
        x = [0.05, 0.2, 0.35, 0.65, 0.8, 0.95]
        diff_num = int(0.3 * len(tran_ids)/2)
        psi[:diff_num] = np.tile(x, np.ceil(diff_num/len(x)))[:diff_num]
    elif options.mode == "Diff2":
        psi = logistic(np.random.normal(0, options.theta, size=len(tran_ids)/2))
        x = [0.95, 0.8, 0.65, 0.35, 0.2, 0.05]
        diff_num = int(0.3 * len(tran_ids)/2)
        psi[:diff_num] = np.tile(x, np.ceil(diff_num/len(x)))[:diff_num]

    fid = open(rpk_file, "w")
    fid.writelines("txid\trpk\n")
    for i in range(len(psi)):
        fid.writelines("%s\t%.4f\n" %(tran_ids[i*2], rpk*psi[i]))
        fid.writelines("%s\t%.4f\n" %(tran_ids[i*2+1], rpk*(1-psi[i])))
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

    ## generate prior
    simu_truth = loadSpanki(out_dir+"/transcript_sims.txt", np.array(tran_ids), 
        np.array(gene_ids))[0][range(0, len(tran_ids), 2)]
    prior = generate_prior(simu_truth, corr=options.priorR, 
        min_sigma=0.1, max_sigma=5, steps=2000)

    fid = open(out_dir + "/prior_fractions.txt", "w")
    fid.writelines("gene_id,feature1\n")
    for i in range(len(prior)):
        fid.writelines("%s,%.3f\n" %(gene_ids[i*2], logit(prior[i])))
        # fid.writelines("%s\t%s\t%.3f\n" %(tran_ids[i*2+1], gene_ids[i*2+1], 
        # 1-prior[i]))
    fid.close()


if __name__ == "__main__":
    main()

