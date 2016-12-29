# This file is to compute Bayes factor for Brie

import os
import sys
import time
import h5py
import numpy as np
import multiprocessing
from optparse import OptionParser, OptionGroup

PROCESSED = 0
TOTAL_GENE = 0
TOTAL_READ = []
START_TIME = time.time()

def logistic(x):
    return np.exp(x) / (np.exp(x) + 1)
    
def show_progress(RV=None):
    global PROCESSED, TOTAL_GENE, START_TIME
    if RV is not None: 
        PROCESSED += 1
        bar_len = 20
        run_time = time.time() - START_TIME
        percents = 100.0 * PROCESSED / TOTAL_GENE
        filled_len = int(bar_len * percents / 100)
        bar = '=' * filled_len + '-' * (bar_len - filled_len)
        
        sys.stdout.write('\r[Brie] [%s] %.1f%% done in %.1f sec.' 
            % (bar, percents, run_time))
        sys.stdout.flush()
    return RV

def main():
    # import warnings
    # warnings.filterwarnings('error')

    # parse command line options
    parser = OptionParser()
    parser.add_option("--cond1_file", "-1", dest="cond1_file", default=None,
        help="Brie output file for condition 1")
    parser.add_option("--cond2_file", "-2", dest="cond2_file", default=None,
        help="Brie output file for condition 2")
    parser.add_option("--bootstrap", "-n", type="int", dest="bootstrap", 
        default=1000, help="Number of bootstrap [default: %default]")
    parser.add_option("--maxBF", "-m", type="float", dest="maxBF", 
        default=100000, help="maximum Bayes factor [default: %default]")
    parser.add_option("--out_file", "-o", dest="out_file", default=None, 
        help="Output files with full path")

    (options, args) = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        print("Welcome to Brie-diff!\n")
        print("use -h or --help for help on argument.")
        sys.exit(1)
    if options.cond1_file is None or options.cond2_file is None:
        print("[Brie-diff] Error: need file on both conditions.")
        sys.exit(1)
    else:
        print options.cond1_file
        print options.cond2_file
        f = h5py.File(options.cond1_file, "r")
        W1 = np.array(f["W_sample"]).mean(axis=1)
        sigma1 = np.array(f["sigma"])[0]
        counts1 = np.array(f["counts"])
        Psi1_all = np.array(f["Psi_sample"])
        features1 = np.array(f["features"])
        tran_ids1 = np.array(f["tran_ids"])
        f.close()

        f = h5py.File(options.cond2_file, "r")
        W2 = np.array(f["W_sample"]).mean(axis=1)
        sigma2 = np.array(f["sigma"])[0]
        counts2 = np.array(f["counts"])
        Psi2_all = np.array(f["Psi_sample"])
        features2 = np.array(f["features"])
        tran_ids2 = np.array(f["tran_ids"])
        f.close()

    maxBF = options.maxBF
    bootstrap = options.bootstrap

    if options.out_file is None:
        out_file = os.path.dirname(os.path.abspath(cond1_file)) + "/brie_BF.tsv"
    else:
        out_file = options.out_file
    try:
        fid = open(out_file, 'w')
    except IOError:
        sys.exit("[Brie-diff] Unable to write: " + out_file)

    global TOTAL_TRAN
    TOTAL_TRAN = len(tran_ids1) / 2

    idx = np.arange(0, len(tran_ids1), 2)
    x1 = Psi1_all[idx,:]
    x2 = Psi2_all[idx,:]
    y1 = np.dot(features1[idx,:], W1)
    y2 = np.dot(features2[idx,:], W2)

    c11 = np.round(counts1[idx])
    c12 = np.round(counts1[idx+1])
    c21 = np.round(counts2[idx])
    c22 = np.round(counts2[idx+1])

    data = np.zeros((len(idx), 11))
    data[:,0] = logistic(y1)
    data[:,1] = logistic(y2)
    data[:,2] = x1.mean(axis=1)
    data[:,3] = x2.mean(axis=1)
    data[:,4] = c11
    data[:,5] = c12
    data[:,6] = c21
    data[:,7] = c22

    # Bayes factor
    for i in range(len(x1)):
        a1 = np.random.normal(y1[i], sigma1, bootstrap)
        a2 = np.random.normal(y2[i], sigma2, bootstrap)
        prior_diff = np.random.permutation(logistic(a1)) - logistic(a2)
        data[i,8] = np.mean(np.abs(prior_diff) <= 0.05)
        
        idx1_perm = np.random.randint(x1.shape[1], size=bootstrap)
        idx2_perm = np.random.randint(x2.shape[1], size=bootstrap)
        post_diff = x1[i,idx1_perm] - x2[i,idx2_perm]
        data[i,9] = np.mean(np.abs(post_diff) <= 0.05)
        data[i,9] = max(data[i,9], data[i,8] / maxBF)
        data[i,10] = data[i,8] / data[i,9]

    labels = ["prior1", "prior2", "pis1", "psi2", "C11", "C12", "C21", "C22", 
              "prior_prob", "post_prob", "Bayes_factor"]
    fid.writelines("tran_id\t" + "\t".join(labels) + "\n")
    for i in range(data.shape[0]):
        aline = "\t".join(["%.2f" %x for x in data[i,:]])
        fid.writelines(tran_ids1[i*2] + "\t" + aline+"\n")
    fid.close()
        

if __name__ == "__main__":
    main()
    