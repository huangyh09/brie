# This file is to compute Bayes factor for Brie

import os
import sys
import time
# import h5py
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
    parser.add_option("--cell1_files", "-1", dest="cell1_files", default=None,
        help="Brie output file for cell (group) 1")
    parser.add_option("--cell2_files", "-2", dest="cell2_files", default=None,
        help="Brie output file for cell group 2")
    parser.add_option("--out_file", "-o", dest="out_file", default=None, 
        help="Output files with full path")
    parser.add_option("--bootstrap", "-n", type="int", dest="bootstrap", 
        default=1000, help="Number of bootstrap [default: %default]")
    parser.add_option("--maxBF", "-m", type="float", dest="maxBF", 
        default=100000, help="maximum Bayes factor [default: %default]")
    

    (options, args) = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        print("Welcome to Brie-diff!\n")
        print("use -h or --help for help on argument.")
        sys.exit(1)

    if options.out_file is None:
        out_file = os.path.dirname(os.path.abspath(cell1_file)) + "/brie_BF.tsv"
    else:
        out_file = options.out_file
    try:
        fid = open(out_file, 'w')
    except IOError:
        sys.exit("[Brie-diff] Unable to write: " + out_file)

    if options.cell1_files is None or options.cell2_files is None:
        print("[Brie-diff] Error: need file on both conditions.")
        sys.exit(1)
    else:
        print("[Brie-diff] detecting differential splicing from files:")
        print("    - %s" %options.cell1_files)
        print("    - %s" %options.cell2_files)

    #######
    cell1_files = options.cell1_files.split(",")
    cell2_files = options.cell2_files.split(",")

    data1 = np.genfromtxt(cell1_files[0], delimiter=",", dtype="str")
    data2 = np.genfromtxt(cell2_files[0], delimiter=",", dtype="str")
    idx = np.arange(0, data1.shape[0], 2)
    tran_ids1 = data1[idx, 0]
    tran_ids2 = data2[idx, 0]
    
    y1 = np.zeros((len(idx), len(cell1_files)))
    y2 = np.zeros((len(idx), len(cell2_files)))
    sigma1 = np.zeros(len(cell1_files))
    sigma2 = np.zeros(len(cell2_files))

    y1[:,0] = data1[idx, 3].astype(float) #prior Y
    y2[:,0] = data2[idx, 3].astype(float) #prior Y
    sigma1[0] = data1[idx, 4].astype(float).mean()
    sigma2[0] = data2[idx, 4].astype(float).mean()

    x1 = data1[idx, 5:].astype(float) #posterior samples
    x2 = data2[idx, 5:].astype(float) #posterior samples
    c11 = np.round(data1[idx, 2].astype(float))
    c21 = np.round(data2[idx, 2].astype(float))
    c12 = np.round(data1[idx+1, 2].astype(float))
    c22 = np.round(data2[idx+1, 2].astype(float))

    for i in range(1, len(cell1_files)):
        data1 = np.genfromtxt(cell1_files[i], delimiter=",", dtype="str")
        y1[:,i] = data1[idx, 3].astype(float)
        sigma1[i] = data1[idx, 4].astype(float).mean()
        x1 = np.append(x1, data1[idx, 5:].astype(float), axis=1)
        c11 += np.round(data1[idx, 2].astype(float))
        c12 += np.round(data1[idx+1, 2].astype(float))

    for i in range(1, len(cell2_files)):
        data2 = np.genfromtxt(cell2_files[i], delimiter=",", dtype="str")
        y2[:,i] = data2[idx, 3].astype(float)
        sigma2[i] = data2[idx, 4].astype(float).mean()
        x2 = np.append(x2, data2[idx, 5:].astype(float), axis=1)
        c21 += np.round(data2[idx, 2].astype(float))
        c22 += np.round(data2[idx+1, 2].astype(float))


    ########
    data = np.zeros((len(idx), 11))
    data[:,0] = logistic(y1).mean(axis=1)
    data[:,1] = logistic(y2).mean(axis=1)
    data[:,2] = x1.mean(axis=1)
    data[:,3] = x2.mean(axis=1)
    data[:,4] = c11
    data[:,5] = c12
    data[:,6] = c21
    data[:,7] = c22

    # Bayes factor #including GDE
    maxBF = options.maxBF
    bootstrap = options.bootstrap
    for i in range(len(x1)):
        a1, a2 = np.array([]), np.array([])
        for j in range(y1.shape[1]):
            a1 = np.append(a1, np.random.normal(y1[i,j], sigma1[j], 
                int(bootstrap/len(cell1_files))))
        for j in range(y2.shape[1]):
            a2 = np.append(a2, np.random.normal(y2[i,j], sigma2[j], 
                int(bootstrap/len(cell1_files))))
        prior_diff = np.random.permutation(logistic(a1)) - logistic(a2)
        data[i,8] = np.mean(np.abs(prior_diff) <= 0.05)
        
        idx1_perm = np.random.randint(x1.shape[1], size=bootstrap)
        idx2_perm = np.random.randint(x2.shape[1], size=bootstrap)
        post_diff = x1[i,idx1_perm] - x2[i,idx2_perm]
        data[i,9] = np.mean(np.abs(post_diff) <= 0.05)
        data[i,10] = data[i,8] / max(data[i,9], data[i,8]/maxBF)

    labels = ["prior1", "prior2", "pis1", "psi2", "C11", "C12", "C21", "C22", 
              "prior_prob", "post_prob", "Bayes_factor"]
    fid.writelines("tran_id\t" + "\t".join(labels) + "\n")
    for i in range(data.shape[0]):
        aline = "\t".join(["%.2f" %x for x in data[i,:]])
        # aline = "\t".join(["%.2f" %x for x in data[i,:-1]])
        # aline += "\t%.2e" %(data[i,-1])
        fid.writelines(tran_ids1[i] + "\t" + aline+"\n")
    fid.close()

    print("[Brie-diff] Finished for %d splicing events." %len(idx))
        

if __name__ == "__main__":
    main()
    