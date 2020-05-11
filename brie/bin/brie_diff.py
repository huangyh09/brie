# This file is to compute Bayes factor for BRIE

import os
import sys
import time
import gzip
import numpy as np
import multiprocessing
from optparse import OptionParser, OptionGroup

FID = None
PROCESSED = 0
TOTAL_GENE = 0
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
        
        sys.stdout.write('\r[Brie-diff] [%s] %.1f%% done in %.1f sec.' 
            % (bar, percents, run_time))
        sys.stdout.flush()

        FID.writelines(RV)
    return RV


def get_prob(x1, x2, method="empirical"):
    """
    Calculate probability of the difference of two distributions.

    method: empirical, GDE

    """
    if method == "empirical":
        diff = x1 - x2
        prob = np.mean(np.abs(diff) <= 0.05)

    return prob

def get_BF(data, cell_names, rand_idx, minBF=0):
    """
    Calculate Bayes factor for an event in multiple cells

    Parameters
    ----------
    data: list of csv string
        each element of the list is a line of csv file for a cell.
    """
    RV_line = ""
    maxBF = rand_idx.shape[0] * 2
    tran_id = data[0][0].split(",")[0]
    gene_id = data[0][0].split(",")[1]
    for i in range(len(data)):
        data1 = data[i][0].split(",")
        c11 = round(float(data1[2])) #counts
        c12 = round(float(data[i][1])) #counts
        u1 = float(data1[3]) #prior Y
        s1 = float(data1[4]) #signma for Y
        x1 = np.array(data1[5:], float)[rand_idx[:,0]] #posterior samples
        y1 = np.random.normal(u1, s1, rand_idx.shape[0])

        for j in range(i+1, len(data)):
            data2 = data[j][0].split(",")
            c21 = round(float(data2[2]))
            c22 = round(float(data[j][1]))
            u2 = float(data2[3])
            s2 = float(data2[4])
            x2 = np.array(data2[5:], float)[rand_idx[:,1]]
            y2 = np.random.normal(u2, s2, rand_idx.shape[0])

            post_prob = get_prob(x1, x2)
            prior_prob = get_prob(logistic(y1), logistic(y2))
            if post_prob == 0:
                bf_val = maxBF
            else:
                bf_val = prior_prob / post_prob

            if bf_val < minBF:
                continue

            RV_line += "%s\t%s\t" %(tran_id, gene_id)
            RV_line += "%s\t%s\t" %(cell_names[i], cell_names[j])
            RV_line += "%.3f\t%.3f\t" %(logistic(u1), logistic(u2))
            RV_line += "%.3f\t%.3f\t" %(np.mean(x1), np.mean(x2))
            RV_line += "%d\t%d\t%d\t%d\t" %(c11, c12, c21, c22)
            RV_line += "%.1e\t%.1e\t%.1e\n" %(prior_prob, post_prob, bf_val)

    return RV_line

def count_BF(BF_file):
    """ Count the number of cell pairs with BF passing the threshold.
    """
    
    # data = np.genfromtxt(out_file+".tsv", dtype="str", delimiter="\t",
    #                      skip_header=1)
    # gene_unique, gene_counts = np.unique(data[:, 1], return_counts=True)
    
    gene_ids = []
    pair_BFs = []
    headers = 1
    with open(BF_file, "r") as f:
        for line in f:
            if headers > 0:
                headers = headers - 1
            else:
                line_val = line.rstrip().split("\t")
                gene_ids.append(line_val[1])
                pair_BFs.append(float(line_val[-1]))
    idx_sorted = np.argsort(gene_ids)
    gene_ids = np.array(gene_ids)[idx_sorted]
    pair_BFs = np.array(pair_BFs)[idx_sorted]
    
    gene_unique, gene_counts, mean_BF, median_BF = [], [], [], []
    for g in range(len(gene_ids)):
        if g == 0:
            last_g = g
            gene_tmp = gene_ids[g]
            continue
        if (g == len(gene_ids)-1) or (gene_ids[g] != gene_tmp):
            gene_unique.append(gene_ids[last_g])
            gene_counts.append(g - last_g)
            mean_BF.append(np.mean(pair_BFs[last_g:g]))
            median_BF.append(np.median(pair_BFs[last_g:g]))
            last_g = g
            gene_tmp = gene_ids[g]
            
    return gene_unique, gene_counts, mean_BF, median_BF


def main():
    # import warnings
    # warnings.filterwarnings('error')

    # parse command line options
    parser = OptionParser()

    parser.add_option("--inFiles", "-i", dest="in_files", default=None,
        help="Input files of Brie samples for multiple cells, comma separated "
        "for each cell, e.g., cell1,cell2,cell3")
    parser.add_option("--outFile", "-o", dest="out_file", default=None, 
        help="Output file with full path")

    group = OptionGroup(parser, "Optional arguments")
    group.add_option("--nproc", "-p", type="int", dest="nproc", default=4,
        help="Number of subprocesses [default: %default]") 
    group.add_option("--bootstrap", "-n", type="int", dest="bootstrap", 
        default=1000, help="Number of bootstrap [default: %default]")
    group.add_option("--minBF", dest="minBF", type="float", default=10, 
        help="Minimum BF for saving out, e.g., 3 or 10. If it is 0, save all "
        "events [default: %default]")
    parser.add_option_group(group)
    
    (options, args) = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        print("Welcome to Brie-diff!\n")
        print("use -h or --help for help on argument.")
        sys.exit(1)

    # Load input files
    if options.in_files is None:
        print("[Brie-diff] Error: need BRIE sample files.")
        sys.exit(1)
    else:
        cell_names = []
        samp_files = []
        samp_files_raw = options.in_files.split(",")
        for i in range(len(samp_files_raw)):
            _file = os.path.abspath(samp_files_raw[i])
            if os.path.isdir(_file):
                if os.path.isfile(os.path.join(_file, "samples.csv.gz")):
                    _file = os.path.join(_file, "samples.csv.gz")
            if os.path.basename(_file) != "samples.csv.gz":
                continue
            samp_files.append(_file)
            cell_names.append(os.path.basename(os.path.split(_file)[0]))

        if len(samp_files) < 2:
            print("[Brie-diff] Error: only %d sample file" %(len(samp_files)))
            sys.exit(1)
        else:
            print("[Brie-diff] detecting differential splicing for %d cells "
                  %(len(samp_files)))

    # Check output file
    global FID, TOTAL_GENE
    if options.out_file is None:
        out_file = os.path.dirname(samp_files[0]) + "/../brie_BF"
    else:
        if options.out_file.endswith(".tsv"):
            out_file = options.out_file[:-4]
        else:
            out_file = options.out_file
    try:
        FID = open(out_file + ".tsv", 'w')
    except IOError:
        sys.exit("[Brie-diff] Unable to write: " + out_file)
    items = ["tran_id", "gene_id", "cell1", "cell2", "prior1", "prior2", 
             "pis1", "psi2", "C1in", "C1out", "C2in", "C2out", "prior_prob", 
             "post_prob", "Bayes_factor"]
    FID.writelines("\t".join(items) + "\n")

    data = np.genfromtxt(samp_files[0], delimiter=",", dtype="str")
    tran_ids = data[:,0]
    samp_num = data.shape[1] - 5
    TOTAL_GENE = int(len(tran_ids) / 2)

    minBF = options.minBF
    rand_idx = np.random.randint(samp_num, size=(options.bootstrap, 2))

    # process many files.
    fid_in = []
    for i in range(len(samp_files)):
        fin = gzip.open(samp_files[i], 'rb')
        tmp = fin.readline()
        fid_in.append(fin)

    # get Bayes factor
    results = []
    if options.nproc <= 1:
        for k in range(TOTAL_GENE):
            line_all = []
            for i in range(len(fid_in)):
                line1 = fid_in[i].readline()
                line2 = fid_in[i].readline()
                line1 = line1.decode('utf8').strip()
                line2 = line2.decode('utf8').strip()
                line_all.append([line1, line2.split(",")[2]])

            result = get_BF(line_all, cell_names, rand_idx, minBF)
            show_progress(result)
    else:
        pool = multiprocessing.Pool(processes=options.nproc)
        for k in range(TOTAL_GENE):
            line_all = []
            for i in range(len(fid_in)):
                line1 = fid_in[i].readline()
                line2 = fid_in[i].readline()
                line1 = line1.decode('utf8').strip()
                line2 = line2.decode('utf8').strip()
                line_all.append([line1, line2.split(",")[2]])

            pool.apply_async(get_BF, (line_all, cell_names, rand_idx, minBF),
                callback=show_progress)
        pool.close()
        pool.join()
    
    FID.close()
    print("")
    print("[Brie-diff] Finished for %d splicing events." %(TOTAL_GENE))
    
    
    #rank genes by the number of cell pairs with differential splicing
    gene_unique, gene_counts, mean_BF, median_BF = count_BF(out_file + ".tsv")
                        
    idx_sorted = np.argsort(gene_counts)[::-1]
    fid = open(out_file + ".rank.tsv", "w")
    fid.writelines("gene_id\tcell_pairs\tmean_BF\tmedian_BF\n")
    for i in idx_sorted:
        fid.writelines("%s\t%d\t%.2f\t%.2f\n" %(gene_unique[i], gene_counts[i], 
                       mean_BF[i], median_BF[i]))
    fid.close()
    

if __name__ == "__main__":
    main()
    