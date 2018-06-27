"""This script simulates data for pseudotime-enhanced BRIE

In all data structures, cells'ids occupy the first column, and genes'ids the
first raw. This script only support the two isoforms per gene scenario.
To see how to use this program, run `python3 pseudotimeSimuPSI.py -h`
BRIE paper can be consulted here:
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1248-5

ex:
annotation_file :
/home/milan/prog/cours/internship/brie-examples/test/gencode.vM12.annotation.gtf

#######
testing
cd /home/milan/prog/cours/internship/brie/simulator
python3 pseudotimeSimuPSI.py /home/milan/prog/cours/internship/pseudotime/notebook/matrix_of_counts.csv feature_matrix_dir destination_dir
"""

import os # to navigate and create directories
import argparse # to parse script arguments (since optparse is deprecated)
import numpy as np # required for pseudotime analysis
import pandas as pd # required for pseudotime analysis
import scanpy.api as sc # to perform pseudotime analysis
import scipy.stats as st # to compute Pearson's r correlation coefficient
from simuPSI import logistic, logit, generate_prior # simulation for simple BRIE
import subprocess

from diceseq import loadgene # to read and load data from annotation file
from diceseq.utils.out_utils import id_mapping
from diceseq.utils.misc_utils import loadresult

def parse_arguments():
    """ parse arguments of this script

    Returns
    -------
    args: Namespace
        Namespace of arguments (object with key-value pairs).
    """
    
    # default will be displayed in help thanks to formater_class:
    parser = argparse.ArgumentParser(description='This script simulates data'
                                     + ' for pseudotime-enhanced BRIE',
                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # mandatory arguments
    parser.add_argument('output_dir',
                        help='output directory.')
    
    parser.add_argument('annotation_file',
                        help='annotation gtf or gff3 file '
                        '(where genes and transcripts to consider are stored).')

    parser.add_argument('reference_genome_file',
                        help='reference genome in fasta format.')

    
    # optionnal arguments
    parser.add_argument('-c', '--nb-cells', type=int, default=30,
                        help='number of cells for data to simulate.')

    parser.add_argument('-g', '--max-gene-nb', type=int, default=None,
                        help='maximum number of genes for data to simulate'
                        + 'from annotaion file.')
    
    parser.add_argument('-t', '--theta', type=float, default=3.0,
                        help='standard deviation for simulated WX~N(0,theta).')
    
    parser.add_argument('-s', '--std-alpha', type=float, default=1.5,
                        help='standard deviation for simulated alpha~N(0,s).')
    
    parser.add_argument('-r', '--pearson-coefficient', type=float, default=0.8,
                        help="target Pearson's correlation coefficient between "
                        + "'true psi' and 'observed psi' simulated values.")

    parser.add_argument('--rpk', type=float, default=100,
                        help='reads per kilo-base.')

    parser.add_argument("--seed", type=int, default=None, 
                        help="Seed for random number.")

    # optional arguments for Spanki:
    group = parser.add_argument_group(title='Spanki arguments',
                                      description=None)

    group.add_argument("-m", "--mismatch_mode", default="random", 
        help="Error mode: random, errorfree, NIST, dm3, flyheads, or custom.")
    group.add_argument("--bp", dest="read_len", type=int, default=76, 
        help="Length of each read.")
    group.add_argument("--frag", dest="frag_len", type=int, default=200, 
        help="Length of fragments.")
    group.add_argument("--ends", dest="ends_num", type=int, default=2, 
        help="Number of reads ends: 1 or 2.")

    return parser.parse_args()


def generate_pseudotime(cell_nb):
    """Generate pseudotimes for cell_nb cells

    Pseudotimes are generated with a uniform distribution between -0.5 and 0.5

    Parameters
    ----------
    cell_nb: int
        number of cells.
        
    Returns
    -------
    pseudotimes: list
        list of pseudotimes (float between -0.5 and 0.5) for each cell.
    """
    pseudotimes = [np.random.uniform(-0.5,0.5) for i in range(cell_nb)]
    return pseudotimes

def generate_psi(cell_nb, gene_nb, pseudotimes, theta, std_alpha):
    """Generate 'true' psi value (exon inclusion ratio == isoform proportion).

    In our model, logit(psi) = WX + alpha.t
    where (in the simulation):
        WX ~ N(0, theta²)
        alpha ~ N(0, std_alpha²)
        t represent the pseudotime and belongs to (-0.5, 0.5)

    Parameters
    ----------
    cell_nb: int
        number of cells.
    gene_nb: int
        number of genes.
    pseudotimes: list
        list of pseudotimes (float between -0.5 and 0.5).
    theta: float
        standard deviation for WX contribution.
    std_alpha: float
        standard deviation for pseudotime contribution.
        
    Returns
    -------
    psi: numpy.array
        array of size (cell_nb, gene_nb) which represent exon inclusion ratio
        for each gene for each cell.
    """
    logit_psi = np.zeros((cell_nb, gene_nb))
    for g in range(gene_nb): # for each gene, create pseudotime contribution
        # ~N(0, std_alpha)
        pseudotime_contribution = np.random.normal(0, std_alpha) * pseudotimes[g]
        for c in range(cell_nb):#for each cell, add non pseudotime related random noise
            logit_psi[c,g] = np.random.normal(0,theta) + pseudotime_contribution
            
    psi = logistic(logit_psi)
    return psi
    

def generate_matrix_prior(psi, corr=0.8, 
                   min_sigma=0.1, max_sigma=5, steps=2000):
    """Generating prior with corr correlation to original psi.

    Finaly useless, generate prior from original SimuPSI.py works.
    
    Parameters
    ----------
    psi: array like
        PSI values, each element ranges in [0,1]
    corr: float
        Pearson's correlation between psi and prior
        
    Returns
    -------
    prior: numpy.array
        Generated prior. If psi was a 2D array, prior is a column vector array.
    """
    #discretised interval of sigmas from min_sigma to max_sigma with steps steps
    sigma_all = np.linspace(min_sigma, max_sigma, steps)#[0.1,..(2000steps)..,5.]
    
    corr_all  = np.zeros(len(sigma_all)) # Pearson's correlation coefficients r

    psi = np.array(psi) #provide compatibility of code for float, list and array

    # reshape 2D array psi to further compute pearson's r:
    psi_is_reshaped = False
    if psi.ndim == 2: # if psi is a two dimensionnal array
        first_dimension, second_dimension = psi.shape
        # reshape psi as a column vector (necessary to compute pearson's r)
        psi = psi.reshape((first_dimension * second_dimension, 1))
        psi_is_reshaped = True
        
    psi_logit = logit(psi, minval=0.0001)
    
    for i in range(len(sigma_all)): # for each value of sigma in sigma_all
        # generate prior with a sigma equal to sigma_all[i]:
        _prior_logit = psi_logit + np.random.normal(0, sigma_all[i],
                                                    size=psi.shape)
        corr_all[i] = st.pearsonr(logistic(_prior_logit), psi)[0] # r at round i
        
    idx = np.argmin(np.abs(corr_all-corr)) # index of r closest to target (corr)
    
    # generate prior with a sigma that gives a prior with target cor. coef. r:
    prior_logit = psi_logit + np.random.normal(0, sigma_all[idx], size=psi.shape)

    # come back to original shape:
    if psi_is_reshaped == True: # if psi has been reshaped
        prior_logit = prior_logit.reshape((first_dimension, second_dimension))
    
    return logistic(prior_logit)

    

def pseudotime(matrix_of_counts):
    """Compute pseudotimes of each cell of file of name matrix_of_counts

    Pseudotime analysis is computed through a diffusion map.

    Args:
    matrix_of_counts (string): name of the matrix of counts file

    Returns:
    pseudotimes (numpy.array): array where first column store cells'ids and
        second column stores corresponding pseudotimes (floats between 0 and 1).
    """
    adata = sc.read(matrix_of_counts) # read matrix of counts
    adata.uns['iroot'] = 1 # choose an arbitrary root cell
    # compute neighbors to prepare for diffusion map:
    sc.pp.neighbors(adata, n_neighbors=5, method='gauss', knn=False, use_rep='X')
    # compute branching and diffusion pseudotime:
    sc.tl.dpt(adata, n_branchings=1)
    # extract pseudotimes (two first columns):
    pseudotimes = adata.obs#[:,:2]
    
    return pseudotimes



def main():
    """Simulate reads for multiple cells with a pseudotime conponent

    Read arguments from sys.argv (thanks to argparse). Store cell specific data
    in a different directory for each cell (called cell_i with i the cell
    number) and shared data at the root of output directory.
    """

    args = parse_arguments()

    # extract genes and isoforms transcrits from annotation file
    anno_file = args.annotation_file
    genes = loadgene(anno_file) # load gene data from annotation file
    # gene_ids[i] store the gene id of transcript given by tran_ids[i]
    # hence, in 2 isoforms per gene scenario, genes[i].geneID==gene_ids[2*i]

    # control number of genes:
    if args.max_gene_nb is not None:
        gene_nb = min(args.max_gene_nb, len(genes))
    else:
        gene_nb = len(genes) # numbers of genes
    
    gene_ids, tran_ids = [], []
    for g in genes[:gene_nb]: # for each gene among the gene_nb first ones
        #tr_nb = 0#
        for t in g.trans:
            tran_ids.append(t.tranID)
            gene_ids.append(g.geneID)
            #tr_nb += 1#
        #print("nb of transcripts for this gene: %s" %(tr_nb))#

    cell_nb = args.nb_cells # number of cells to simulate

    ref_file = args.reference_genome_file

    out_dir = args.output_dir
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    rpk = args.rpk # reads per kilobases

    # seed for random numbers
    if args.seed is not None:
        np.random.seed(args.seed)

    pseudotimes = generate_pseudotime(cell_nb)

    psi_matrix = generate_psi(cell_nb, gene_nb, pseudotimes, args.theta,
                              args.std_alpha)

    for i in range(cell_nb): # for each simulated cell
        psi = psi_matrix[i,:] # psi values for cell i

        cell_dir = os.path.join(out_dir,"cell_" + str(i)) # directory for cell i
        if not os.path.exists(cell_dir):
            os.makedirs(cell_dir)

        rpk_file = os.path.join(cell_dir, "tran_rpk.txt")

        #########################################################
        #print("psi: %s\ntran_ids: %s" %(len(psi), len(tran_ids)))

        with open(rpk_file, "w") as fid:
            fid.writelines("txid\trpk\n")
            for i in range(len(psi)):
                fid.writelines("%s\t%.4f\n" %(tran_ids[i*2], rpk*psi[i]))
                fid.writelines("%s\t%.4f\n" %(tran_ids[i*2+1], rpk*(1-psi[i])))

        # generation of transcripts:
        bashCommand = "spankisim_transcripts -o %s -g %s -f %s -t %s" %(cell_dir,
                                                   anno_file, ref_file, rpk_file)
        bashCommand += " -bp %d -frag %d -ends %d -m %s" %(args.read_len,
                  args.frag_len, args.ends_num, args.mismatch_mode) 
        print(bashCommand)
        pro = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        output = pro.communicate()[0]

        bashCommand = "gzip %s/sim_1.fastq %s/sim_2.fastq" %(cell_dir, cell_dir) 
        pro = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        output = pro.communicate()[0]

        bashCommand = "rm -rf %s/tmp %s/log %s/sim.*" %(cell_dir, cell_dir, cell_dir) 
        pro = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        output = pro.communicate()[0]

        ## generate prior
    
        simu_truth = loadresult(cell_dir+"/transcript_sims.txt",
                                np.array(tran_ids), np.array(gene_ids),
                                method="spanki")[0][range(0, len(tran_ids), 2)]
        # loadresult(...)[0] is a numpy.array of fractions of reads per isoform
        # here we consider one over two isoforms, since we have only 2 isoforms

        # compute prior that is likely to have lead to simu_truth reads in our model
        prior = generate_prior(simu_truth, corr=args.pearson_coefficient, 
                               min_sigma=0.1, max_sigma=5, steps=2000)

        with open(cell_dir + "/prior_fractions.txt", "w") as fid:
            fid.writelines("gene_id,feature1\n")
            for i in range(len(prior)):
                fid.writelines("%s,%.3f\n" %(gene_ids[i*2], logit(prior[i])))
                # ith gene id is gene_ids[2*i] because
                # genes[i].geneID==gene_ids[2*i] in 2 isoforms per gene scenario
    
if __name__ == "__main__":
    main()
