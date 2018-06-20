"""This script simulates data for pseudotime-enhanced BRIE

In all data structures, cells'ids occupy the first column, and genes'ids the
first raw.
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
from simuPSI import logistic, logit # simulation script for simple BRIE

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
    
    parser.add_argument('-g', '--nb-genes', type=int, default=20,
                        help='number of genes for data to simulate.')
    
    parser.add_argument('-t', '--theta', type=float, default=3.0,
                        help='standard deviation for simulated WX~N(0,theta).')
    
    parser.add_argument('-s', '--std-alpha', type=float, default=1.5,
                        help='standard deviation for simulated alpha~N(0,s).')
    
    parser.add_argument('-r', '--pearson-coefficient', type=float, default=0.8,
                        help="target Pearson's correlation coefficient between "
                        + "'true psi' and 'observed psi' simulated values.")

    parser.add_argument('--rpk', type=float, default=100,
                        help='reads per kilo-base.')

    return parser.parse_args()

## useless because a cell is identified by its index
# def generate_cell_ids(nb_cell):
#     """Generate a list of nb_cell cell ids

#     Parameters
#     ----------
#     nb_cell: int
#         number of cells.
        
#     Returns
#     -------
#     cell_ids: list
#         list of cell ids (strings).
#     """
#     cell_ids = ["cell" + str(i) for i in range(nb_cell)]
#    return cell_ids

def generate_pseudotime(cell_nb):
    """Generate pseudotimes for cell_nb cells

    Pseudotimes are generated with a uniform distribution between -0.5 and 0.5

    Parameters
    ----------
    cell_ids: list
        list of cells'ids (strings).
        
    Returns
    -------
    pseudotimes: list
        list of pseudotimes (float between -0.5 and 0.5)
    """
    pseudotimes = [np.random.uniform(-0.5,0.5) for i in range(cell_nb)]
    return pseudotimes

def generate_logit_psi(cell_nb, gene_nb, pseudotimes, theta, std_alpha):
    """TODO: fill this and look closer at logit...
    """
    psi = np.zeros((cell_nb, gene_nb))
    for g in gene_nb: # for each gene
        pseudotime_contribution = np.random.normal(0,std_alpha) * pseudotimes[g]
        for c in cell_nb: # for each cell
            psi[c,g] = np.random.normal(0,theta) + pseudotime_contribution
    return psi
    

def generate_prior(psi, theta, std_alpha,
                   corr=0.8, min_sigma=0.1, max_sigma=5, steps=2000):
    """Generating prior with corr correlation to original psi.
    
    Parameters
    ----------
    psi: array like
        PSI values, each element ranges in [0,1]
    corr: float
        Pearson's correlation between psi and prior
        
    Returns
    -------
    prior: array like
        Generated prior
    """
    return
    

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

    args = parse_arguments()
    
    cell_ids = generate_cell_ids(args.nb_cells)
    
    pseudotimes = generate_pseudotime(cell_ids)

if __name__ == "__main__":
    main()
