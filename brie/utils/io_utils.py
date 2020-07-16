# Containing API to load the count matrix data

import anndata
import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix

from anndata import read_h5ad
from .gtf_utils import load_genes as read_gff


def convert_to_annData(Rmat_dict, effLen_tensor, cell_note, gene_note):
    """Convert matrices and annotation to annData
    """
    Rmat = {}
    for _key in Rmat_dict:
        Rmat[_key] = Rmat_dict[_key].astype(np.float32)#.toarray()
    Rmat.keys()
    
    X = Rmat['1'] + Rmat['2'] + Rmat['3']
    layers = {}
    layers['isoform1']  = Rmat['1']
    layers['isoform2']  = Rmat['2']
    layers['ambiguous'] = Rmat['3']
    layers['poorQual']  = Rmat['0']
    
    obs = pd.DataFrame(cell_note[1:, :],
                       index = cell_note[1:, 0],
                       columns = cell_note[0, :])
    
    var = pd.DataFrame(gene_note[1:, :],
                       index = gene_note[1:, 0],
                       columns = gene_note[0, :])
    
    Prob_tensor = effLen_tensor / effLen_tensor.sum(2, keepdims=True)
    
    varm = {}
    varm['effLen'] = np.append(effLen_tensor[:, 0, :],
                               effLen_tensor[:, 1, :], axis=1)
    varm['p_ambiguous'] = Prob_tensor[:, :, 2]
    
    adata = anndata.AnnData(X=X, obs=obs, var=var, varm=varm,
                            layers=layers, dtype='float32')
    return adata


def read_npz(path):
    """Read count data in the npz format into anaData
    """
    brie_dat = np.load(path, allow_pickle=True)
    cell_note = brie_dat['cell_note']
    gene_note = brie_dat['gene_note']
    Rmat_dict = brie_dat['Rmat_dict'].item()
    effLen_tensor = brie_dat['effLen_tensor']
    
    adata = convert_to_annData(Rmat_dict, effLen_tensor, cell_note, gene_note)
    return adata


def read_brieMM(path):
    """Read brie count generated Market martrix: dictionary-format 
    sparse count matrix
    """
    fid = open(path, 'r')
    lines = fid.readlines()
    fid.close()
    
    # check mtx file format
    n_gene, n_cell, size = lines[1].strip().split("\t")
    n_gene, n_cell, size = int(n_gene), int(n_cell), int(size)

    dat_dict = {}
    for _line in lines[2:]:
        i, j, _str = _line.strip().split("\t")
        _dat = eval(_str)
        for _key in _dat:
            if _key not in dat_dict:
                dat_dict[_key] = []
            dat_dict[_key].append([i, j, _dat[_key]])
        
    mat_dict = {}
    for _key in dat_dict:
        _mat = np.array(dat_dict[_key], dtype='int')
        _mat[:, :2] -= 1 # 0-based index
        mat_dict[_key] = csc_matrix(
            (_mat[:, 2], (_mat[:, 0], _mat[:, 1])), 
            shape=(n_gene, n_cell)
        )
        
    return mat_dict


def fetch_gene_info(genes, fraglen=None, out_file=None):
    """
    Extract the isoform information from a list of Gene
    """
    out_all = []
    for g in genes:
        tran_ids, tran_lens = [], []
        for t in g.trans:
            tran_ids.append(t.tranID)
            tran_lens.append(str(t.tranL))
        out_list = [g.geneID, g.geneName, ",".join(tran_lens), 
                    ",".join(tran_ids)]
        out_all.append(out_list)

    if out_file is not None:
        fid = open(out_dir + "/gene_note.tsv", "w")
        fid.writelines("GeneID\tGeneName\tTranLens\tTranIDs\n")
        for _line_val in out_all:
            fid.writelines("\t".join(_line_val) + "\n")
        fid.close()

    return out_all
