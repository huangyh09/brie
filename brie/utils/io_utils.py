# Containing API to load the count matrix data

import numpy as np
from scipy.sparse import csc_matrix

def load_brie_count(path):
    """
    Load dictionary-format sparse count matrix
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
        mat_dict[_key] = csc_matrix((_mat[:, 2], (_mat[:, 0], _mat[:, 1])), 
                                    shape=(n_gene, n_cell))
        
    return mat_dict


def load_gene_info():
    pass

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
