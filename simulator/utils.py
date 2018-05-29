# utils file for simulation

import numpy as np

def id_mapping(IDs1, IDs2):
    """
    Mapping IDs2 to IDs1, both of which should only contain unique ids.
    
    Parameters
    ----------
    IDs1 : array_like or list
        ids for reference.
    IDs2 : array_like or list
        ids waiting to map.
        
    Returns
    -------
    RV_idx : array_like, the same length of IDs1
        The index for IDs2 mapped to IDs1. If an id in IDs1 does not exist 
        in IDs2, then return a None for that id.
    """
    idx1 = np.argsort(IDs1)
    idx2 = np.argsort(IDs2)
    RV_idx1, RV_idx2 = [], []
    
    i, j = 0, 0
    while i < len(idx1):
        if j == len(idx2) or IDs1[idx1[i]] < IDs2[idx2[j]]:
            RV_idx1.append(idx1[i])
            RV_idx2.append(None)
            i += 1
        elif IDs1[idx1[i]] == IDs2[idx2[j]]:
            RV_idx1.append(idx1[i])
            RV_idx2.append(idx2[j])
            i += 1
            j += 1
        elif IDs1[idx1[i]] > IDs2[idx2[j]]:
            j += 1
            
    origin_idx = np.argsort(RV_idx1)
    RV_idx = np.array(RV_idx2)[origin_idx]
    return RV_idx


def get_fraction(gene_ids, FPKM, ignoreNan=False):
    """Get the fraction from FPKM"""
    idx0 = 0
    frac = np.zeros(len(FPKM))
    for i in range(len(FPKM)+1):
        if i >= len(FPKM) or gene_ids[idx0] != gene_ids[i]:
            FPKM_sum = float(np.sum(FPKM[idx0:i]))
            if FPKM_sum == 0:
                if ignoreNan == True:
                    frac[idx0:i] = None
                else:
                    frac[idx0:i] = 1.0/(i-idx0)
            else:
                frac[idx0:i] = FPKM[idx0:i] / FPKM_sum
            idx0 = i
    return frac



def loadSpanki(file_name, tran_ids, gene_ids=None,
    ignoreNan=False):
    frac = np.zeros(len(tran_ids))
    FPKM = np.zeros(len(tran_ids))
    FPKM[:] = None

    data = np.genfromtxt(file_name, skip_header=1, dtype="S100")
    idxT = id_mapping(tran_ids, data[:,1])
    idx_no_None = idxT >= 0
    FPKM[idx_no_None] = data[idxT[idx_no_None].astype(int), 3].astype(float)
    frac = get_fraction(gene_ids, FPKM, ignoreNan)
    
    return frac, FPKM





