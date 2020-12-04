# functions for preprocessing

import numpy as np

def filter_genes(data, min_counts=0, min_cells=0, 
                 min_counts_uniq=0,  min_cells_uniq=0, min_MIF_uniq=0.001,
                 uniq_layers=['isoform1', 'isoform2'],
                 ambg_layers=['ambiguous'], copy=False):
    """Filter genes based on number of cells or counts.
    Keep genes that have at least `min_counts` counts or are expressed in at
    least `min_cells` cells at the total counts of uniquely counts. 
    
    Parameters
    ----------
    data : :class:`~anndata.AnnData`, `np.ndarray`, `sp.spmatrix`
        The (annotated) data matrix of shape `n_obs` × `n_vars`. Rows correspond
        to cells and columns to genes.
    min_counts : `int`, optional (default: 0)
        Minimum number of counts required for a gene to pass filtering.
    min_cells : `int`, optional (default: 0)
        Minimum number of cells expressed required for a gene to pass filtering.
    min_counts_uniq : `int`, optional (default: 0)
        Minimum number of unique counts required for a gene to pass filtering.
    min_cells_uniq : `int`, optional (default: 0)
        Minimum number of cells with unique reads required for a gene to pass filtering.
    uniq_layers : list, optional (default: ['isoform1', 'isoform 2'])
        Layers to aggregate for unique counts
    ambg_layers : list, optional (default: ['ambiguous'])
        Layers to aggregate with unique layers for total counts
    copy : `bool`, optional (default: `False`)
        Determines whether a copy is returned.
        
    Returns
    -------
    Filters the object and adds `n_counts` and `n_counts_uniq` to `adata.var`.
    """
    adata = data.copy() if copy else data
    
    unique_counts = np.zeros(adata.shape)
    for _layer in uniq_layers:
        unique_counts += adata.layers[_layer]
    
    total_counts = unique_counts.copy()
    for _layer in ambg_layers:
        total_counts += adata.layers[_layer]
    
    # filtering
    gene_subset = np.ones(adata.n_vars, dtype=bool)
    
    gene_subset &= np.array(total_counts.sum(0)).reshape(-1) >= min_counts
    gene_subset &= np.array((total_counts > 0).sum(0)).reshape(-1) >= min_cells
    
    gene_subset &= np.array(unique_counts.sum(0)).reshape(-1) >= min_counts_uniq
    gene_subset &= np.array((unique_counts > 0).sum(0)).reshape(-1) >= min_cells_uniq
    
    # limiting the isoform fractions
    unique_counts1 = adata.layers[uniq_layers[0]]
    unique_counts2 = adata.layers[uniq_layers[1]]
    gene_subset &= (np.array(unique_counts1.sum(0)).reshape(-1) >= 
                    min_MIF_uniq * np.array(unique_counts.sum(0)).reshape(-1))
    gene_subset &= (np.array(unique_counts2.sum(0)).reshape(-1) >= 
                    min_MIF_uniq * np.array(unique_counts.sum(0)).reshape(-1))
    
    adata._inplace_subset_var(gene_subset)
    adata.var['n_counts'] = np.array(total_counts.sum(0)).reshape(-1)[gene_subset]
    adata.var['n_counts_uniq'] = np.array(unique_counts.sum(0)).reshape(-1)[gene_subset]

    s = np.sum(~gene_subset)
    if s > 0:
        terms = []
        if min_cells > 0:
            terms.append('%d cells with any count' %(min_cells))
        if min_counts > 0:
            terms.append('%d total counts' %(min_counts))
        if min_cells_uniq > 0:
            terms.append('%d cells with unique counts' %(min_cells_uniq))
        if min_counts_uniq > 0:
            terms.append('%d unique counts' %(min_counts_uniq))  
        if min_MIF_uniq > 0:
            terms.append('%.4f minor isoform frequency' %(min_MIF_uniq))  
        print('Filtered out %d genes with less than ' %(s) + " or ".join(terms))

    return adata if copy else None
