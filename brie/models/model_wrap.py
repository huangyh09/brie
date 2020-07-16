# Tools for BRIE package
# Wrap functions for running BRIE2

import numpy as np
import concurrent
import multiprocessing
from .model_TFProb import BRIE2
from ..settings import verbosity

from scipy.stats import chi2
from scipy.sparse import csc_matrix
from statsmodels.stats.multitest import multipletests
        

def fit_BRIE_matrix(data, Xc=None, Xg=None, effLen=None, intercept=None, 
                    intercept_mode='gene', LRT_index=None, **keyargs):
    """Fit a BRIE model with cell and/or gene features
    
    Parameters
    ----------
    data : :class:`np.ndarray`, `sp.spmatrix`
        List of matrices (2 or 3) for isoform1, isoform2, [ambiguous].
        Rows correspond to cells and columns to genes.
    Xc : `numpy.array`, (n_feature, n_cell), np.float32, optional
        Cell features.
    Xg : `numpy.array`, (n_gene, n_feature), np.float32, optional
        Gene features.
    **keyargs : keyargs for BRIE2.fit()
        
    Returns
    -------
    BRIE2 model object.
    """
    if Xc is None:
        Xc = np.ones((_count_layers[0].shape[0], 0), np.float32)
    if Xg is None:
        Xg = np.ones((_count_layers[0].shape[1], 0), np.float32)
                
    model = BRIE2(Nc=Xc.shape[0], Ng=Xg.shape[0],
                  Kc=Xc.shape[1], Kg=Xg.shape[1], 
                  intercept=intercept, intercept_mode=intercept_mode, 
                  effLen=effLen)
    
    losses = model.fit(data, Xc = Xc, Xg = Xg, **keyargs)
    
    ## Perform ELBO gain in the analogy to likelihood ratio
    if LRT_index is None:
        LRT_index = range(Xc.shape[1])
        
    if len(LRT_index) == 0:
        return model
        
    ELBO_gain = np.zeros((data[0].shape[1], len(LRT_index)), dtype=np.float32)
    for ii, idx in enumerate(LRT_index):
        if verbosity == 3:
            print("Fitting null model without feature %d" %(idx))
        Xc_del = np.delete(Xc, idx, 1)
        model_test = BRIE2(Nc=Xc_del.shape[0], Ng=data[0].shape[1],
                           Kc=Xc_del.shape[1], Kg=Xg.shape[1],
                           intercept = intercept, effLen=effLen)
        
        losses = model_test.fit(data, Xc = Xc_del, Xg = Xg, **keyargs)
        ELBO_gain[:, ii] = model_test.loss_gene - model.loss_gene
            
    model.ELBO_gain = ELBO_gain # H1 vs NUll
    model.pval = chi2.sf(2 * ELBO_gain, df = 1)
    # model_real.pval_log10 = chi2.logsf(-2 * LR, df = 1) * np.log10(np.exp(1))
    
    fdr = np.zeros(ELBO_gain.shape)
    for i in range(fdr.shape[1]):
        fdr[:, i] = multipletests(model.pval[:, i], method="fdr_bh")[1]
    model.fdr = fdr
    
    # return BRIE model
    return model


def fitBRIE(adata, Xc=None, Xg=None, intercept=None, intercept_mode='gene', 
            LRT_index=[], layer_keys=['isoform1', 'isoform2', 'ambiguous'], 
            **keyargs):
    """Fit a BRIE model from AnnData with cell and/or gene features
    
    Parameters
    ----------
    adata : :class:`~anndata.AnnData`, `np.ndarray`, `sp.spmatrix`
        The (annotated) data matrix of shape `n_obs` × `n_vars`. Rows correspond
        to cells and columns to genes.
    Xc : `numpy.array`, (n_feature, n_cell), np.float32, optional
        Cell features.
    Xg : `numpy.array`, (n_gene, n_feature), np.float32, optional
        Gene features.
    add_intercept : bool
        If add intecept for cell feature
    do_LRT : bool
        If perform likelihood ratio test for each gene
    layer_keys : list of `str` with length == 2 or 3
        isoform1, isoform2, [ambiguous]. Ambiguous count is optional.
        
    **keyargs : keyargs for BRIE2.fit()
        
    Returns
    -------
    BRIE2 model object. Also, adds 
        `Xc` and `weight_g` to `adata.obsm`; 
        `Xg` and `weight_c` to `adata.varm`;
        `Psi` and `Psi_var` to `adata.layers`.
    """
    #TODO: support sparse matrix!!
    _count_layers = [adata.layers[_key].toarray() for _key in layer_keys]
    
    if Xc is None:
        Xc = np.ones((_count_layers[0].shape[0], 0), np.float32)
    if Xg is None:
        Xg = np.ones((_count_layers[0].shape[1], 0), np.float32)
            
    _effLen = adata.varm['effLen'] if 'effLen' in adata.varm else None
    
    model = fit_BRIE_matrix(_count_layers, Xc=Xc, Xg=Xg, intercept=intercept,
                            intercept_mode=intercept_mode, LRT_index=LRT_index, 
                            effLen=_effLen, **keyargs)
    
    # update adata
    if Xc.shape[0] > 0:
        adata.obsm['Xc'] = Xc
        adata.varm['cell_coeff'] = model.Wc_loc.numpy().T
        
    if Xg.shape[1] > 0:
        adata.varm['Xg'] = Xg
        adata.obsm['gene_coeff'] = model.Wg_loc.numpy()
        
    if model.intercept_mode == 'gene':
        adata.varm['intercept'] = model.intercept.numpy().T
    elif model.intercept_mode == 'cell':
        adata.obsm['intercept'] = model.intercept.numpy()
        
    adata.varm['sigma'] = model.sigma.numpy().T
        
    # introduce sparse matrix for this
    adata.layers['Psi'] = model.Psi.numpy()
    adata.layers['Z_std'] = np.exp(model.Z_std.numpy())
    
    # losses
    adata.uns['brie_losses'] = model.losses.numpy()
    
    if LRT_index is None or len(LRT_index)>1:
        adata.varm['fdr'] = model.fdr
        adata.varm['pval'] = model.pval
        adata.varm['ELBO_gain'] = model.ELBO_gain
    
    # return BRIE model
    return model
