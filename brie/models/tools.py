# Tools for BRIE package

import numpy as np
from .model_TFProb import BRIE2
from .simulator import simulator

from scipy.stats import chi2
from statsmodels.stats.multitest import multipletests

def fitBRIE(adata, Xc=None, Xg=None, add_intercept=True, **keyargs):
    """Fit a BRIE model with cell and/or gene features
    
    Parameters
    ----------
    data : :class:`~anndata.AnnData`, `np.ndarray`, `sp.spmatrix`
        The (annotated) data matrix of shape `n_obs` × `n_vars`. Rows correspond
        to cells and columns to genes.
    Xc : `numpy.array`, (n_feature, n_cell), np.float32, optional
        Cell features.
    Xg : `numpy.array`, (n_gene, n_feature), np.float32, optional
        Gene features.
    **keyargs : keyargs for BRIE2.fit()
        
    Returns
    -------
    BRIE2 model object. Also, adds 
        `Xc` and `weight_g` to `adata.obsm`; 
        `Xg` and `weight_c` to `adata.varm`;
        `Psi` and `Psi_var` to `adata.layers`.
    """
    if Xc is None:
        Xc = np.ones((0, adata.shape[0]), np.float32)
    if Xg is None:
        Xg = np.ones((adata.shape[1], 0), np.float32)
        
    _intecept = None if add_intercept else 0
        
    model = BRIE2(Nc=adata.shape[0], Ng=adata.shape[1],
                  Kc=Xc.shape[0], Kg=Xg.shape[1], 
                  intercept = _intecept,
                  p_ambiguous=adata.varm['p_ambiguous'])
    
#     _layers = {}
#     _layers['1'] = adata.layers['isoform1']
#     _layers['2'] = adata.layers['isoform2']
#     _layers['3'] = adata.layers['ambiguous']
#     losses = model.fit(_layers, Xc = Xc, Xg = Xg, **keyargs)
    
    losses = model.fit(adata.layers, Xc = Xc, Xg = Xg, **keyargs)
    
    # update adata
    if Xc.shape[0] > 0:
        adata.obsm['Xc'] = Xc.transpose()
        adata.varm['cell_coeff'] = model.Wc_loc.numpy()
        
    if Xg.shape[1] > 0:
        adata.varm['Xg'] = Xg
        adata.obsm['gene_coeff'] = model.Wg_loc.numpy().transpose()
        
    adata.varm['intercept'] = model.intercept.numpy()
    adata.varm['sigma'] = model.sigma.numpy()
        
    adata.layers['Psi'] = model.Psi.numpy().transpose()
    adata.layers['Z_std'] = np.exp(model.Z_std.numpy()).transpose()
    
    # return BRIE model
    return model



def LRTest(adata, Xc, Xg=None, index=None, add_intercept=True, **kwargs):
    """likelihood ratio test
    
    layer_keys: ['isoform1', 'isoform2', 'ambiguous']
    target: marginLik, ELBO
    """
    ## Fit the full model
    model_real = fitBRIE(adata, Xc=Xc, Xg=Xg, add_intercept=add_intercept, 
                         **kwargs)
        
    ## Fit model with one removed feature
    if index is None:
        index = range(Xc.shape[0])
        
    _intecept = None if add_intercept else 0
        
    LR = np.zeros((adata.shape[1], len(index)), dtype=np.float32)
    for ii, idx in enumerate(index):
        print("Fitting null model with feature %d" %(idx))
        Xc_del = np.delete(Xc, idx, 0)
        model_test = BRIE2(Nc=Xc_del.shape[1], Ng=adata.shape[1],
                           Kc=Xc_del.shape[0], Kg=0, 
                           intercept = _intecept,
                           p_ambiguous=adata.varm['p_ambiguous'])
        
        losses = model_test.fit(adata.layers, Xc = Xc_del, **kwargs)
        LR[:, ii] = model_real.loss_gene - model_test.loss_gene
            
    model_real.LR = LR # NUll vs H1
    model_real.pval = chi2.sf(-2 * LR, df = 1)
    # model_real.pval_log10 = chi2.logsf(-2 * LR, df = 1) * np.log10(np.exp(1))
    
    fdr = np.zeros(LR.shape)
    for i in range(fdr.shape[1]):
        fdr[:, i] = multipletests(model_real.pval[:, i], method="fdr_bh")[1]
    model_real.fdr = fdr
    
    adata.varm['LR'] = model_real.LR
    adata.varm['fdr'] = model_real.fdr
    adata.varm['pval'] = model_real.pval
    
    return model_real


def perm_test(adata, Xc, Xg=None, index=None, n_sample=500, add_intercept=True, 
              **kwargs):
    """likelihood ratio test
    
    layer_keys: ['isoform1', 'isoform2', 'ambiguous']
    target: marginLik, ELBO
    """
    ## Fit the full model
    model_real = fitBRIE(adata, Xc=Xc, Xg=Xg, add_intercept=add_intercept, 
                         **kwargs)
    
    ## Fit model with permuting one feature
    if index is None:
        index = range(Xc.shape[0])
        
    adata_copy = adata.copy()
    Wc_perm = np.zeros((n_sample, adata.shape[1], len(index)), dtype=np.float32)
    for ii, idx in enumerate(index):
        print("Fitting null model with feature %d" %(idx))
        Xc_perm = Xc.copy()
        for ir in range(n_sample):
            Xc_perm[idx, :] = np.random.permutation(Xc_perm[idx, :])
            model_test = fitBRIE(adata_copy, Xc=Xc_perm, Xg=Xg, 
                                 add_intercept=add_intercept, **kwargs)
            Wc_perm[ir, :, ii] = model_test.Wc_loc.numpy()[:, ii]
            
    Wc_real = np.expand_dims(model_real.Wc_loc.numpy()[:, index], 0)
    quantile = np.mean((Wc_real - Wc_perm) > 0, axis=0)
    pval = (0.5 - np.abs(quantile - 0.5)) * 2
    
    fdr = np.zeros(pval.shape, dtype=np.float32)
    for i in range(fdr.shape[1]):
        fdr[:, i] = multipletests(pval[:, i], method="fdr_bh")[1]
    
    model_real.fdr_perm = fdr
    model_real.pval_perm = pval
    model_real.quantile = quantile
    
    adata.varm['quantile'] = model_real.quantile
    adata.varm['fdr_perm'] = model_real.fdr_perm
    adata.varm['pval_perm'] = model_real.pval_perm
    
    return model_real
