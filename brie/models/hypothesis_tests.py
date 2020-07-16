import numpy as np
from scipy.stats import chi2
from statsmodels.stats.multitest import multipletests

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

