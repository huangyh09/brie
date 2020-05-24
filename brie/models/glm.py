import numpy as np
from scipy.stats import chi2
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

class glm_res():
    def __int__(self):
        pass

def GLMfit(adata, Xc=None, index=None, add_intercept=True, family='Binomial', 
           min_non_zero=1, **keyargs):
    """Fit a BRIE model with cell and/or gene features
    
    Parameters
    ----------
    data : :class:`~anndata.AnnData`, `np.ndarray`, `sp.spmatrix`
        The (annotated) data matrix of shape `n_obs` × `n_vars`. Rows correspond
        to cells and columns to genes.
    Xc : `numpy.array`, (n_feature, n_cell), np.float32, optional
        Cell features.
    **keyargs : keyargs for BRIE2.fit()
        
    Returns
    -------
    BRIE2 model object. Also, adds 
        `Xc` and `weight_g` to `adata.obsm`; 
        `Xg` and `weight_c` to `adata.varm`;
        `Psi` and `Psi_var` to `adata.layers`.
    """
    if Xc is None:
        Xc = np.ones((adata.shape[0], 0), np.float32)
    if add_intercept and sum((Xc == 1).mean(0) == 1) == 0:
        Xc = np.append(Xc, np.ones((adata.shape[0], 1), np.float32), axis=1)
        
    if family == 'Binomial':
        _family = sm.families.Binomial()
    else:
        print("%s not supported yet" %(family))
        
    if index is None:
        index = range(Xc.shape[0])
        
    LR_mat = np.zeros((adata.shape[1], Xc.shape[1]))
    pv_mat = np.ones((adata.shape[1], Xc.shape[1]))
    fdr_mat = np.ones((adata.shape[1], Xc.shape[1]))
    coeff_mat = np.zeros((adata.shape[1], Xc.shape[1]))
    
    for i in range(adata.shape[1]):
        
        yy = np.append(adata.layers['isoform1'][:, i:(i+1)],
                       adata.layers['isoform2'][:, i:(i+1)], axis=1)
        _idx = np.sum(yy, axis=1) > 0
        xx = Xc[_idx, :]
        yy = yy[_idx, :]
        
        if (sum(np.min(yy, axis=1) > 0) < min_non_zero): continue
            
        glm_binom1 = sm.GLM(yy, xx, family=_family)
        try:
            res1 = glm_binom1.fit()
        except :
            continue

        coeff_mat[i, :] = res1.params
                
        # Likelihood ratio test
        for j in index:
            xx_del = np.delete(xx.copy(), j, 1)
            
            glm_binom0 = sm.GLM(yy, xx_del, family=_family)
            res0 = glm_binom0.fit()
            LR_mat[i, j] = res0.llf - res1.llf
            pv_mat[i, j] = chi2.sf(-2 * LR_mat[i, j], df = 1)
    
    for i in range(fdr_mat.shape[1]):
        fdr_mat[:, i] = multipletests(pv_mat[:, i], method="fdr_bh")[1]
    
    # return results
    RV = glm_res()
    RV.LR = LR_mat
    RV.pv = pv_mat
    RV.fdr = fdr_mat
    RV.coeff_mat = coeff_mat
    return RV
