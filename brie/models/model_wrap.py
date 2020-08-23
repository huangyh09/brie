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


class BRIE_RV():
    """Return value object for BRIE2 model
    """
    def __init__(self, model):
        """Initialization with BRIE2 model
        """
        self.Nc = model.Nc
        self.Ng = model.Ng
        self.Kc = model.Kc
        self.Kg = model.Kg
        self.shape = (self.Nc, self.Ng)
        
        self.Xc = model.Xc
        self.Xg = model.Xg
        self.sigma = model.sigma.numpy()
        self.intercept = model.intercept.numpy()
        self.cell_coeff = model.Wc_loc.numpy()
        self.gene_coeff = model.Wg_loc.numpy()
        
        self.Psi = model.Psi.numpy()
        self.Psi95CI = model.Psi95CI
        self.Z_std = model.Z_std.numpy()
        self.losses = model.losses.numpy()
        self.intercept_mode = model.intercept_mode
        
        
    def __str__(self):
        return "BRIE2 results for %d cells and %d genes" %(self.Nc, self.Ng)
        
    def concate(self, new_RV, axis=1):
        """
        """
        if axis != 1:
            print("Warning: only suppoting gene level concate!")
            return None
        
        self.Ng += new_RV.Ng
        self.losses = np.append(self.losses, new_RV.losses)
        
        self.sigma = np.append(self.sigma, new_RV.sigma, axis=1)
        self.intercept = np.append(self.intercept, new_RV.intercept, axis=1)
        self.cell_coeff = np.append(self.cell_coeff, new_RV.cell_coeff, axis=1)
        self.Psi = np.append(self.Psi, new_RV.Psi, axis=1)
        self.Psi95CI = np.append(self.Psi95CI, new_RV.Psi95CI, axis=1)
        self.Z_std = np.append(self.Z_std, new_RV.Z_std, axis=1)
        
        if hasattr(new_RV, 'ELBO_gain'):
            self.fdr = np.append(self.fdr, new_RV.fdr, axis=0)
            self.pval = np.append(self.pval, new_RV.pval, axis=0)
            self.ELBO_gain = np.append(self.ELBO_gain, new_RV.ELBO_gain, axis=0)
        

def concate(BRIE_RV_list):
    """Concate a list of BRIE results
    """
    res_merge = BRIE_RV_list[0]
    for _res in BRIE_RV_list[1:]:
        res_merge.concate(_res)
        
    return res_merge


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
        Xc = np.ones((data[0].shape[0], 0), np.float32)
    if Xg is None:
        Xg = np.ones((data[0].shape[1], 0), np.float32)
    
    # if len(data) > 2:
    #     from .model_TFProb import BRIE2
    # else:
    #    print("BRIE2-Beta in use")
    #    from .model_Beta import BRIE2_Beta as BRIE2
    
    model = BRIE2(Nc=Xc.shape[0], Ng=Xg.shape[0],
                  Kc=Xc.shape[1], Kg=Xg.shape[1], 
                  effLen=effLen, intercept=intercept,
                  intercept_mode=intercept_mode)
    
    losses = model.fit(data, Xc = Xc, Xg = Xg, **keyargs)
    
    brie_results = BRIE_RV(model)
    
    ## Perform ELBO gain in the analogy to likelihood ratio
    if LRT_index is None:
        LRT_index = range(Xc.shape[1])
        
    if len(LRT_index) == 0:
        return brie_results
        
    ELBO_gain = np.zeros((data[0].shape[1], len(LRT_index)), dtype=np.float32)
    for ii, idx in enumerate(LRT_index):
        if verbosity == 3:
            print("[BRIE2] fitting null model without feature %d" %(idx))
        Xc_del = np.delete(Xc, idx, 1)
        model_test = BRIE2(Nc=Xc_del.shape[0], Ng=data[0].shape[1],
                           Kc=Xc_del.shape[1], Kg=Xg.shape[1],
                           intercept = intercept, effLen=effLen)
        
        losses = model_test.fit(data, Xc = Xc_del, Xg = Xg, **keyargs)
        ELBO_gain[:, ii] = model_test.loss_gene - model.loss_gene
            
    brie_results.ELBO_gain = ELBO_gain # H1 vs NUll
    brie_results.pval = chi2.sf(2 * ELBO_gain, df = 1)
    # model_real.pval_log10 = chi2.logsf(-2 * LR, df = 1) * np.log10(np.exp(1))
    
    fdr = np.zeros(ELBO_gain.shape)
    for i in range(fdr.shape[1]):
        fdr[:, i] = multipletests(brie_results.pval[:, i], method="fdr_bh")[1]
    brie_results.fdr = fdr
    
    # return BRIE model
    return brie_results


def fitBRIE(adata, Xc=None, Xg=None, intercept=None, intercept_mode='gene', 
            LRT_index=[], layer_keys=['isoform1', 'isoform2', 'ambiguous'], 
            batch_size=500000, **keyargs):
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
    if Xc is None:
        Xc = np.ones((adata.shape[0], 0), np.float32)
    if Xg is None:
        Xg = np.ones((adata.shape[1], 0), np.float32)
        
    if (Xg is None or Xg.shape[1] == 0) and intercept_mode.upper() != 'CELL':
        _n_gene = int(np.ceil(batch_size / adata.shape[0]))
        _n_batch = int(np.ceil(adata.shape[1] / _n_gene))
        res_list = []
        for i in range(_n_batch):
            _idx = range(_n_gene * i, min(_n_gene * (i+1), adata.shape[1]))
            _count_layers = [adata.layers[_key][:, _idx] for _key in layer_keys]
            _effLen = adata.varm['effLen'][_idx, :] if 'effLen' in adata.varm else None

            _ResVal = fit_BRIE_matrix(
                _count_layers, Xc=Xc, Xg=Xg[_idx, :], effLen=_effLen, 
                intercept=intercept, intercept_mode=intercept_mode, 
                LRT_index=LRT_index, **keyargs)
            
            res_list.append(_ResVal)
            print("[BRIE2] %d out %d genes done" 
                  %(min(_n_gene * (i + 1), adata.shape[1]), adata.shape[1]))
            
        ResVal = concate(res_list)
    else:
        _count_layers = [adata.layers[_key] for _key in layer_keys]
        _effLen = adata.varm['effLen'] if 'effLen' in adata.varm else None

        ResVal = fit_BRIE_matrix(
            _count_layers, Xc=Xc, Xg=Xg, effLen=_effLen, intercept=intercept, 
            intercept_mode=intercept_mode, LRT_index=LRT_index, **keyargs)
    
    # update adata
    if Xc.shape[0] > 0:
        adata.obsm['Xc'] = Xc
        adata.varm['cell_coeff'] = ResVal.cell_coeff.T
        
    if Xg.shape[1] > 0:
        adata.varm['Xg'] = Xg
        adata.obsm['gene_coeff'] = ResVal.gene_coeff
        
    if ResVal.intercept_mode == 'gene':
        adata.varm['intercept'] = ResVal.intercept.T
    elif ResVal.intercept_mode == 'cell':
        adata.obsm['intercept'] = ResVal.intercept
        
    adata.varm['sigma'] = ResVal.sigma.T
        
    # TODO: introduce sparse matrix for this
    adata.layers['Psi'] = ResVal.Psi
    adata.layers['Z_std'] = ResVal.Z_std
    adata.layers['Psi_95CI'] = ResVal.Psi95CI
    
    # losses
    adata.uns['brie_losses'] = ResVal.losses
        
    if LRT_index is None or len(LRT_index) >= 1:
        adata.varm['fdr'] = ResVal.fdr
        adata.varm['pval'] = ResVal.pval
        adata.varm['ELBO_gain'] = ResVal.ELBO_gain
    
    # return BRIE model result value
    return ResVal
