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

class Result_list:
    def __init__(self, model_list, idx_list=None, do_LRT=True):
        fdr, pval, ELBO_gain = [], [], []
        cell_coeff, intercept, sigma = [], [], []
        Psi, Z_std, row_idx, col_idx = [], [], [], []
        for ii, res in enumerate(model_list):
            cell_coeff.append(list(res.Wc_loc.numpy()[0, :]))
            intercept.append(res.intercept[0, 0])
            sigma.append(res.sigma[0, 0])
            if do_LRT:
                fdr.append(list(res.fdr[0, :]))
                pval.append(list(res.pval[0, :]))
                ELBO_gain.append(list(res.ELBO_gain[0, :]))

            Psi += list(res.Psi.numpy().transpose()[0, :])
            Z_std += list(res.Z_std.numpy().transpose()[0, :])
            col_idx += [ii] * len(res.Wc_loc.shape[1])
            row_idx += idx_list[ii]

        self.cell_coeff = np.array(cell_coeff, dtype=np.float32)
        self.intercept = np.array(intercept, dtype=np.float32)
        self.sigma = np.array(sigma, dtype=np.float32)
        self.Psi = csc_matrix((Psi, (row_idx, col_idx)), shape=data[0].shape)
        self.Z_loc = csc_matrix((Z_loc, (row_idx, col_idx)), shape=data[0].shape)
        if do_LRT:
            self.fdr = np.array(fdr, dtype=np.float32)
            self.pval = np.array(pval, dtype=np.float32)
            self.ELBO_gain = np.array(ELBO_gain, dtype=np.float32)
            
        

def fit_BRIE_matrix(data, Xc=None, Xg=None, intercept=None, p_ambiguous=None, 
                    do_LRT=True, LRT_index=None, **keyargs):
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
        Xc = np.ones((0, data[0].shape[0]), np.float32)
    if Xg is None:
        Xg = np.ones((data[0].shape[1], 0), np.float32)
                
    model = BRIE2(Nc=Xc.shape[1], Ng=Xg.shape[0],
                  Kc=Xc.shape[0], Kg=Xg.shape[1], 
                  intercept = intercept, p_ambiguous=p_ambiguous)
    
    losses = model.fit(data, Xc = Xc, Xg = Xg, **keyargs)
    
    if do_LRT == False:
        return model
    
    ## Perform ELBO gain in the analogy to likelihood ratio
    if LRT_index is None:
        LRT_index = range(Xc.shape[0])
        
    ELBO_gain = np.zeros((data[0].shape[1], len(LRT_index)), dtype=np.float32)
    for ii, idx in enumerate(LRT_index):
#         if verbosity == 3:
#             print("Fitting null model with feature %d" %(idx))
        print("Fitting null model with feature %d" %(idx))
        Xc_del = np.delete(Xc, idx, 0)
        model_test = BRIE2(Nc=Xc_del.shape[1], Ng=data[0].shape[1],
                           Kc=Xc_del.shape[0], Kg=Xg.shape[1],
                           intercept = intercept, p_ambiguous=p_ambiguous)
        
        losses = model_test.fit(data, Xc = Xc_del, Xg = Xg, **keyargs)
        ELBO_gain[:, ii] = model.loss_gene - model_test.loss_gene
            
    model.ELBO_gain = ELBO_gain # H1 vs NUll
    model.pval = chi2.sf(2 * ELBO_gain, df = 1)
    # model_real.pval_log10 = chi2.logsf(-2 * LR, df = 1) * np.log10(np.exp(1))
    
    fdr = np.zeros(ELBO_gain.shape)
    for i in range(fdr.shape[1]):
        fdr[:, i] = multipletests(model.pval[:, i], method="fdr_bh")[1]
    model.fdr = fdr
    
    # return BRIE model
    return model


def show_progress(RV):
    return RV


def fit_BRIE_mproc(data, Xc=None, Xg=None, intercept=None, p_ambiguous=None, 
                   do_LRT=True, LRT_index=None, nproc=10, **keyargs):
    """Fit BRIE gene wise, which doesn't support gene feature
    """
    if Xc is None:
        Xc = np.ones((0, data[0].shape[0]), np.float32)
        
    Xg = np.ones((data[0].shape[1], 0), np.float32)
    
    results = []
    idx_list = []
    #pool = multiprocessing.Pool(processes=nproc)
    executor = concurrent.futures.ProcessPoolExecutor(nproc)
    for g in range(data[0].shape[1]):
        _p_ambiguous = None if p_ambiguous is None else p_ambiguous[g:g+1, :]
        
        _count = data[0][:, g] + data[1][:, g]
        # if len(data) > 2 and _p_ambiguous[0] != 0.5:
        #     _count += data[2][:, g]
            
        _idx = _count > 0
        idx_list.append(_idx)
        # if sum(idx) <= 10:
        #     results.append(None)
        #     continue
            
        _data = [x[_idx, g:g+1] for x in data]
        
        
        print(g, sum(_idx))
        
#         results.append(fit_BRIE_matrix(_data, Xc=Xc[:, _idx], Xg=None, 
#                        intercept=intercept, 
#                        p_ambiguous=_p_ambiguous, do_LRT=do_LRT, 
#                        LRT_index=LRT_index, verbose=False, **keyargs))
        
        
        results.append(executor.submit(
            fit_BRIE_matrix, _data, Xc=Xc[:, _idx], Xg=None, 
            intercept=intercept, p_ambiguous=_p_ambiguous, 
            do_LRT=do_LRT, LRT_index=LRT_index, verbose=False, **keyargs
        ))#
        
    concurrent.futures.wait(results)

#         results.append(pool.apply_async(fit_BRIE_matrix,
#                                         (_data, Xc[:, _idx], None, intercept, 
#                                          _p_ambiguous, do_LRT, LRT_index), 
#                                         #dict(**keyargs), 
#                                         callback=show_progress))
#         pool.close()
#         pool.join()
        
#     results = [res.get() for res in results]

    results = [res.result() for res in results]
    RV = BRIEresult(resutls, idx_list, do_LRT)
    RV.Xc = Xc
    
    return RV
    

#         _model = fit_BRIE_model(_data, Xc=Xc[:, _idx], Xg=None, 
#                                 intercept=intercept, 
#                                 p_ambiguous=_p_ambiguous, do_LRT=do_LRT, 
#                                 LRT_index=LRT_index, **keyargs)