import time
import numpy as np
import tensorflow as tf
from scipy.stats import chi2
from statsmodels.stats.multitest import multipletests

from .model_TFProb import BRIE2

def LikRatio_test(count_layers, Xc, p_ambiguous=None, index=None,
                  layer_ids=['isoform1', 'isoform2', 'ambiguous'], 
                  learn_steps=[100, 200, 300], 
                  learn_rates=[0.02, 0.5, 0.1], **kwargs):
    """likelihood ratio test
    """
    start_time = time.time()
    
    if isinstance(learn_steps, int):
        learn_steps = [learn_steps]
    if isinstance(learn_rates, float) or isinstance(learn_rates, int):
        learn_rates = [learn_rates]
    
    if (len(learn_steps) != len(learn_rates)):
        print("Error: len(learn_steps) != len(learn_rates).")
        return None
    
    ## load the dict
#     _count_layers = {}
#     _keys = _count_layers.keys()
#     for i in range(3):
#         _count_layers[_keys[i]] = count_layers[layer_ids[i]]
    _count_layers = count_layers.copy()

    ## Fit model with real data
    n_genes = count_layers[layer_ids[0]].shape[1]
    
    model_real = BRIE2(Nc=Xc.shape[1], Ng=n_genes,
                       Kc=Xc.shape[0], Kg=0, 
                       p_ambiguous=p_ambiguous)
    
#     losses = model_real.fit(_count_layers, Xc = Xc, **kwargs)
#     model_real.losses = losses

    loss_list = []
    for i in range(len(learn_steps)):
        loss_list.append(
            model_real.fit(_count_layers, Xc = Xc, max_iter=learn_steps[i], 
                           learn_rate=learn_rates[i], **kwargs)
        )
    model_real.losses = tf.concat(loss_list, axis=0)
    
    weight_real = model_real.Wc_loc.numpy()
    
    print("Fit real data: %.2f min" % ((time.time() - start_time)/60))
    print(model_real.losses[-1])
    
    ## Fit model with one removed feature
    start_time = time.time()
    if index is None:
        index = range(Xc.shape[0])
        
    LR = np.zeros((weight_real.shape[0], len(index)), dtype=np.float32)
    for ii, idx in enumerate(index):
        print("Fitting null model with feature %d" %(idx))
        Xc_del = np.delete(Xc, idx, 0)
        model_test = BRIE2(Nc=Xc_del.shape[1], Ng=n_genes,
                           Kc=Xc_del.shape[0], Kg=0, 
                           p_ambiguous=p_ambiguous)
        
#         losses = model_test.fit(_count_layers, Xc = Xc_del, **kwargs)

        for i in range(len(learn_steps)):
            model_test.fit(_count_layers, Xc = Xc_del, max_iter=learn_steps[i], 
                           learn_rate=learn_rates[i], **kwargs)
        
        LR[:, ii] = model_real.loss_gene - model_test.loss_gene
        
    print("Fit deleted data: %.2f min" % ((time.time() - start_time)/60))
    
    model_real.LR = LR # NUll vs H1
    model_real.pval_log10 = chi2.logsf(-2 * LR, df = 1) * np.log10(np.exp(1))
    
    fdr = np.zeros(LR.shape)
    for i in range(fdr.shape[1]):
        fdr[:, i] = multipletests(10**(model_real.pval_log10[:, i]), 
                                  method="fdr_bh")[1]
    model_real.fdr = fdr
    
    return model_real
