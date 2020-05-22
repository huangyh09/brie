import time
import numpy as np
import tensorflow as tf
from .model_TFProb import BRIE2


def permutation_test(count_layers, Xc, p_ambiguous=None, index=None,
                     layer_ids=['isoform1', 'isoform2', 'ambiguous'], 
                     learn_steps=[100, 200, 300], learn_rates=[0.02, 0.5, 0.1], 
                     n_permute=100, **kwargs):
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
    _count_layers = {}
    for i in range(3):
        _count_layers[str(i + 1)] = count_layers[layer_ids[i]]

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
        
    pval_perm = np.ones((weight_real.shape[0], len(index)))
    weigt_perm = np.zeros((weight_real.shape[0], len(index), n_permute))
    loss_perm = np.zeros((len(index), n_permute))
    for ii, idx in enumerate(index):
        print("Fitting null model with feature %d" %(idx))
        
        for ir in range(n_permute):
            model_perm = BRIE2(Nc=Xc.shape[1], Ng=P_iso1.shape[0],
                               Kc=Xc.shape[0], Kg=0)
            _idx_perm = np.random.permutation(Xc.shape[1])
            Xc_perm = Xc.copy()
            Xc_perm[idx, :] = Xc_perm[idx, _idx_perm]

            # losses = model_perm.fit(_count_layers, Xc = Xc_del, **kwargs)

            for i in range(len(learn_steps)):
                model_perm.fit(_count_layers, Xc = Xc_perm, 
                               max_iter=learn_steps[i], 
                               learn_rate=learn_rates[i], **kwargs)

            weigt_perm[:, ii, ir] = model_perm.cell_weight.mean().numpy()[:, idx]
            loss_perm[ii, ir] = model_perm.losses[-1]
            
        pval_perm[:, ii] = np.mean(weight_real[:, idx:(idx+1)] < 
                                   weigt_perm[:, ii, :])
    pval_perm = (0.5 - np.abs(pval_perm - 0.5)) * 2
        
    print("Fit deleted data: %.2f min" % ((time.time() - start_time)/60))
    
    model_real.weight_perm = weight_perm
    model_real.loss_perm = loss_perm
    model_real.pval_perm = pval_perm
    
    fdr = np.zeros(pval_perm.shape)
    for i in range(fdr.shape[1]):
        fdr[:, i] = multipletests(model_real.pval_perm[:, i], 
                                  method="fdr_bh")[1]
    model_real.fdr = fdr
    
    return model_real
