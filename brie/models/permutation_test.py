import time
import numpy as np
import tensorflow as tf
from .model_TFProb import BRIE2

def permutation_test(Xc, Rmat, P_iso1, P_iso2, 
                     test_index=[0], n_permute=100):
    """Permutation test
    """
    start_time = time.time()

    ## Fit model with real data
    model_real = BRIE2(Nc=Xc.shape[1], Ng=P_iso1.shape[0],
                       Kc=Xc.shape[0], Kg=0)
    
    losses = tf.concat([
            model_real.fit(Xc = Xc, Rs = Rmat,
                           P_iso1 = P_iso1, P_iso2 = P_iso2,
                           num_steps=100, learn_rate=0.02, size=10),
            model_real.fit(Xc = Xc, Rs = Rmat,
                           P_iso1 = P_iso1, P_iso2 = P_iso2,
                           num_steps=100, learn_rate=0.05, size=10),
            model_real.fit(Xc = Xc, Rs = Rmat,
                           P_iso1 = P_iso1, P_iso2 = P_iso2,
                           num_steps=300, learn_rate=0.1, size=10)
        ], axis=0)
    loss_real = losses[-1]
    weight_real = model_real.cell_weight.mean().numpy()
    
    print("Fit real data: %.2f min" % ((time.time() - start_time)/60))
    print(loss_real)
    
    
    ## Fit model with permulated data
    start_time = time.time()
    weigt_perm = np.zeros((weight_real.shape[0], 
                           weight_real.shape[1], n_permute))
    loss_perm = np.zeros(n_permute)
    for ir in range(n_permute):
        model_perm = BRIE2(Nc=Xc.shape[1], Ng=P_iso1.shape[0],
                           Kc=Xc.shape[0], Kg=0)
        _idx = np.random.permutation(Xc.shape[1])
        Xc_perm = Xc.copy()
        Xc_perm[test_index[0], _idx]
        
        losses = tf.concat([
            model_perm.fit(Xc = Xc, Rs = Rmat,
                           P_iso1 = P_iso1, P_iso2 = P_iso2,
                           num_steps=100, learn_rate=0.02, size=10),
            model_perm.fit(Xc = Xc, Rs = Rmat,
                           P_iso1 = P_iso1, P_iso2 = P_iso2,
                           num_steps=100, learn_rate=0.05, size=10),
            model_perm.fit(Xc = Xc, Rs = Rmat,
                           P_iso1 = P_iso1, P_iso2 = P_iso2,
                           num_steps=300, learn_rate=0.1, size=10)
        ], axis=0)

        weigt_perm[:, :, ir] = model_perm.cell_weight.mean().numpy()
        loss_perm[ir] = losses[-1]
    print("Fit permutated data: %.2f min" % ((time.time() - start_time)/60))
    
    return weight_real, loss_real, weight_perm, loss_perm
