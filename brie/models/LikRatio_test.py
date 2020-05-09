import time
import numpy as np
import tensorflow as tf
from scipy.stats import chi2
from statsmodels.stats.multitest import multipletests

from .model_TFProb import BRIE2

def pval_LR_test(LR, df=1):
    """Return the p values of loglikelihood ratio test
    """
    pval = chi2.logsf(2 * LR, df = df)

def LR_test(Xc, Rmat, P_iso1, P_iso2):
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
    
    model_real.losses = losses
    loss_real = losses[-1]
    weight_real = model_real.cell_weight.mean().numpy()
    
    print("Fit real data: %.2f min" % ((time.time() - start_time)/60))
    print(loss_real)
    
    
    ## Fit model with one removed feature
    start_time = time.time()
    LR = np.zeros(weight_real.shape)
    for i in range(Xc.shape[0]):
        model_test = BRIE2(Nc=Xc.shape[1], Ng=P_iso1.shape[0],
                           Kc=Xc.shape[0] - 1, Kg=0)
        Xc_del = np.delete(Xc, i, 0)
        
        losses = tf.concat([
            model_test.fit(Xc = Xc_del, Rs = Rmat,
                           P_iso1 = P_iso1, P_iso2 = P_iso2,
                           num_steps=100, learn_rate=0.02, size=10),
            model_test.fit(Xc = Xc_del, Rs = Rmat,
                           P_iso1 = P_iso1, P_iso2 = P_iso2,
                           num_steps=100, learn_rate=0.05, size=10),
            model_test.fit(Xc = Xc_del, Rs = Rmat,
                           P_iso1 = P_iso1, P_iso2 = P_iso2,
                           num_steps=300, learn_rate=0.1, size=10)
        ], axis=0)

        LR[:, i] = model_real.loss_gene - model_test.loss_gene
    print("Fit deleted data: %.2f min" % ((time.time() - start_time)/60))
    
    model_real.LR = LR # NUll vs H1
    model_real.pval_log = chi2.logsf(-2 * LR, df = 1) #* log10(exp(1))
    
    fdr = np.zeros(LR.shape)
    for i in range(fdr.shape[1]):
        fdr[:, i] = multipletests(np.exp(model_real.pval_log[:, i]))[1]
    model_real.fdr = fdr
    
    return model_real

