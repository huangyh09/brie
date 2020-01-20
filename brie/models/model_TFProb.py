## BRIE model

import numpy as np
import tensorflow as tf
import tensorflow_probability as tfp
from tensorflow_probability import distributions as tfd

class BRIE2(tf.keras.Model):
    """
    Ng : number of genes
    Nc : number of cells
    Kg : number of gene features
    Kc : number of cell features
    """
    def __init__(self, Nc, Ng, Kg, Kc, name=None):
        super(BRIE2, self).__init__(name=name)
        self.Wg_loc = tf.Variable(tf.random.normal([Kg, Nc]), name='Wg_loc')
        self.Wg_std = tf.Variable(tf.random.normal([Kg, Nc]), name='Wg_std')
        self.Wc_loc = tf.Variable(tf.random.normal([Ng, Kc]), name='Wc_loc')
        self.Wc_std = tf.Variable(tf.random.normal([Ng, Kc]), name='Wc_std')
        # self.var_s1 = tf.Variable(tf.exp(tf.random.normal([1])), name='var_s1')
        # self.var_s2 = tf.Variable(tf.exp(tf.random.normal([1])), name='var_s2')
    
    @property
    def cell_weight(self):
        """Variational posterior for the weight"""
        return tfd.Normal(self.Wc_loc, tf.exp(self.Wc_std))
    
    @property
    def gene_weight(self):
        """Variational posterior for the weight"""
        return tfd.Normal(self.Wg_loc, tf.exp(self.Wg_std))

    # @property
    # def std(self):
    #     """Variational posterior for the noise standard deviation"""
    #     return tfd.InverseGamma(tf.exp(self.var_s1), tf.exp(self.var_s2))
    
    @property
    def KLdiverg(self):
        """Sum of KL divergences between posteriors and priors"""
        prior = tfd.Normal(0, 1)
        return (tf.reduce_sum(tfd.kl_divergence(self.cell_weight, prior)) +
                tf.reduce_sum(tfd.kl_divergence(self.gene_weight, prior)))
    
    def logLik(self, Xc, Xg, Rs, P_iso1, P_iso2, sampling=False):
        """Get logLikelihood
        TODO: add sampling to Y, by introducing the variance term for y
        """
        sample = lambda x: x.sample() if sampling else x.mean()
        _loc = (tf.matmul(sample(self.cell_weight), Xc) + 
                tf.matmul(Xg, sample(self.gene_weight)))
        # _std = tf.sqrt(sample(self.std))
        
        # Y = tfd.Normal(_loc, _std)
        # Psi = tf.sigmoid(sample(Y))
        Psi = tf.sigmoid(_loc)
        
        logLik = (Rs['1'] * tf.math.log(P_iso1[:, 0:1] * Psi) + 
                  Rs['2'] * tf.math.log(P_iso2[:, 1:2] * (1 - Psi)) + 
                 (Rs['3'] * tf.math.log(P_iso1[:, 2:3] * Psi + 
                                        P_iso2[:, 2:3] * (1 - Psi))))
        
        return logLik

    def fit(self, Xc, Xg, Rs, P_iso1, P_iso2, num_steps=100, optimizer=None,
            learn_rate=0.05):
        """Fit the model's parameters"""
        if optimizer is None:
            optimizer = tf.optimizers.Adam(learning_rate=learn_rate)
            
        loss_fn = lambda: (self.KLdiverg - 
                           tf.reduce_sum(self.logLik(Xc, Xg, Rs, P_iso1, P_iso2)))
        
        losses = tfp.math.minimize(loss_fn, 
                                   num_steps=num_steps, 
                                   optimizer=optimizer)
        return losses

