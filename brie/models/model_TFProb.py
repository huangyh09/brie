## BRIE model

import numpy as np
import tensorflow as tf
import tensorflow_probability as tfp
from tensorflow_probability import distributions as tfd

class BRIE2():
    """
    Ng : number of genes
    Nc : number of cells
    Kg : number of gene features
    Kc : number of cell features
    """
    # TODO: support without gene and / or cell feature
    def __init__(self, Nc, Ng, Kg, Kc, name=None):
        self.Nc = Nc
        self.Ng = Ng
        self.Kc = Kc
        self.Kg = Kg

        self.Wg_loc = tf.Variable(tf.random.normal([Kg, Nc]), name='Wg_loc')
        self.Wg_std = tf.Variable(tf.random.normal([Kg, Nc]), name='Wg_std')
        self.Wc_loc = tf.Variable(tf.random.normal([Ng, Kc]), name='Wc_loc')
        self.Wc_std = tf.Variable(tf.random.normal([Ng, Kc]), name='Wc_std')
        self.var_s1 = tf.Variable(tf.exp(tf.random.normal([1])), name='var_s1')
        self.var_s2 = tf.Variable(tf.exp(tf.random.normal([1])), name='var_s2')
    
    @property
    def cell_weight(self):
        """Variational posterior for the weight"""
        return tfd.Normal(self.Wc_loc, tf.exp(self.Wc_std))
    
    @property
    def gene_weight(self):
        """Variational posterior for the weight"""
        return tfd.Normal(self.Wg_loc, tf.exp(self.Wg_std))

    @property
    def std(self):
        """Variational posterior for the noise standard deviation"""
        return tfd.InverseGamma(tf.exp(self.var_s1), tf.exp(self.var_s2))
    
    @property
    def KLdiverg(self):
        """Sum of KL divergences between posteriors and priors"""
        prior = tfd.Normal(0, 1)
        return (tf.reduce_sum(tfd.kl_divergence(self.cell_weight, prior)) +
                tf.reduce_sum(tfd.kl_divergence(self.gene_weight, prior)))
    
    def logLik(self, Xc, Xg, Rs, P_iso1, P_iso2, sampling=True, size=10):
        """Get marginal logLikelihood via variaitonal posterior
        """
        # TODO: introduce sparse tensor here
        # There will be two layers of sampling: weights and z
        # tensor shapes: (n_z, n_w, n_gene, n_cell)
        if sampling:
            _W_cell = self.cell_weight.sample(size) # (n_w, n_g, k_c)
            _W_gene = self.gene_weight.sample(size) # (n_w, k_g, n_c)
            _Z_var = tf.sqrt(self.std.sample((size, 1)))
            
            _Z_loc = (
                tf.tensordot(_W_cell, Xc, [2, 0]) + 
                tf.transpose(tf.tensordot(Xg, _W_gene, [1, 1]), [1, 0, 2]))
            _Z = tfd.Normal(_Z_loc, _Z_var).sample(size)
        else:
            _W_cell = tf.reshape(self.cell_weight.mean(), ())
            _W_gene = tf.reshape(self.gene_weight.mean(), ())
            _Z = tf.tensordot(_W_cell, Xc) + tf.tensordot(Xg, _W_gene)

        # Transform Z to Psi and log Psi for isoform 1 and 2
        Psi1 = tf.sigmoid(_Z)
        Psi2 = 1 - Psi1
        Psi1_log = tf.math.log_sigmoid(_Z)
        Psi2_log = tf.math.log_sigmoid(0 - _Z)
                
        # Calculate element wise logLikelihood
        def _re2(x):
            return tf.reshape(x, (1, 1, self.Ng, self.Nc))
        def _re3(x): 
            return tf.reshape(x, (1, 1, self.Ng, 1))
        
        _logLik = (
            _re2(Rs['1']) * (tf.math.log(_re3(P_iso1[:, 0])) + Psi1_log) + 
            _re2(Rs['2']) * (tf.math.log(_re3(P_iso2[:, 1])) + Psi2_log) + 
            _re2(Rs['3']) * (tf.math.log(
                _re3(P_iso1[:, 2]) * Psi1 + _re3(P_iso2[:, 2]) * Psi2)))
        
        return tf.math.reduce_mean(
            tf.math.reduce_logsumexp(_logLik, axis=0), axis=1)

    def fit(self, Xc, Xg, Rs, P_iso1, P_iso2, num_steps=100, optimizer=None,
            learn_rate=0.05, **kwargs):
        """Fit the model's parameters"""
        if optimizer is None:
            optimizer = tf.optimizers.Adam(learning_rate=learn_rate)
            
        loss_fn = lambda: (self.KLdiverg - 
                           tf.reduce_sum(self.logLik(Xc, Xg, Rs, P_iso1, 
                                                     P_iso2, **kwargs)))
        
        losses = tfp.math.minimize(loss_fn, 
                                   num_steps=num_steps, 
                                   optimizer=optimizer)
        return losses

