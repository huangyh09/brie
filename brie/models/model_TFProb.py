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
    def __init__(self, Nc, Ng, Kc, Kg=0, name=None):
        self.Nc = Nc
        self.Ng = Ng
        self.Kc = Kc
        self.Kg = Kg
        
        self.Z_loc = tf.Variable(tf.random.normal([Ng, Nc]), name='Z_loc')
        self.Z_std = tf.Variable(tf.random.normal([Ng, Nc]), name='Z_var')
        self.Wg_loc = tf.Variable(tf.random.normal([Kg, Nc]), name='Wg_loc')
        self.Wg_std = tf.Variable(tf.random.normal([Kg, Nc]), name='Wg_std')
        self.Wc_loc = tf.Variable(tf.random.normal([Ng, Kc]), name='Wc_loc')
        self.Wc_std = tf.Variable(tf.random.normal([Ng, Kc]), name='Wc_std')
        self.var_s1 = tf.Variable(tf.ones([Ng, 1]) * 3, name='var_s1')
        self.var_s2 = tf.Variable(tf.ones([Ng, 1]) * 1, name='var_s2')

        self.set_prior()

    def set_prior(self, 
            w_prior=tfd.Normal(0, 1), 
            std_prior=tfd.Gamma(7, 3)):
        """Set priors distributions"""
        self.w_prior = w_prior
        self.std_prior = std_prior
    
    @property
    def cell_weight(self):
        """Variational posterior for the weight"""
        return tfd.Normal(self.Wc_loc, tf.exp(self.Wc_std))
    
    @property
    def gene_weight(self):
        """Variational posterior for the weight"""
        return tfd.Normal(self.Wg_loc, tf.exp(self.Wg_std))

    @property
    def Z(self):
        """Variational posterior for the logit Psi"""
        return tfd.Normal(self.Z_loc, tf.exp(self.Z_std))
    
    @property
    def std(self):
        """Variational posterior for the variance of regression residue"""
        return tfd.Gamma(tf.exp(self.var_s1), tf.exp(self.var_s2))
    
    @property
    def KLdiverg(self):
        """Sum of KL divergences between posteriors and priors"""
        return (
            tf.reduce_sum(tfd.kl_divergence(self.cell_weight, self.w_prior)) +
            tf.reduce_sum(tfd.kl_divergence(self.gene_weight, self.w_prior)) +
            tf.reduce_sum(tfd.kl_divergence(self.std, self.std_prior)))
    
    def logLik(self, Xc, Xg, Rs, P_iso1, P_iso2, sampling=True, size=10):
        """Get marginal logLikelihood via variaitonal posterior
        """
        # TODO: introduce sparse tensor here
        # tensor shapes: (size, n_gene, n_cell)
        _Z = self.Z.sample(size)                         # (size, n_g, n_c)
        _W_cell = self.cell_weight.sample(size)          # (size, n_g, k_c)
        _W_gene = self.gene_weight.sample(size)          # (size, k_g, n_c)
        _zz_var  = tf.sqrt(self.std.sample((size)))      # (size, n_g, 1)
        
        # LogLike part 1: Psi from regression
        _zz_loc = (
            tf.tensordot(_W_cell, Xc, [2, 0]) + 
            tf.transpose(tf.tensordot(Xg, _W_gene, [1, 1]), [1, 0, 2]))
        _normal = tfd.Normal(_zz_loc, _zz_var)
        _logLik_Z = _normal.log_prob(_Z)

        # LogLike part 1: reads counts in regard to Z
        Psi1 = tf.sigmoid(_Z)                           # fraction of isoform 1
        Psi2 = 1 - Psi1                                 # fraction of isoform 2
        Psi1_log = tf.math.log_sigmoid(_Z)
        Psi2_log = tf.math.log_sigmoid(0 - _Z)
                
        # Calculate element wise logLikelihood
        def _re1(x):
            return tf.reshape(x, (1, self.Ng, self.Nc))
        def _re2(x): 
            return tf.reshape(x, (1, self.Ng, 1))
        
        _logLik_S = (
            _re1(Rs['1']) * (tf.math.log(_re2(P_iso1[:, 0])) + Psi1_log) + 
            _re1(Rs['2']) * (tf.math.log(_re2(P_iso2[:, 1])) + Psi2_log) + 
            _re1(Rs['3']) * (tf.math.log(
                _re2(P_iso1[:, 2]) * Psi1 + _re2(P_iso2[:, 2]) * Psi2)))
                
        return tf.reduce_mean(_logLik_S + _logLik_Z, axis=0)

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
