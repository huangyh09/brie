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
        self.tau_s1 = tf.Variable(tf.ones([Ng, 1]) * 3, name='var_s1')
        self.tau_s2 = tf.Variable(tf.ones([Ng, 1]) * 1, name='var_s2')
        
        self.Wc_loc = tf.Variable(tf.random.normal([Ng, Kc]), name='Wc_loc')
        self.Wc_std = tf.Variable(tf.random.normal([Ng, Kc]), name='Wc_std')
        if self.Kg == 0:
            self.Wg_loc = tf.constant(tf.zeros([Kg, Nc]), name='Wg_loc')
            self.Wg_std = tf.constant(tf.zeros([Kg, Nc]), name='Wg_std')
        else:
            self.Wg_loc = tf.Variable(tf.random.normal([Kg, Nc]), name='Wg_loc')
            self.Wg_std = tf.Variable(tf.random.normal([Kg, Nc]), name='Wg_std')

        self.set_prior()

    def set_prior(self, 
            w_prior=tfd.Normal(0, 1), 
            tau_prior=tfd.Gamma(7, 3)):
        """Set priors distributions"""
        self.w_prior = w_prior
        self.tau_prior = tau_prior
    
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
    def tau(self):
        """Variational posterior for the precision of regression residue"""
        return tfd.Gamma(tf.exp(self.tau_s1), tf.exp(self.tau_s2))
    
    @property
    def KLdiverg(self):
        """Sum of KL divergences between posteriors and priors"""
        return (
            tf.reduce_sum(tfd.kl_divergence(self.cell_weight, self.w_prior), axis=1) +
            tf.reduce_sum(tfd.kl_divergence(self.tau, self.tau_prior), axis=1) +
            tf.reduce_sum(tfd.kl_divergence(self.gene_weight, self.w_prior)))
    
    def Expect_Z(self, Xc, Xg=None):
        """Get the expectation of z analytically
        """
        # # Monte Carlo Based expectation
        # _Z = self.Z.sample(size)                         # (size, n_g, n_c)
        # _W_cell = self.cell_weight.sample(size)          # (size, n_g, k_c)
        # _W_gene = self.gene_weight.sample(size)          # (size, k_g, n_c)
        # _zz_tau  = tf.sqrt(self.tau.sample((size)))      # (size, n_g, 1)
        # _zz_loc = (
        #     tf.tensordot(_W_cell, Xc, [2, 0]) + 
        #     tf.transpose(tf.tensordot(Xg, _W_gene, [1, 1]), [1, 0, 2]))
        # _normal = tfd.Normal(_zz_loc, 1 / _zz_tau)
        # _logLik_Z = tf.reduce_mean(_normal.log_prob(_Z), axis=0)
        
        _zz_loc = tf.matmul(self.Wc_loc, Xc)
        if self.Kg > 0 and Xg is not None:
            _zz_loc += tf.matmul(Xg, self.Wg_loc)
        _zz_tau = self.tau.mean()
        _normal = tfd.Normal(_zz_loc, 1 / tf.sqrt(_zz_tau))
        # _logLik_Z1 = - self.Z.cross_entropy(_normal)
        _logLik_Z1 = tfd.kl_divergence(self.Z, _normal)
        
        # LogLik part 1.2: analytical Other terms
#         _logLik_Z2 = - 0.5 * (
#             tf.math.digamma(tf.exp(self.tau_s1)) * self.Nc -
#             tf.math.log(tf.exp(self.tau_s1)) * 2  * self.Nc + 
#             tf.math.log(tf.exp(self.tau_s2)) * self.Nc - 
#             tf.matmul(tf.exp(self.Wc_std)**2, Xc**2) * _zz_tau**2)

        _logLik_Z2 = - 0.5 * (
            tf.math.digamma(tf.exp(self.tau_s1)) * self.Nc -
            tf.math.log(tf.exp(self.tau_s1)) * self.Nc -
            tf.matmul(tf.exp(self.Wc_std)**2, Xc**2) * _zz_tau)
        
        if self.Kg > 0 and Xg is not None:
            _logLik_Z2 += - 0.5 * (
                -tf.matmul(Xg**2, tf.exp(self.Wg_std)**2) * _zz_tau)
            # _logLik_Z2 += 0.5 * tf.matmul(Xg**2, tf.exp(self.Wg_std)**2)
            
        return _logLik_Z1 + _logLik_Z2
    
    
    def logLik(self, Xc, Rs, P_iso1, P_iso2, Xg=None, sampling=True, size=10):
        """Get marginal logLikelihood on link distribution
        """
        # TODO: introduce sparse tensor here
        # LogLike part 2: reads counts in regard to Z
        _Z = self.Z.sample(size)                # (size, n_g, n_c)
        Psi1 = tf.sigmoid(_Z)                   # fraction of isoform 1
        Psi2 = 1 - Psi1                         # fraction of isoform 2
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
                
        return tf.reduce_mean(_logLik_S, axis=0)
    

    def fit(self, Xc, Rs, P_iso1, P_iso2, Xg=None, num_steps=100, optimizer=None,
            learn_rate=0.05, **kwargs):
        """Fit the model's parameters"""
        if optimizer is None:
            optimizer = tf.optimizers.Adam(learning_rate=learn_rate)
            
        loss_fn = lambda: (
            #tf.reduce_sum(self.KLdiverg) + 
            tf.reduce_sum(self.Expect_Z(Xc, Xg)) -
            tf.reduce_sum(self.logLik(Xc, Rs, P_iso1, P_iso2,
                                      Xg, **kwargs)))
        
        losses = tfp.math.minimize(loss_fn, 
                                   num_steps=num_steps, 
                                   optimizer=optimizer)
        
        self.loss_gene = (
            #tf.reduce_sum(self.KLdiverg) + 
            tf.reduce_sum(self.Expect_Z(Xc, Xg), axis=1) -
            tf.reduce_sum(self.logLik(Xc, Rs, P_iso1, P_iso2,
                                      Xg, **kwargs), axis=1))
        return losses
