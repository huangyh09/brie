## BRIE model

import time
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
    def __init__(self, Nc, Ng, Kc, Kg=0, intercept=None, 
                 p_ambiguous=None, name=None):
        self.Nc = Nc
        self.Ng = Ng
        self.Kc = Kc
        self.Kg = Kg
        
        self.Z_loc = tf.Variable(tf.random.normal([Ng, Nc]), name='Z_loc',
            constraint=lambda t: tf.clip_by_value(t, -9, 9))
        self.Z_std = tf.Variable(tf.random.normal([Ng, Nc]), name='Z_var')

        self.sigma_log = tf.Variable(tf.ones([Ng, 1]), name='sigma_log')
        
        self.Wc_loc = tf.Variable(tf.random.normal([Ng, Kc]), name='Wc_loc')
        if self.Kg == 0:
            self.Wg_loc = tf.constant(tf.zeros([Kg, Nc]), name='Wg_loc')
        else:
            self.Wg_loc = tf.Variable(tf.random.normal([Kg, Nc]), name='Wg_loc')
        
        if intercept is None:
            self.intercept = tf.Variable(tf.random.normal([Ng, 1]), name='biase')
        else:
            _intercept = tf.ones([Ng, 1]) * intercept
            self.intercept = tf.constant(_intercept, name='biase')
            
        if p_ambiguous is None:
            self.p_ambiguous = tf.ones([Ng, 2]) * 0.5
        else:
            self.p_ambiguous = p_ambiguous
        self.rho = self.p_ambiguous[:, 0] / self.p_ambiguous.sum(1)

#         if prob_ambiguous is None:
#             self.rho_logit = tf.Variable(tf.random.normal([Ng, 1]), 
#                                          name='rho_logit')
#         else:
#             rho = prob_ambiguous[:, 0] / prob_ambiguous.sum(axis=1)
#             self.rho_logit = tf.constant(tf.math.log(rho / (1 - rho)), 
#                                          name='rho_logit')
           
#     @property
#     def rho(self):
#         """Ratio of ambiguous reads between isoform 1 and 2"""
#         return tf.sigmoid(self.rho_logit)

    @property
    def Z(self):
        """Variational posterior for the logit Psi"""
        return tfd.Normal(self.Z_loc, tf.exp(self.Z_std))
        
    @property
    def sigma(self):
        """Standard deviation of predicted Z"""
        return tf.exp(self.sigma_log)
        
    
    def regression_KL(self, Xc, Xg=None):
        """Get the expectation of z analytically
        """
        _zz_loc = tf.matmul(self.Wc_loc, Xc) + self.intercept
        if self.Kg > 0 and Xg is not None:
            _zz_loc += tf.matmul(Xg, self.Wg_loc)
            
        _normal = tfd.Normal(_zz_loc, self.sigma)
        _regression_KL = tfd.kl_divergence(self.Z, _normal)
        
        return _regression_KL
    
    
    def logLik_VB(self, count_layers, sampling=True, size=10):
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
            return tf.expand_dims(x, 0) #(1, self.Ng, self.Nc)
        def _re2(x): 
            return tf.expand_dims(tf.expand_dims(x, 0), 2) #(1, self.Ng, 1)
        
        _logLik_S = (
            _re1(count_layers['1'].transpose()) * Psi1_log + 
            _re1(count_layers['2'].transpose()) * Psi2_log + 
            _re1(count_layers['3'].transpose()) * tf.math.log(
                _re2(self.rho) * Psi1 + _re2(1 - self.rho) * Psi2))

#         _prob1_log = _re2(tf.math.log(1 - self.p_ambiguous[:, 0]))
#         _prob2_log = _re2(tf.math.log(1 - self.p_ambiguous[:, 1]))
#         _logLik_S = (
#             _re1(count_layers['1']) * (_prob1_log + Psi1_log) + 
#             _re1(count_layers['2']) * (_prob2_log + Psi2_log) + 
#             _re1(count_layers['3']) * tf.math.log(
#                 _re2(self.p_ambiguous[:, 0]) * Psi1 + 
#                 _re2(self.p_ambiguous[:, 1]) * Psi2))
                
        return tf.reduce_mean(_logLik_S, axis=0)
    

    def fit(self, count_layers, Xc, Xg=None, optimizer=None,
            learn_rate=0.05, min_iter=200, max_iter=5000, 
            add_iter=100, epsilon_conv=1e-2, verbose=True, **kwargs):
        """Fit the model's parameters"""
        start_time = time.time()
        
        if optimizer is None:
            optimizer = tf.optimizers.Adam(learning_rate=learn_rate)
            
        loss_fn = lambda: (
            tf.reduce_sum(self.regression_KL(Xc, Xg)) -
            tf.reduce_sum(self.logLik_VB(count_layers, **kwargs)))
        
        losses = tfp.math.minimize(loss_fn, 
                                   num_steps=min_iter, 
                                   optimizer=optimizer)
        
        n_iter = min_iter + 0
        while ((losses[-20:-10].numpy().mean() - losses[-10:].numpy().mean() > 
                epsilon_conv) and  n_iter < max_iter):
            n_iter += add_iter
            losses = tf.concat([
                losses,
                tfp.math.minimize(loss_fn, 
                                  num_steps=add_iter, 
                                  optimizer=optimizer)
            ], axis=0)
            
        if verbose:
            print("BRIE2 model fit with %d steps in %.2f min, loss: %.2f" %(
                n_iter, (time.time() - start_time) / 60, losses[-1]))
        
        self.loss_gene = (
            tf.reduce_sum(self.regression_KL(Xc, Xg), axis=1) -
            tf.reduce_sum(self.logLik_VB(count_layers, **kwargs), axis=1))
        
        self.losses = losses
        return losses
