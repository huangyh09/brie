## BRIE-Beta model

import time
import numpy as np
from scipy. special import betaln
from scipy.sparse import issparse

import tensorflow as tf
import tensorflow_probability as tfp
from tensorflow.math import digamma, polygamma
from tensorflow_probability import distributions as tfd


# def Binomial_constant_ln(X_a, X_b, is_sparse):
#     """Calculate the constant difference from Binomial to Beta
#     """
#     if is_sparse:
#         RV = X_a.copy() * 0
#         idx = (X_a > 0).multiply(X_b > 0)
#         _A = np.array(X_a[idx]).reshape(-1)
#         _B = np.array(X_b[idx]).reshape(-1)
        
#         RV[idx] = (gammaln(_A + _B + 1) - gammaln(_A + 1) - gammaln(_B + 1) - 
#                    betaln(_A + 1, _B + 1))
#     else:
#         RV = (gammaln(X_a + X_b + 1) - gammaln(X_a + 1) - gammaln(X_b + 1) - 
#               betaln(X_a + 1, X_b + 1))
        
#     return RV


def entropy_Beta_LogitNormal(Z_a, Z_b, Y_mu, Y_std):
    """Calculate KL divergence between Beta distribution and LogitNormal
    
    https://en.wikipedia.org/wiki/Beta_distribution#Other_moments
    https://en.wikipedia.org/wiki/Logit-normal_distribution
    """    
    E_logit = digamma(Z_a) - digamma(Z_b)
    E_logit_sqr = (
        digamma(Z_a)**2 + polygamma(tf.constant(3, tf.float32), Z_a) + 
        digamma(Z_b)**2 + polygamma(tf.constant(3, tf.float32), Z_b) -
        digamma(Z_a) * digamma(Z_b) * 2
    )
    
    _RV_part1 = - tf.constant(0.5*np.log(2 * np.pi), tf.float32) - tf.math.log(Y_std)
    _RV_part2 = - digamma(Z_a) - digamma(Z_b) + 2 * digamma(Z_a + Z_b)
    _RV_part3 = - (E_logit_sqr - 2 * Y_mu * E_logit + Y_mu **2) / (2 * Y_std**2)
    
    return _RV_part1 + _RV_part2 + _RV_part3


def KL_Beta_Binomial(Z_a, Z_b, X_a, X_b):
    """Calculate KL divergence between Beta distribution and Binomial likelihood
    See the relationship between Beta function and binomial coefficient:
    https://en.wikipedia.org/wiki/Beta_function#Properties
    """
    # TODO: introduce sparse matrix
    _KL = tfd.kl_divergence(tfd.Beta(Z_a, Z_b), tfd.Beta(X_a + 1, X_b + 1))
    _diff_binomLik_to_beta = -tf.math.log(X_a + X_b + 1)
    return _KL + _diff_binomLik_to_beta


class BRIE2_Beta():
    """
    Ng : number of genes
    Nc : number of cells
    Kg : number of gene features
    Kc : number of cell features
    
    Note, as the prior is LogitNormal, but posterior is Beta. The KL can never 
    be zero, so may introduce unwanted non-zero values for empty genes.
    """
    def __init__(self, Nc, Ng, Kc=0, Kg=0, effLen=None,
                 intercept=None, intercept_mode='gene',
                 sigma=None, name=None):
        self.Nc = Nc
        self.Ng = Ng
        self.Kc = Kc
        self.Kg = Kg
        self.effLen = effLen # (Ng, 3 * 2)
        self.intercept_mode = intercept_mode

        self.Z_a_log = tf.Variable(tf.random.uniform([Nc, Ng]), name='Z_alpha')
        self.Z_b_log = tf.Variable(tf.random.uniform([Nc, Ng]), name='Z_beta')
        
        self.Wc_loc = tf.Variable(tf.random.normal([Kc, Ng]), name='Wc_loc')
        self.Wg_loc = tf.Variable(tf.random.normal([Nc, Kg]), name='Wg_loc')
        
        if intercept_mode.upper() == 'GENE':
            _intercept_shape = (1, Ng)
        elif intercept_mode.upper() == 'CELL':
            _intercept_shape = (Nc, 1)
        else:
            # print("[BIRE2] Error: intercept_mode only supports gene or cell")
            _intercept_shape = (1, Ng)
            
        if intercept is None:
            self.intercept = tf.Variable(tf.random.normal(_intercept_shape), 
                name='bias', constraint=lambda t: tf.clip_by_value(t, -9, 9))
        else:
            _intercept = tf.ones(_intercept_shape) * intercept
            self.intercept = tf.constant(_intercept, name='bias')
            
        if sigma is None:
            self.sigma_log = tf.Variable(tf.ones([1, Ng]), name='sigma_log')
        else:
            _sigma = tf.ones([1, Ng]) * sigma
            self.sigma_log = tf.constant(tf.math.log(_sigma), name='sigma_log')

    @property
    def Z_a(self):
        return tf.math.exp(self.Z_a_log)
    
    @property
    def Z_b(self):
        return tf.math.exp(self.Z_b_log)
    
    @property
    def Z_std(self):
        return 1 / (self.Z_a + self.Z_b)
    
    @property
    def Psi(self):
        """Mean value of Psi in variational posterior"""
        return self.Z_a / (self.Z_a + self.Z_b)
    
    @property
    def PsiDist(self):
        """Variational Beta distribution for Psi"""
        return tfd.Beta(self.Z_a, self.Z_b)
    
    @property
    def Psi95CI(self):
        """95% confidence interval around mean"""
        #return self.PsiDist.quantile(0.975) - self.PsiDist.quantile(0.025)
        
        from scipy.stats import beta
        return (beta.ppf(0.975, self.Z_a.numpy(), self.Z_b.numpy()) - 
                beta.ppf(0.025, self.Z_a.numpy(), self.Z_b.numpy()))
    
    @property
    def sigma(self):
        """Standard deviation of predicted Z"""
        return tf.exp(self.sigma_log)
    
    @property
    def Z_prior(self):
        """Predicted informative prior for Z"""
        _zz_loc = tf.zeros((self.Nc, self.Ng))
        if self.Kc > 0 and self.Xc is not None:
            _zz_loc = tf.matmul(self.Xc, self.Wc_loc) 
        if self.Kg > 0 and self.Xg is not None:
            _zz_loc += tf.matmul(self.Wg_loc, self.Xg.T)
        _zz_loc += self.intercept
        return _zz_loc
    

    def get_loss(self, count_layers, target="ELBO", axis=None):
        """Loss function per gene (axis=0) or all genes
        
        Please be careful: for loss function, you should reduce_sum of each 
        module first then add them up!!! Otherwise, it doesn't work propertly
        by adding modules first and then reduce_sum.
        """
        for i in range(len(count_layers)):
            if issparse(count_layers[i]):
                count_layers[i] = count_layers[i].toarray()
        
        ## target function
        loss = (
            tf.reduce_sum(KL_Beta_Binomial(self.Z_a, self.Z_b,
                                           count_layers[0], count_layers[1]),
                          axis=axis) -
            tf.reduce_sum(entropy_Beta_LogitNormal(self.Z_a, self.Z_b, 
                                              self.Z_prior, self.sigma), 
                          axis=axis)
        )
        return loss

    
    def fit(self, count_layers, Xc=None, Xg=None, target="ELBO", optimizer=None, 
            learn_rate=0.05, min_iter=200, max_iter=5000, add_iter=100, 
            epsilon_conv=1e-2, verbose=True):
        """Fit the model's parameters"""
        start_time = time.time()
        
        self.Xc = Xc  #(Nc, Kc)
        self.Xg = Xg  #(Ng, Kg)
        self.target = target
        
        ## target function
        loss_fn = lambda: self.get_loss(count_layers, target)
            
        ## optimization
        if optimizer is None:
            optimizer = tf.optimizers.Adam(learning_rate=learn_rate)
        
        losses = tfp.math.minimize(loss_fn, 
                                   num_steps=min_iter, 
                                   optimizer=optimizer)
        
        n_iter = min_iter + 0
        d1 = min(50, add_iter / 2)
        d2 = d1 * 2
        while ((losses[-d2:-d1].numpy().mean() - losses[-d1:].numpy().mean() > 
                epsilon_conv) and  n_iter < max_iter):
            n_iter += add_iter
            losses = tf.concat([
                losses,
                tfp.math.minimize(loss_fn, 
                                  num_steps=add_iter, 
                                  optimizer=optimizer)
            ], axis=0)
            
        
        self.loss_gene = self.get_loss(count_layers, target, axis=0).numpy()
        
        self.losses = losses
        
        if verbose:
            print("[BRIE2] model fit with %d steps in %.2f min, loss: %.2f" %(
                n_iter, (time.time() - start_time) / 60, 
                tf.reduce_sum(self.loss_gene)))
            
        return losses
