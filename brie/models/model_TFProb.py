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
    def __init__(self, Nc, Ng, Kc, Kg=0, p_ambiguous=None, 
                 intercept=None, sigma=None, name=None):
        self.Nc = Nc
        self.Ng = Ng
        self.Kc = Kc
        self.Kg = Kg
        
        self.Z_loc = tf.Variable(tf.random.normal([Ng, Nc]), name='Z_loc',
            constraint=lambda t: tf.clip_by_value(t, -9, 9))
        self.Z_std = tf.Variable(tf.random.normal([Ng, Nc]), name='Z_var')
        
        self.Wc_loc = tf.Variable(tf.random.normal([Ng, Kc]), name='Wc_loc')
        self.Wg_loc = tf.Variable(tf.random.normal([Kg, Nc]), name='Wg_loc')
        
        if intercept is None:
            self.intercept = tf.Variable(tf.random.normal([Ng, 1]), name='bias',
                constraint=lambda t: tf.clip_by_value(t, -9, 9))
        else:
            _intercept = tf.ones([Ng, 1]) * intercept
            self.intercept = tf.constant(_intercept, name='bias')
            
        if sigma is None:
            self.sigma_log = tf.Variable(tf.ones([Ng, 1]), name='sigma_log')
        else:
            _sigma = tf.ones([Ng, 1]) * sigma
            self.sigma_log = tf.constant(tf.math.log(_sigma), name='sigma_log')
            
        if p_ambiguous is None:
            self.p_ambiguous = np.ones([Ng, 2], dtype=np.float32) * 0.5
        else:
            self.p_ambiguous = p_ambiguous
        self.rho = self.p_ambiguous[:, 0] / self.p_ambiguous.sum(1)

    @property
    def Psi(self):
        """Logistic of Z
        """
        return tf.sigmoid(self.Z_loc)

    @property
    def sigma(self):
        """Standard deviation of predicted Z"""
        return tf.exp(self.sigma_log)

    @property
    def Z(self):
        """Variational posterior for the logit Psi"""
        return tfd.Normal(self.Z_loc, tf.exp(self.Z_std))
    
    @property
    def Z_prior(self):
        """Predicted informative prior for Z"""
        _zz_loc = tf.matmul(self.Wc_loc, self.Xc) + self.intercept
        if self.Kg > 0 and self.Xg is not None:
            _zz_loc += tf.matmul(self.Xg, self.Wg_loc)
        return tfd.Normal(_zz_loc, self.sigma)
    
        
    def logLik_MC(self, count_layers, mode="post", size=10):
        """Get marginal logLikelihood on variational or prior distribution
        with Monte Carlo sampling
        """
        # TODO: introduce sparse tensor here
        
        # Reshape the tensors
        def _re1(x):
            return tf.expand_dims(x, 0) #(1, self.Ng, self.Nc)
        def _re2(x): 
            return tf.expand_dims(tf.expand_dims(x, 0), 2) #(1, self.Ng, 1)
        
        
        ## Manual re-parametrization (VAE) - works similarly well as build-in
        ## https://gregorygundersen.com/blog/2018/04/29/reparameterization/
        # _zzz = tfd.Normal(0, 1).sample((size, self.Ng, self.Nc)) #(size, 1, 1)
        # if mode == "prior":
        #     _Z = _re1(self.sigma) * _zzz + _re1(self.Z_prior.parameters['loc'])
        # else:
        #     _Z = _re1(tf.exp(self.Z_std)) * _zzz + _re1(self.Z_loc)
        
        
        ## Build-in re-parametrized: Gaussian is FULLY_REPARAMETERIZED
        if mode == "prior":
            _Z = self.Z_prior.sample(size)      # (size, n_g, n_c)
        else:
            _Z = self.Z.sample(size)            # (size, n_g, n_c)
        

        ## Calculate element wise logLikelihood
        Psi1 = tf.sigmoid(_Z)                   # fraction of isoform 1
        Psi2 = 1 - Psi1                         # fraction of isoform 2
        Psi1_log = tf.math.log_sigmoid(_Z)
        Psi2_log = tf.math.log_sigmoid(0 - _Z)
        
        _logLik_S = (
            _re1(count_layers[0].transpose()) * Psi1_log + 
            _re1(count_layers[1].transpose()) * Psi2_log)
        
        if len(count_layers) > 2 and np.mean(self.rho == 0.5) < 1:
            _logLik_S += _re1(count_layers[2].transpose()) * tf.math.log(
                _re2(self.rho) * Psi1 + _re2(1 - self.rho) * Psi2)
                
        ## return the mean over the sampling
        if mode == "prior":
            return tfp.math.reduce_logmeanexp(_logLik_S, axis=0)
        else:
            return tf.reduce_mean(_logLik_S, axis=0)
    

    def get_loss(self, count_layers, target="ELBO", axis=None, **kwargs):
        """Loss function per gene
        
        Please be careful: for loss function, you should reduce_sum of each 
        module first then add them up!!! Otherwise, it doesn't work propertly
        by adding modules first and then reduce_sum.
        """
        ## target function
        if target == "marginLik":
            return -tf.reduce_sum(
                self.logLik_MC(count_layers, mode="prior", **kwargs), axis=axis)
        else:
            return (
                tf.reduce_sum(tfd.kl_divergence(self.Z, self.Z_prior), 
                              axis=axis) -
                tf.reduce_sum(self.logLik_MC(count_layers, mode="post", 
                                             **kwargs), axis=axis))

    
    def fit(self, count_layers, Xc, Xg=None, target="ELBO", optimizer=None,
            learn_rate=0.05, min_iter=200, max_iter=5000, add_iter=100, 
            epsilon_conv=1e-2, verbose=True, **kwargs):
        """Fit the model's parameters"""
        start_time = time.time()
        
        self.Xc = Xc
        self.Xg = Xg
        self.target = target
        
        ## target function
        loss_fn = lambda: self.get_loss(count_layers, target, **kwargs)
            
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
            
        
        self.loss_gene = self.get_loss(count_layers, target, axis=1, 
                                       size=1000, **kwargs)
        
        self.losses = losses
        
        if verbose:
            print("BRIE2 model fit with %d steps in %.2f min, loss: %.2f" %(
                n_iter, (time.time() - start_time) / 60, 
                tf.reduce_sum(self.loss_gene)))
            
        return losses
