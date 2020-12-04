import numpy as np
from scipy.stats import multinomial
from scipy.special import logit, expit, gammaln

from scipy.special import logit
from scipy.stats import norm, rv_continuous

class LogitNormal(rv_continuous):
    """LogitNormal distribution based on scipy.stats.rv_continuous
    """
    def __init__(self, loc=0, scale=1):
        super().__init__(self)
        self.loc = loc
        self.scale = scale
        
    def _pdf(self, x):
        return norm.pdf(logit(x), loc=self.loc, scale=self.scale)/(x*(1-x))


def BRIE_base_lik(psi, counts, lengths):
    """Base likelihood function of BRIE model
    """
    size_vect = np.array([psi, (1 - psi), 1]) * lengths
    prob_vect = size_vect / np.sum(size_vect)
    
    rv = multinomial(np.sum(counts), prob_vect)
    return rv.pmf(counts)

def get_CI95(Psi, Z_std):
    """Get the 95% confidence interval
    """
    Z = logit(Psi)
    Z_low = Z - 1.96 * Z_std
    Z_high = Z + 1.96 * Z_std
    
    return expit(Z_low), expit(Z_high)


def logbincoeff(n, k, is_sparse=False):
    """
    Ramanujan's approximation of log [n! / (k! (n-k)!)]
    This is mainly for convinience with pen. Please use betaln or gammaln
    """
    if is_sparse:
        RV_sparse = n.copy() * 0
        idx = (k > 0).multiply(k < n)
        n = np.array(n[idx]).reshape(-1)
        k = np.array(k[idx]).reshape(-1)
        
    RV = gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)
    
    if is_sparse:
        RV_sparse[idx] += RV
        RV = RV_sparse
    return RV
