import numpy as np
from scipy.stats import multinomial
from scipy.special import logit, expit

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