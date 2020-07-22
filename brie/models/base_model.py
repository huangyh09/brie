import numpy as np
from scipy.stats import multinomial

def BRIE_base_lik(psi, counts, lengths):
    """Base likelihood function of BRIE model
    """
    size_vect = np.array([psi, (1 - psi), 1]) * lengths
    prob_vect = size_vect / np.sum(size_vect)
    
    rv = multinomial(np.sum(counts), prob_vect)
    return rv.pmf(counts)