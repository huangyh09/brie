## simulator for BRIE2

import numpy as np
import tensorflow as tf
import tensorflow_probability as tfp
from tensorflow_probability import distributions as tfd

def simulator(model, adata):
    """Simulate read counts for BRIE model
    """
    # generating Psi
#     _Psi = tf.sigmoid(model.Z_prior.sample()).numpy()
#     _Psi_tensor = np.expand_dims(_Psi, 2), 
#     Prob_iso = np.zeros((n, 3, 2))

    adata = adata.copy()
    total_counts = (adata.layers['isoform1'] + 
                    adata.layers['isoform2'] + 
                    adata.layers['ambiguous']).transpose()
        
    _Psi = tf.sigmoid(model.Z_prior.sample()).numpy()
    _binom = tfd.Binomial(total_counts, probs=_Psi) ## accounting for tranL
    N1 = _binom.sample().numpy()
    N2 = total_counts - N1
    
    _binom = tfd.Binomial(N1, probs=model.p_ambiguous[:, :1])
    N13 = _binom.sample().numpy()
    N11 = N1 - N13
    
    _binom = tfd.Binomial(N2, probs=model.p_ambiguous[:, 1:])
    N23 = _binom.sample().numpy()
    N22 = N2 - N23
    
    adata.layers['isoform1'] = N11.transpose()
    adata.layers['isoform2'] = N22.transpose()
    adata.layers['ambiguous'] = (N13 + N23).transpose()
        
    # update adata
    if model.Xc is not None and model.Xc.shape[0] > 0:
        adata.obsm['Xc_sim'] = model.Xc.transpose()
        adata.varm['cell_coeff_sim'] = model.Wc_loc.numpy()
        
    if model.Xg is not None and model.Xg.shape[1] > 0:
        adata.varm['Xg_sim'] = Xg
        adata.obsm['gene_coeff_sim'] = model.Wg_loc.numpy().transpose()
        
    adata.varm['intercept_sim'] = model.intercept.numpy()
    adata.varm['sigma_sim'] = model.sigma.numpy()
        
    adata.layers['Psi_sim'] = _Psi.transpose()
    
    return adata
    