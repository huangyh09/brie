## simulator for BRIE2

import numpy as np
from scipy.special import expit
from tensorflow_probability import distributions as tfd

def simulator(adata, Psi=None, effLen=None, mode="posterior",
              layer_keys=['isoform1', 'isoform2', 'ambiguous'],
              prior_sigma=None):
    """Simulate read counts for BRIE model
    """
    # Check Psi
    if Psi is None and "Psi" not in adata.layers:
        print("Error: no Psi available in adata.layers.")
        exit()
    elif Psi is None:
        if mode == "posterior":
            Psi = adata.layers['Psi'].copy()
        else:
            Psi = np.zeros((adata.shape), np.float32)
            if 'Xc' in adata.obsm and adata.obsm['Xc'].shape[1] > 0:
                Psi += np.dot(adata.obsm['Xc'], adata.varm['cell_coeff'].T)
            if 'Xg' in adata.varm and adata.varm['Xg'].shape[1] > 0:
                Psi += np.dot(adata.obsm['gene_coeff'], adata.varm['Xg'].T)
            if 'intercept' in adata.varm and adata.varm['intercept'].shape[1] > 0:
                Psi += adata.varm['intercept'].T
            if 'intercept' in adata.obsm and adata.obsm['intercept'].shape[1] > 0:
                Psi += adata.obsm['intercept']
                
            adata.layers['Psi_sim_noNoise'] = expit(Psi)
            
            if prior_sigma is None:
                _sigma = adata.varm['sigma'].T
            else:
                _sigma = np.ones([1, adata.shape[1]]) * prior_sigma
            _noise = np.random.normal(loc=0.0, scale=_sigma, size=None) 
            Psi += _noise
            
            Psi[Psi > 9] = 9
            Psi[Psi < -9] = -9
            Psi = expit(Psi)
    adata.layers['Psi_sim'] = Psi

    # Check effective length for isoform specific positions
    if effLen is None and 'effLen' not in adata.varm:
        print("Error: no effLen available in adata.varm.")
        exit()
    elif effLen is None:
        effLen = adata.varm['effLen'][:, [0, 4, 5]]
    else:
        effLen = effLen[:, [0, 4, 5]].copy()
    effLen = np.expand_dims(effLen, 0)
    
    Psi_tensor = np.concatenate((
        np.expand_dims(Psi, 2), 
        1-np.expand_dims(Psi, 2), 
        np.ones((Psi.shape[0], Psi.shape[1], 1), np.float32)
    ), axis=2)
    
    Phi = Psi_tensor * effLen
    Phi = Phi / np.sum(Phi, axis=2, keepdims=True)

    adata = adata.copy()
    total_counts = np.zeros(adata.shape, np.float32)
    for i in range(len(layer_keys)):
        total_counts += adata.layers[layer_keys[i]]
        
    model = tfd.Multinomial(total_counts, probs=Phi)
    count_sim = model.sample().numpy()
        
    adata.layers[layer_keys[0]] = count_sim[:, :, 0]
    adata.layers[layer_keys[1]] = count_sim[:, :, 1]
    adata.layers[layer_keys[2]] = count_sim[:, :, 2]
    
    return adata
    