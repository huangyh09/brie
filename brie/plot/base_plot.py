import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import issparse

def loss(losses, last=200):
    plt.figure(figsize=(8, 3.5))
    plt.subplot(1, 2, 1)
    plt.plot(losses)
    plt.xlabel("iterations")
    plt.ylabel("loss")

    plt.subplot(1, 2, 2)
    plt.plot(range(len(losses)-last, len(losses)), losses[-last:])
    plt.xlabel("iterations")
    plt.ylabel("loss")

    plt.tight_layout()
    plt.show()

    
def counts(adata, genes, size='Psi', color=None, gene_key='index',
           layers=['isoform1', 'isoform2'], nrow=None, ncol=None, 
           show_key="index", add_val=None, noise_scale=0.1, **keyargs):
    """Plot the counts for isoforms
    """
    import pandas as pd
    import seaborn as sns
    if noise_scale > 0:
        np.random.seed = 0
    
    # support a single gene input
    if type(genes) == str:
        genes = [genes]
            
    if ncol is None:
        ncol = min(4, len(genes))
    if nrow is None:
        nrow = np.ceil(len(genes) / ncol)
    
    if color is None:
        color_use = None
    else:
        try:
            if len(color) == adata.shape[0]:
                color_use = color
            else:
                color_use = adata.obs[color]
        except ValueError:
            color_use = None
            #print("Not a valid value for color: %s" %(color))
    
    for i in range(len(genes)):
        plt.subplot(nrow, ncol, i + 1)
        if gene_key is None or gene_key == 'index':
            idx = adata.var.index == genes[i]
        else:
            idx = adata.var[gene_key] == genes[i]
            
        adata_use = adata[:, idx]
        
        x_mat = adata_use.layers[layers[0]]
        y_mat = adata_use.layers[layers[1]]
        if issparse(x_mat):
            x_mat = x_mat.toarray()
        if issparse(y_mat):
            y_mat = y_mat.toarray()
        
        df_tmp = pd.DataFrame({
            "x": x_mat[:, 0] + np.random.normal(loc = 0.0, scale = noise_scale, 
                                                size = len(x_mat[:, 0])),
            "y": y_mat[:, 0] + np.random.normal(loc = 0.0, scale = noise_scale, 
                                                size = len(x_mat[:, 0])),
            color: color_use,
            size: adata_use.layers[size][:, 0]})
        
        ax = sns.scatterplot(x="x",  y="y", hue=color, size=size, 
                             data=df_tmp, **keyargs)

        plt.xlabel("n_reads: %s" %(layers[0]))
        plt.ylabel("n_reads: %s" %(layers[1]))
        if show_key is None or show_key == 'index':
            _title = adata_use.var.index[0]
        else:
            _title = adata_use.var[show_key][0]
        if add_val in adata_use.varm:
            _title += "; %s: %s" %(add_val, adata_use.varm[add_val][0, 0])
        plt.title(_title)
        #print(-res_md.pval_log10[sig_idx[i], 0])

    plt.tight_layout()

    
def vi_post(brie_res):
    """Plot variational posterior
    """
    import seaborn as sns
    from scipy.special import logit
    from ..models import LogitNormal
    
    sns.heatmap(brie_res.Psi, cmap="GnBu", annot=False, fmt=".3f", 
                vmin=0, vmax=1, cbar=False, linewidths=0.1, alpha=0.5)

    # sns.heatmap(brie_res.Psi * 0, cmap="Pastel2", annot=False, fmt=".3f", 
    #             vmin=0, vmax=1, cbar=False, linewidths=0.1, alpha=1)

    for i in range(brie_res.Psi.shape[1]):
        for j in range(brie_res.Psi.shape[0]):
            x = np.linspace(0.01, 0.99, 100)
            y = LogitNormal(scale=brie_res.Z_std[j, i], 
                            loc=logit(brie_res.Psi[j, i])).pdf(x)
            y = y / np.max(y) * 0.9
            xx = i + x
            yy = j - y + 1
            plt.plot(xx, yy, c='black', linewidth=1.5)
            plt.plot([i + brie_res.Psi[j, i], i + brie_res.Psi[j, i]],
                     [j + 0.05, j + 1], c='firebrick', linewidth=0.5)
            
