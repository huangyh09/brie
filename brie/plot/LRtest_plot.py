import numpy as np
import matplotlib.pyplot as plt

def LRT_basic(LRT_res):
    fig = plt.figure(figsize=(10, 3))
    plt.subplot(1, 3, 1)
    plt.hist(LRT_res.Wc_loc.numpy()[:, 0], 30)
    plt.xlabel("weight")
    plt.ylabel("#genes")

    plt.subplot(1, 3, 2)
    plt.hist(LRT_res.LR[:, 0], 30)
    plt.xlabel("log likelihood ratio")

    plt.subplot(1, 3, 3)
    plt.hist(-LRT_res.pval_log[:, 0] * np.log10(np.exp(1)), 30)
    plt.xlabel("-log10(p value)")

    plt.tight_layout()
    plt.show()


def volcano(adata, index, x="cell_coeff", y="pval", pval_red=0.0001):
    """Volcano plot for p values and weights
    """
    yval = adata.varm[y][:, index]
    idx = yval < pval_red
    
    #pval_log10[pval_log10 < -10] = -10
    plt.scatter(adata.varm[x][:, index][~idx], 
                -np.log10(adata.varm[y][:, index][~idx]),
               color="gray")
    plt.scatter(adata.varm[x][:, index][idx], 
                -np.log10(adata.varm[y][:, index][idx]), 
                color="firebrick")
    plt.xlabel(x)
    plt.ylabel("-log10(%s)" %(y))
    
    
def qqplot(pval):
    """QQ plot between expected uniformed p values (NULL) and observed p values 
    """
    pval_obs = np.sort(pval.reshape(-1))
    pval_exp = np.linspace(0, 1, len(pval_obs) + 2)[1:-1]
    plt.plot(-np.log10(pval_exp), -np.log10(pval_exp), color="darkgrey")
    plt.scatter(-np.log10(pval_exp), -np.log10(pval_obs), 
                facecolors='none', edgecolors='dimgrey')
    plt.xlabel("expected -log10(p)")
    plt.ylabel("observed -log10(p)")
    
    
# def power_plot(score, group, threshold=0.05):
#     """Plot the power for each categories
#     """
    
#     prop_detect = np.zeros(len(group))
#     for i in range(len(group)):
#         idx = adata_sim_pow.varm['cell_coeff_sim'][:, 0] == group[i]
#         prop_detect[i] = np.mean(score[idx] < 0.05)
        
#     plt.bar(group.astype(str), prop_detect, width=0.5)
#     plt.xlabel("abs(effect size)")
#     plt.ylabel("power: fdr < 0.05")
#     plt.show()
    
    
# def volcano(LRT_res, x="weight", y="pval", highlight_min=4):
#     """Volcano plot for p values and weights
#     """
#     pval_log10 = LRT_res.pval_log10.copy()
#     idx = pval_log10 < -highlight_min
    
#     pval_log10[pval_log10 < -10] = -10
#     plt.scatter(LRT_res.Wc_loc.numpy()[~idx, 0], -pval_log10[~idx],
#                color="gray")
#     plt.scatter(LRT_res.Wc_loc.numpy()[idx, 0], -pval_log10[idx], 
#                 color="firebrick")
#     plt.xlabel("weight")
#     plt.ylabel("-log10(p value)")
