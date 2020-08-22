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


def volcano(adata, x="cell_coeff", y="pval", index=0, score_red=0.0001, 
            n_anno=10, anno_id='index', log_y=True, adjust=True):
    """Volcano plot for p values and weights
    """
    xval = adata.varm[x][:, index]
    yval = adata.varm[y][:, index]
    idx = yval < score_red
    idx_anno = np.argsort(yval)[:n_anno]
    if log_y:
        yval = -np.log10(yval)
        
    plt.scatter(xval[~idx], yval[~idx], color="gray")
    plt.scatter(xval[idx], yval[idx], color="firebrick")
    
    lable = adata.var.index if anno_id is 'index' else adata.var[anno_id]
    
    texts = []
    for i in idx_anno:
        _label = lable[i]
        _xx = xval[i]
        _yy = yval[i]
        texts.append(plt.text(_xx, _yy, _label, size=8))
    
    if adjust and n_anno > 0:
        from adjustText import adjust_text
        adjust_text(texts, arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    
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
    plt.xlabel("-log10(p), expected")
    plt.ylabel("-log10(p), observed")
    
    
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
