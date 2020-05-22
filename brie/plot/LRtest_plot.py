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

    
def volcano(LRT_res, x="weight", y="pval", highlight_min=4):
    """Volcano plot for p values and weights
    """
    pval_log10 = LRT_res.pval_log10.copy()
    idx = pval_log10 < -highlight_min
    
    pval_log10[pval_log10 < -10] = -10
    plt.scatter(LRT_res.Wc_loc.numpy()[~idx, 0], -pval_log10[~idx],
               color="gray")
    plt.scatter(LRT_res.Wc_loc.numpy()[idx, 0], -pval_log10[idx], 
                color="firebrick")
    plt.xlabel("weight")
    plt.ylabel("-log10(p value)")
