
import sys
import time
import numpy as np
import multiprocessing


class BayesianRegress:
    def __init__(self, lambda_=0.1, sigma_=2.5, coef_=None, intercept_=None,
                 fitsigma=True, fitlambda=False):
        """Initialize a Bayesian regressor.
        _lambda: L2 constrains of weights
        _sigma: square root of variance.
        """
        self.sigma_     = sigma_
        self.lambda_    = lambda_
        self.fitsigma   = fitsigma
        self.fitlambda  = fitlambda
        self.coef_      = coef_
        self.intercept_ = intercept_
        
    def fit(self, X, Y):
        """Fit the model.
        """
        X = np.array(X)
        Y = np.array(Y)
        if len(X.shape) == 1:
            X = X.reshape(-1,1)
        self.X = np.append(X, np.ones((X.shape[0],1)), axis=1)
        
        self.updateW(Y=Y)
        if self.fitsigma is True:
            self.sigma_ = np.std(Y-self.predict(X))
            self.updateW(Y=Y)
            
    def updateW(self, Y, non_W=False):
        """
        Update weights with Y
        This is useful for fixed parameters.
        If non_W is True, don't update self.W_prefix.
        """
        #TODO: check whether updating sigma_
        if non_W is False:
            self.W_matrix = self.lambda_ * np.identity(self.X.shape[1])
            self.W_prefix = np.dot(np.linalg.inv(np.dot(self.X.T, self.X)
                                + self.W_matrix*self.sigma_**2), self.X.T)
        
        weights = np.dot(self.W_prefix, Y)
        self.coef_ = weights[:-1]
        self.intercept_ = weights[-1]
        
    def predict(self, X):
        """Predict Y"""
        if len(X.shape) == 1:
            X = X.reshape(-1,1)
        RV = np.dot(self.coef_, X.T) + self.intercept_
        return RV.reshape(-1)


def normal_pdf(x, mu, cov, log=True):
    """
    Calculate the probability density of Gaussian (Normal) distribution.

    Parameters
    ----------
    x : float or 1-D array_like, (K, )
        The variable for calculating the probability density.
    mu : float or 1-D array_like, (K, )
        The mean of the Gaussian distribution.
    cov : float or 2-D array_like, (K, K)
        The variance or the covariance matrix of the Gaussian distribution.
    log : bool
        If true, the return value is at log scale.

    Returns
    -------
    pdf : numpy float
        The probability density of x.
    """
    x = np.array(x) - np.array(mu)
    if len(np.array(cov).shape) < 2:
        cov = np.array(cov).reshape(-1,1)
    if len(np.array(x).shape) < 1:
        x = np.array(x).reshape(-1)
    cov_inv = np.linalg.inv(cov)
    cov_det = np.linalg.det(cov)
    if cov_det < 0:
        print("The det of covariance is negative, please check!")
        return None
    pdf = (-0.5*np.log(2*np.pi*cov_det) - 0.5*np.dot(np.dot(x, cov_inv), x))
    if log == False: pdf = np.exp(pdf)
    return pdf


def Geweke_Z(X, first=0.1, last=0.5):
    """
    Geweke diagnostics for MCMC chain convergence.
    See Geweke J. Evaluating the accuracy of sampling-based approaches to the 
    calculation of posterior moments[M]. Minneapolis, MN, USA: Federal Reserve 
    Bank of Minneapolis, Research Department, 1991.
    and https://pymc-devs.github.io/pymc/modelchecking.html#formal-methods

    Parameters
    ----------
    X : 1-D array_like, (N, )
        The uni-variate MCMC sampled chain for convergence diagnostic.
    first : float
        The proportion of first part in Geweke diagnostics.
    last : float
        The proportion of last part in Geweke diagnostics.

    Returns
    -------
    Z : float
        The Z score of Geweke diagnostics.
    """
    N = X.shape[0]
    A = X[:int(first*N)]
    B = X[int(last*N):]
    if np.sqrt(np.var(A) + np.var(B)) == 0: 
        Z = None
    else:
        Z = abs(A.mean() - B.mean()) / np.sqrt(np.var(A) + np.var(B))
    return Z


def Iso_read_check(R_mat, len_isos, prob_isos):
    """
    Check the input data for isoform quantification.

    Parameters
    ----------
    R_mat : 2-D array_like, (N, K)
        N reads identities of belonging to K isoforms.
    prob_isos : 2-D array_like, (N, K)
        N reads probablity of belonging to K isoforms.
    len_isos : 1-D array_like, (K,)
        The effective length of each isoform.

    Returns
    -------
    prob_isos : 2-D array_like, (N, K)
        N reads probablity of belonging to K isoforms.
    len_isos : 1-D array_like, (K,)
        The effective length of each isoform.
    """
    idx = (len_isos != len_isos)
    len_isos[idx] = 0.0
    prob_isos[:,idx] = 0.0
    R_mat[:,idx] = False

    idx = np.where(R_mat != R_mat)
    R_mat[idx] = False

    idx = np.where(prob_isos != prob_isos)
    prob_isos[idx] = 0.0

    idx = (R_mat.sum(axis=1) > 0) * (prob_isos.sum(axis=1) > 0)
    R_mat = R_mat[idx,:]
    prob_isos = prob_isos[idx,:]

    return R_mat, prob_isos, len_isos


def MH_propose(Y_now, Y_cov, prob_isos, len_isos, gene_Cnt=None, 
    total_count=10**6, F_pre=None, F_sigma=None, M=1, ftype="RPK"):
    """
    A Matroplis-Hasting sampler with multivariate Gaussian proposal.
    """
    cnt = 0.0
    Y_try = np.zeros(Y_now.shape[0])
    Y_all = np.zeros((Y_now.shape[0], M))
    Psi_all = np.zeros((Y_now.shape[0], M))
    Cnt_all = np.zeros((Y_now.shape[0], M))

    Psi_now = np.exp(Y_now) / np.sum(np.exp(Y_now))
    Fsi_now = len_isos*Psi_now / np.sum(len_isos*Psi_now)
    Cnt_now = gene_Cnt*Fsi_now

    if ["RPK", "RPKM", "FPKM", "rpk", "rpkm", "fpkm"].count(ftype) == 1:
        # F_now = np.log2(Cnt_now* 10**9 / total_count +0.00001)
        F_now = Cnt_now / len_isos / total_count * 10**9
        F_now = np.log10(F_now +0.01)
    elif ftype == "Y" or ftype == "y":
        F_now = Y_now
    else:
        F_now = Psi_now

    P_now = np.log(np.dot(prob_isos, Fsi_now)).sum()
    for k in range(F_now.shape[0]):
        if F_pre[k] is None or  F_pre[k] != F_pre[k]:
            continue
        P_now += normal_pdf(F_now[k], F_pre[k], F_sigma**2)

    for m in range(M):
        # propose
        Y_try[:-1] = np.random.multivariate_normal(Y_now[:-1], Y_cov)
        Y_try[Y_try < -700] = -700
        Y_try[Y_try > 700 ] = 700
        Q_now = normal_pdf(Y_now[:-1], Y_try[:-1], Y_cov)
        Q_try = normal_pdf(Y_try[:-1], Y_now[:-1], Y_cov)

        Psi_try = np.exp(Y_try) / np.sum(np.exp(Y_try))
        Fsi_try = len_isos*Psi_try / np.sum(len_isos*Psi_try)
        Cnt_try = gene_Cnt*Fsi_try

        if ["RPK", "RPKM", "FPKM", "rpk", "rpkm", "fpkm"].count(ftype) == 1:
            F_try = Cnt_try / len_isos / total_count * 10**9
            F_try = np.log10(F_try +0.01)
        elif ftype == "Y" or ftype == "y":
            F_try = Y_try
        else:
            F_try = Psi_try

        _prob = np.dot(prob_isos, Fsi_try)
        if len(np.where(_prob <= 0)[0]) >= 1: 
            P_try = -np.inf
        else:
            P_try = np.log(_prob).sum()
            for k in range(F_try.shape[0]):
                if F_pre[k] is None or  F_pre[k] != F_pre[k]: 
                    continue
                P_try += normal_pdf(F_try[k], F_pre[k], F_sigma**2)

        # step 2: accept or reject the proposal
        alpha = np.exp(min(P_try+Q_now-P_now-Q_try, 0))
        if alpha is None:
            print("alpha is none!")
        elif np.random.rand(1) < alpha:
            #print alpha
            cnt += 1
            P_now = P_try + 0.0
            Y_now = Y_try + 0.0
            Psi_now = Psi_try + 0.0
            Cnt_now = Cnt_try + 0.0
        
        Y_all[:, m] = Y_now
        Psi_all[:, m] = Psi_now
        Cnt_all[:, m] = Cnt_now

    # print("Accept ratio: %.3f for %d isoforms with %d reads." 
    #     %(cnt/(m+1.0), prob_isos.shape[1], prob_isos.shape[0]))
    return Y_all, Psi_all, Cnt_all


def brie_MH_Heuristic(R_mat, len_isos, prob_isos, feature_all, idxF, 
    weights_in=None, _sigma=None, _lambda=2.4, ftype="Y", 
    total_count=10**6, M=10000, Mmin=1000, gap=10, nproc=1): # _lambda=0.134
    """
    Matroplis-Hasting sampler to infer Brie parameters in a heuristic way,
    namely, sampling isoform proportion and optimize the bayesian regression
    in turn.

    Parameters
    ----------
    R_mat: a list of array of bool
        each element in the list is a gene matrix (N*m) of reads identity
    len_isos: a list of array of int
        each element in the list is an array for length of m isoforms
    prob_isos: a list of array of float
        each element in the list is a gene matrix (N*m) for isoforms specific
        probability
    feature_all: an array 
        size (tranNum * K) for multiple isoforms or (tranNum/2, K) splicing 
        events
    weights_in: array, (K+1)
        input weights for Bayesian regression, which may be learned from 
        better data sets.
    _sigma: float
        the variance of residue between target and fitted prior
    _lambda: float
        the L2 constrain of Baysian regression
    ftype: string
        the type of target for regression: FPKM, Y or Psi. Psi=softmax(Y)
    total_count: int
        the total number of reads in the whole bam file(s)
    M: int
        the maximum iterations
    Mmin: int
        the minmum iterations
    gap: int
        the gap iterations between updates of weights in regression
    nproc: int
        the number of processors to use

    Returns
    -------
    Psi_all: samplers of Psi
    Y_all: samplers of Y
    FPKM_all: samplers of FPKM
    Cnt_all: samplers if reads counts
    W_all: sampled weights of Bayesian regression
    """
    START_TIME = time.time()

    tranLen = np.array([])
    geneNum = len(len_isos)
    for t in range(len(len_isos)):
        R_mat[t], prob_isos[t], len_isos[t] = Iso_read_check(R_mat[t], 
            len_isos[t], prob_isos[t])
        prob_isos[t] = R_mat[t] * prob_isos[t]
        tranLen = np.append(tranLen, len_isos[t])
    tranNum = len(tranLen)

    if _sigma is None or _sigma != _sigma:
        sigma_in = 1.5
    else:
        sigma_in = _sigma

    X = feature_all[idxF, :]
    W_mat = _lambda*np.identity(X.shape[1])
    W_pt1 = np.dot(np.linalg.inv(np.dot(X.T, X) + W_mat*sigma_in**2), X.T)

    F_pre = np.zeros(tranNum)
    Y_now = np.zeros(tranNum)
    Y_all = np.zeros((tranNum, M))
    W_all = np.zeros((X.shape[1], M/gap))
    Psi_now = np.zeros(tranNum)
    Psi_all = np.zeros((tranNum, M))
    Cnt_now = np.zeros(tranNum)
    Cnt_all = np.zeros((tranNum, M))
    gCounts = np.zeros(geneNum)
    Idx_all = np.zeros(geneNum+1, "int")
    for g in range(len(len_isos)):
        Idx_all[g+1] = Idx_all[g] + len(len_isos[g])
        idxG = range(Idx_all[g], Idx_all[g+1])

        _psi = np.exp(Y_now[idxG]) / np.sum(np.exp(Y_now[idxG]))
        _fsi = (len_isos[g]*_psi) / np.sum(len_isos[g]*_psi)

        gCounts[g] = prob_isos[g].shape[0]
        Psi_now[idxG] = _psi
        Cnt_now[idxG] = _fsi * gCounts[g]

    if ["RPK", "RPKM", "FPKM", "rpk", "rpkm", "fpkm"].count(ftype) == 1:
        F_now = Cnt_now / tranLen / total_count * 10**9
        F_now = np.log10(F_now +0.01)
    elif ftype == "Y" or ftype == "y":
        F_now = Y_now
    else:
        F_now = Psi_now
    W_sub = np.dot(W_pt1, F_now[idxF]) # suboptimal
    if weights_in is not None:
        W_sub = weights_in
    F_pre[:] = None
    F_pre[idxF] = np.dot(W_sub, X.T)

    # nproc = 1 #!!!
    CONVERG = np.zeros(tranNum, "bool")
    for m in range(M/gap):
        idxT = range(m*gap, (m+1)*gap)
        # step 1: propose a value (original)
        if nproc == 1:
            for g in range(len(len_isos)):
                idxG = range(Idx_all[g], Idx_all[g+1])
                # adaptive MCMC
                if m*gap >= 11:
                    Y_cov = np.cov(Y_all[idxG[:-1], :m*gap])
                else:
                    Y_cov = 1.5 * np.diag(np.ones(len(idxG)-1))
                Y_cov  = Y_cov + np.diag(np.ones(len(idxG)-1))*0.001
                Y_cov *= 5.0/(len(idxG)-1)/(1+prob_isos[g].shape[0]/5000.0)

                _Y, _Psi, _Cnt = MH_propose(Y_now[idxG], Y_cov, prob_isos[g], 
                    len_isos[g], gCounts[g], total_count, F_pre[idxG], 
                    sigma_in, gap, ftype)
                for k in idxG:
                    Y_all[k, idxT] = _Y[k-idxG[0], :]
                    Psi_all[k, idxT] = _Psi[k-idxG[0], :]
                    Cnt_all[k, idxT] = _Cnt[k-idxG[0], :]

        # step 1: propose a value (parallise)
        else:
            pool = multiprocessing.Pool(processes=nproc)
            result = []
            for g in range(len(len_isos)):
                idxG = range(Idx_all[g], Idx_all[g+1])
                # adaptive MCMC
                if m*gap >= 11:
                    Y_cov = np.cov(Y_all[idxG[:-1], :m*gap])
                else:
                    Y_cov = 1.5 * np.diag(np.ones(len(idxG)-1))
                Y_cov  = Y_cov + np.diag(np.ones(len(idxG)-1))*0.001
                Y_cov *= 5.0/(len(idxG)-1)/(1+prob_isos[g].shape[0]/5000.0)

                result.append(pool.apply_async(MH_propose, (Y_now[idxG], 
                    Y_cov, prob_isos[g], len_isos[g], gCounts[g],
                    total_count, F_pre[idxG], sigma_in, gap, ftype)))
            pool.close()
            pool.join()

            for g in range(len(result)):
                idxG = range(Idx_all[g], Idx_all[g+1])
                _Y, _Psi, _Cnt = result[g].get()
                for k in idxG:
                    Y_all[k, idxT] = _Y[k-idxG[0], :]
                    Psi_all[k, idxT] = _Psi[k-idxG[0], :]
                    Cnt_all[k, idxT] = _Cnt[k-idxG[0], :]

        Y_now = Y_all[:,idxT[-1]]

        # Update weight and prediction suboptimally
        if ["RPK", "RPKM", "FPKM", "rpk", "rpkm", "fpkm"].count(ftype) == 1:
            F_now = Cnt_all[:,idxT[-1]] / tranLen / total_count * 10**9
            F_now = np.log10(F_now +0.01)
        elif ftype == "Y" or ftype == "y":
            F_now = Y_all[:,idxT[-1]]
            # F_now = Y_all[:,idxT[-1]/4: idxT[-1]].mean(axis=1)
        else:
            F_now = Psi_all[:,idxT[-1]]
        W_sub = np.dot(W_pt1, F_now[idxF]) # suboptimal
        if weights_in is not None:
            W_sub = weights_in
        F_pre[idxF] = np.dot(W_sub, X.T)
        W_all[:,m] = W_sub

        if _sigma is None or _sigma != _sigma:
            sigma_in = np.std(F_now[idxF] - F_pre[idxF])

            # sigma_in = np.std(Y_all[:,idxT[-1]/4: idxT[-1]]
            #.mean(axis=1)[F_idx] - F_pre[F_idx])
        else:
            sigma_in = _sigma

        #step 2. convergence diagnostics
        for k in range(Psi_all.shape[0]):
            Z = Geweke_Z(Psi_all[k, :(m+1)*gap])
            if Z is not None and Z <= 2:
                CONVERG[k] = True

        # show progress
        bar_len = 20
        run_time = time.time() - START_TIME
        percents = 100.0 * np.mean(CONVERG)
        filled_len = int(bar_len * percents / 100)
        bar = '=' * filled_len + '-' * (bar_len - filled_len)
        # sys.stdout.write('\r[%s] %.2f%% converged in %d run %.1f sec.' 
        #     % (bar, np.mean(CONVERG)*100, (m+1)*gap, run_time))
        sys.stdout.write('\r[Brie] [%s] %.1f%% converged in %d run %.1f sec. %.2f' 
            % (bar, np.mean(CONVERG)*100, (m+1)*gap, run_time, sigma_in))
        sys.stdout.flush()

        if sum(CONVERG) == len(CONVERG) and m*gap >= Mmin:
            W_all = W_all[:, :m]
            Y_all = Y_all[:, :(m+1)*gap]
            Psi_all = Psi_all[:, :(m+1)*gap]
            Cnt_all = Cnt_all[:, :(m+1)*gap]
            break
    print("")

    FPKM_all = Cnt_all / tranLen.reshape(-1, 1) / total_count * 10**9
    return Psi_all, Y_all, FPKM_all, Cnt_all, W_all, sigma_in
