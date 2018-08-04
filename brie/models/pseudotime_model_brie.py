
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
    #print("cov: ", cov)
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
    """A Matroplis-Hasting sampler with multivariate Gaussian proposal.


    """
    cnt = 0.0 # counter
    Y_try = np.zeros(Y_now.shape[0])
    Y_all = np.zeros((Y_now.shape[0], M)) # will contain all instances of Y for this round
    Psi_all = np.zeros((Y_now.shape[0], M)) # idem for psi
    Cnt_all = np.zeros((Y_now.shape[0], M)) # idem for counts ?

    Psi_now = np.exp(Y_now) / np.sum(np.exp(Y_now))
    Fsi_now = len_isos*Psi_now / np.sum(len_isos*Psi_now)
    Cnt_now = gene_Cnt*Fsi_now

    if ["RPK", "RPKM", "FPKM", "rpk", "rpkm", "fpkm"].count(ftype) == 1:
        # F_now = np.log2(Cnt_now* 10**9 / total_count +0.00001)
        # convert reads into a log of normalized reads to look like a Gaussian:
        F_now = Cnt_now / len_isos / total_count * 10**9
        F_now = np.log10(F_now +0.01)
    elif ftype == "Y" or ftype == "y":
        F_now = Y_now
    else:
        F_now = Psi_now

    P_now = np.log(np.dot(prob_isos, Fsi_now)).sum() # P(y*|R)
    #print("F_sigma : ", F_sigma)
    #print("F_sigma**2: ", F_sigma**2)
    for k in range(F_now.shape[0]):
        if F_pre[k] is None or  F_pre[k] != F_pre[k]:
            continue
        P_now += normal_pdf(F_now[k], F_pre[k], F_sigma**2)  # evaluate likelihood P(y|R) (denominator, l 12)

    for m in range(M): # for each baby step
        # propose
        Y_try[:-1] = np.random.multivariate_normal(Y_now[:-1], Y_cov) # Q(y*|y,eta), eta == Y_cov (line 10)
        Y_try[Y_try < -700] = -700 # min values treshold
        Y_try[Y_try > 700 ] = 700 # max values treshold
        Q_now = normal_pdf(Y_now[:-1], Y_try[:-1], Y_cov) # evaluate proposal (numerator of Q, l 12)
        Q_try = normal_pdf(Y_try[:-1], Y_now[:-1], Y_cov) # evaluate proposal (denominator of Q, l 12)

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
                P_try += normal_pdf(F_try[k], F_pre[k], F_sigma**2)  # evaluate likelihood P(y*|R) (numerator, l 12)

        # step 2: accept or reject the proposal
        alpha = np.exp(min(P_try+Q_now-P_now-Q_try, 0)) # 0=log(1) ( log of l 12 multiplication becomes sumation )
        if alpha is None:
            print("alpha is none!")
        elif np.random.rand(1) < alpha: # if mu < alpha, then accept (l 13)
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


def brie_MH_Heuristic(cell, feature_all, idxF, weights_in=None, _sigma=None,
                      _lambda=2.4, ftype="Y", M=10000, Mmin=1000, gap=10,
                      nproc=1): # _lambda=0.134
    """
    Matroplis-Hasting sampler to infer Brie parameters in a heuristic way,
    namely, sampling isoform proportion and optimize the bayesian regression
    in turn.

    Parameters
    ----------
    cell : dict
        keys are cell ids, values are dict carying values of each cell, such as:
        "R_mat": a list of array of bool
            each element in the list is a gene matrix (N*m) of reads identity
        "len_isos": a list of array of int
            each element in the list is an array for length of m isoforms
        "prob_isos": a list of array of float
            each element in the list is a gene matrix (N*m) for isoforms
            specific probability
        "total_count": int
            the total number of reads in the whole bam file(s)
    feature_all: an array 
        size (tranNum, K) for multiple isoforms or (tranNum/2, K) splicing 
        events
    idxF: int
        id of feature (Y here) to consider ? == id of relevent transcripts
    weights_in: array, (K+1)
        input weights for Bayesian regression, which may be learned from 
        better data sets. (+1 for the offset)
    _sigma: float
        the variance of residue between target and fitted prior
    _lambda: float
        the L2 constrain of Baysian regression
    ftype: string
        the type of target for regression: FPKM, Y or Psi. Psi=softmax(Y)
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

    ### initialisation of variables:
    cell_id_list = []
    for _id in cell: # for each cell
        # build list of cells' ids for matrix reference order
        cell_id_list.append(_id)

    ## define transcript and gene numbers from first cell data:
    _id = cell_id_list[0] # id of first cell
    
    R_mat = cell[_id]['R_mat']
    prob_isos = cell[_id]['prob_isos']
    len_isos = cell[_id]['len_isos']

    tranLen = np.array([])
    geneNum = len(len_isos) # number of considered genes
    for t in range(geneNum):
        R_mat[t], prob_isos[t], len_isos[t] = Iso_read_check(R_mat[t], 
            len_isos[t], prob_isos[t])
        tranLen = np.append(tranLen, len_isos[t])
    tranNum = len(tranLen)

    ## initialise global (across all cells) variables:
    WX_brie_all = np.zeros((len(cell_id_list), geneNum))#brie prediction for psi
    # collection of Ys for all cells (pseudotime component only):
    Y_t_NOW = np.zeros((len(cell_id_list), geneNum)) # each line is a cell
    T = np.zeros((len(cell_id_list), 1)) # column vector of pseudotimes
    for i in range(len(cell_id_list)): # for each cell
        _id = cell_id_list[i] # id of current cell
        WX_brie_all[i] = cell[_id]['WX']
        T[i,0] = cell[_id]['t']

    Y_NOW = WX_brie_all # complete collection of Y_now for all cells
    #if ftype == "Y" or ftype == "y":
    F_NOW = Y_NOW
    #else:
    #    raise InputError("brie_MH_Heuristic was given a value of ftype "
    #                     "different from Y or y: impossible in current "
    #                     "implementation.")
    
    #Psi_NOW = np.zeros((len(cell_id_list), geneNum)) # each line is a cell
    # actual feature vector (Y in original brie):
    # F_PRE = np.zeros((len(cell_id_list), geneNum)) # each line is a cell

    # define sigma:
    if _sigma == None:
        sigma_in = 1 # do not make it equal to zero
    else:
        sigma_in = _sigma # do not make it equal to zero
        
    # part one of W_t:
    W_t_part1 = np.dot(np.linalg.inv(np.dot(T.T, T) + _lambda*sigma_in**2), T.T)
    if weights_in is not None:
        W_t = weights_in
    else:
        W_t = W_t_part1.dot(Y_t_NOW) # weights are equal to zero at this stage

    W_t_all = [W_t] # list that store every W_t at the end of each iteration

    # compute new "Y bar":
    F_PRE = WX_brie_all + T.dot(W_t)

    for i in range(len(cell_id_list)): # for each cell
        _id = cell_id_list[i] # id of current cell

        # define relevant values:
        R_mat = cell[_id]['R_mat']
        prob_isos = cell[_id]['prob_isos']
        len_isos = cell[_id]['len_isos']
        total_count = cell[_id]['total_count']
        
        tranLen = np.array([])
        geneNum = len(len_isos) # number of considered genes
        for t in range(geneNum):
            R_mat[t], prob_isos[t], len_isos[t] = Iso_read_check(R_mat[t], 
                len_isos[t], prob_isos[t])
            prob_isos[t] = R_mat[t] * prob_isos[t]
            tranLen = np.append(tranLen, len_isos[t])
        tranNum = len(tranLen)

        Y_now = np.zeros(tranNum)
        Y_now[idxF] = Y_NOW[i,:] # 1 transcript over 2 # before: [::2]
        # (both monodimentional line vectors)
        Y_all = np.zeros((tranNum, M))
        Psi_now = np.zeros(tranNum)
        Psi_all = np.zeros((tranNum, M))
        Cnt_now = np.zeros(tranNum)
        Cnt_all = np.zeros((tranNum, M))
        gCounts = np.zeros(geneNum)
        Idx_all = np.zeros(geneNum+1, "int") # +1 because first id is for offset
        for g in range(len(len_isos)): # for each gene
            Idx_all[g+1] = Idx_all[g] + len(len_isos[g]) # id of a gene is id of
            # precedent gene plus the number of transcripts that precedent gene
            # has (so that the first transcript has the same id has its gene) 
            idxG = range(Idx_all[g], Idx_all[g+1]) # current transcript id range
            #print("Y_now: ", Y_now, "\nidxG: ", idxG)
            #print("\nY_now.shape: ", Y_now.shape)
            _psi = np.exp(Y_now[idxG]) / np.sum(np.exp(Y_now[idxG]))# psi vector
            _fsi = (len_isos[g]*_psi) / np.sum(len_isos[g]*_psi)

            gCounts[g] = prob_isos[g].shape[0]
            Psi_now[idxG] = _psi # psi for each transcript of idxG
            Cnt_now[idxG] = _fsi * gCounts[g]

        # initialize F_now
        if ["RPK", "RPKM", "FPKM", "rpk", "rpkm", "fpkm"].count(ftype) == 1:
            F_now = Cnt_now / tranLen / total_count * 10**9
            F_now = np.log10(F_now +0.01)
        elif ftype == "Y" or ftype == "y":
            F_now = np.zeros(tranNum)
            F_now[idxF] = F_NOW[i,:] # both monodimentional line vectors
            # (F_now[idxF] do not work)
        else:
            F_now = Psi_now

        # initialize F_pre
        F_pre = np.zeros(tranNum)
        F_pre[:] = None
        F_pre[idxF] = F_PRE[i,:] # for considered transcripts (1/2)
        
        # store relevant values:
        cell[_id]['R_mat'] = R_mat
        cell[_id]['prob_isos'] = prob_isos
        cell[_id]['len_isos'] = len_isos
        cell[_id]['F_pre'] = F_pre
        cell[_id]['Y_now'] = Y_now
        cell[_id]['Y_all'] = Y_all
        cell[_id]['Psi_now'] = Psi_now
        cell[_id]['Psi_all'] = Psi_all
        cell[_id]['Cnt_now'] = Cnt_now
        cell[_id]['Cnt_all'] = Cnt_all
        cell[_id]['gCounts'] = gCounts
        cell[_id]['Idx_all'] = Idx_all
        cell[_id]['F_now'] = F_now
            
        
    ### initialisation ends
    print("End of initialisation\n")

    # phase of proposals
    CONVERG = np.zeros(tranNum, "bool") # table of convergence
    for m in range(int(M/gap)): # for each giant step (i from 1 to n)
        print(f"\rbig step number {m}")
        idxT = range(m*gap, (m+1)*gap) # idxt = range of next baby steps
        # index for transcripts (T)
        # step 1: propose a value (original)
        for i in range(len(cell_id_list)): # for each cell
            _id = cell_id_list[i] # id of current cell
            
            # get useful variables
            prob_isos = cell[_id]['prob_isos']
            len_isos = cell[_id]['len_isos']
            F_pre = cell[_id]['F_pre']
            Y_now = cell[_id]['Y_now']
            Y_all = cell[_id]['Y_all']
            Psi_all = cell[_id]['Psi_all']
            Cnt_all = cell[_id]['Cnt_all']
            gCounts = cell[_id]['gCounts']
            Idx_all = cell[_id]['Idx_all']
            total_count = cell[_id]['total_count']
            
            if nproc == 1: # one process case
                for g in range(len(len_isos)): # for each gene (k from 1 to K)
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
                    for k in idxG: # corrections
                        Y_all[k, idxT] = _Y[k-idxG[0], :]
                        Psi_all[k, idxT] = _Psi[k-idxG[0], :]
                        Cnt_all[k, idxT] = _Cnt[k-idxG[0], :]

            # step 1: propose a value (parallelise)
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
            # update Y_NOW:
            Y_NOW[i,:] = Y_now[idxF] # monodimentional line vectors ([::2])

            # update F_now
            if ["RPK", "RPKM", "FPKM", "rpk", "rpkm", "fpkm"].count(ftype) == 1:
                F_now = Cnt_all[:,idxT[-1]] / tranLen / total_count * 10**9
                F_now = np.log10(F_now +0.01)
            elif ftype == "Y" or ftype == "y":
                F_now = Y_all[:,idxT[-1]] # == Y_now
            else:
                F_now = Psi_all[:,idxT[-1]]

            # upgrade F_NOW:
            F_NOW[i,:] = F_now[idxF] # monodimentional line vectors

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

            # store relevant values
            cell[_id]['Y_now'] = Y_now
            cell[_id]['Y_all'] = Y_all
            cell[_id]['Psi_all'] = Psi_all
            cell[_id]['Cnt_all'] = Cnt_all
            cell[_id]['gCounts'] = gCounts
            cell[_id]['Idx_all'] = Idx_all
            cell[_id]['F_now'] = F_now

            # wait for convergence and a minimum of steps
            if sum(CONVERG) == len(CONVERG) and m*gap >= Mmin:
                Y_all = Y_all[:, :(m+1)*gap]
                Psi_all = Psi_all[:, :(m+1)*gap]
                Cnt_all = Cnt_all[:, :(m+1)*gap]
                # store values
                cell[_id]['Y_all'] = Y_all
                cell[_id]['Cnt_all'] = Cnt_all
                cell[_id]['Psi_all'] = Psi_all
                break

        # update weights and Y bar:
        Y_t_NOW = Y_NOW - WX_brie_all # update pseudotemporal component of Y
        W_t = W_t_part1.dot(Y_t_NOW) # weights are equal to zero at this stage
        W_t_all.append(W_t) # store W_t
        F_PRE = WX_brie_all + T.dot(W_t)
        
        # update F_pre for each cell
        for i in range(len(cell_id_list)): # for each cell
            _id = cell_id_list[i] # id of current cell
            cell[_id]['F_pre'][idxF] = F_PRE[i,:]

            

        #     # store relevant values:
        #     cell[_id]['W_mat'] = W_mat
        #     cell[_id]['W_pt1'] = W_pt1
        #     cell[_id]['F_pre'] = F_pre
        #     cell[_id]['Y_now'] = Y_now
        #     cell[_id]['Y_all'] = Y_all
        #     cell[_id]['W_all'] = W_all
        #     cell[_id]['Psi_now'] = Psi_now
        #     cell[_id]['Psi_all'] = Psi_all
        #     cell[_id]['Cnt_now'] = Cnt_now
        #     cell[_id]['Cnt_all'] = Cnt_all
        #     cell[_id]['gCounts'] = gCounts
        #     cell[_id]['Idx_all'] = Idx_all
        #     cell[_id]['F_now'] = F_now
        #     cell[_id]['W_sub'] = W_sub
                   
    print("")

    # compute FPKM_all:
    for _id in cell:
        Cnt_all = cell[_id]['Cnt_all']
        total_count = cell[_id]['total_count']

        cell[_id]['FPKM_all'] = Cnt_all / tranLen.reshape(-1, 1) / total_count * 10**9
        
    return cell, W_t_all, sigma_in # Psi_all, Y_all, FPKM_all, Cnt_all, W_all, sigma_in

        # cell[id]['fractions.tsv']
