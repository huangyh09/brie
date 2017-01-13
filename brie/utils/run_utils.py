# bases functions supporting Brie

import os
import subprocess
import numpy as np

from .sam_utils import load_samfile
from .bias_utils import BiasFile, FastaFile
from .tran_utils import TranUnits, TranSplice


def set_info(g, sam_file, bias_mode, ref_file, bias_file, FLmean, FLstd,
    mate_mode, auto_min):
    RV = {}
    g = TranSplice(g)
    for ss in sam_file.split(","):
        _sam = load_samfile(ss)
        g.set_reads(_sam)
    if bias_mode != "unif":
        biasFile  = BiasFile(bias_file)
        fastaFile = FastaFile(ref_file)
        g.set_sequence(fastaFile)
        g.set_bias(biasFile, "seq") #under development
        if FLmean is None and biasFile.flen_mean != 0: 
            FLmean = biasFile.flen_mean
        if FLstd is None and biasFile.flen_std != 0:
            FLstd = biasFile.flen_std

    g.get_ready(bias_mode, FLmean, FLstd, mate_mode, auto_min)
    Rmat = g.Rmat
    if bias_mode == "unif": 
        len_iso  = g.efflen_unif
        prob_iso = g.proU
    else: 
        #len_iso = g.efflen_bias
        len_iso  = g.efflen_unif #under development
        prob_iso = g.proB
    RV["Rmat"] = Rmat
    RV["len_iso"] = len_iso
    RV["prob_iso"] = prob_iso
    return RV


def map_data(feature_file, tran_ids, log_out=False, add_intercept=True):
    """
    Format of feature file: genen_id tran_id feature_1 ...
    The feature file could contain only part of the transcriptome.
    """
    if ["hdf5", "h5", "HDF5", "H5"].count(feature_file.split(".")[-1]) == 1:
        ##Note: hdf5 is only supported with Python2
        import h5py
        f = h5py.File(feature_file, "r")
        feature = np.array(f["features"])
        feature_ids = np.array(f["factors"])

        # ids = np.array([x+".in" for x in f["gene_ids"]])
        ids = np.array([str(x.decode("utf-8"))+".in" for x in f["gene_ids"]])
        f.close()
    else:
        data = np.genfromtxt(feature_file, delimiter=",", dtype="str")
        ids = np.array([x+".in" for x in data[1:,0]])
        feature = data[1:, 1:].astype("float")
        feature_ids = data[0, 1:]

    feature_all = np.ones((len(tran_ids), feature.shape[1]))
    feature_all[:,:] = None

    idxF = []
    i, j = 0, 0
    idx1 = np.argsort(ids)
    idx2 = np.argsort(tran_ids)
    while j < len(idx2):
        if i >= len(idx1) or ids[idx1[i]] > tran_ids[idx2[j]]:
            feature_all[idx2[j], :] = None
            j += 1
        elif ids[idx1[i]] == tran_ids[idx2[j]]:
            idxF.append(j)
            feature_all[idx2[j], :] = feature[idx1[i], :]
            i += 1
            j += 1
        elif ids[idx1[i]] < tran_ids[idx2[j]]:
            i += 1

    if log_out is True:
        feature_all = np.log(feature_all)

    if add_intercept is True:
        feature_ids = np.append(feature_ids, "intercept")
        feature_all = np.append(feature_all, 
                                np.ones((feature_all.shape[0], 1)), axis=1)

    return feature_all, feature_ids, np.array(idxF, "int")


def get_CI(data, percent=0.95):
    """calculate the confidence intervals
    """
    if len(data.shape) <= 1:
        data = data.reshape(-1,1)
    RV = np.zeros((data.shape[1],2))
    CI_idx = int(data.shape[0] * (1-percent)/2)
    for k in range(data.shape[1]):
        temp = np.sort(data[:,k])
        RV[k,:] = [temp[-CI_idx], temp[CI_idx]]
    return RV


def save_data(out_dir, sample_num, gene_ids, tran_ids, tran_len, 
    feature_all, feature_ids, Psi_all, RPK_all, Cnt_all, W_all, sigma_):

    m1 = int(Psi_all.shape[1]*3/4)
    m2 = int(W_all.shape[1]*3/4)

    # save weights
    fid = open(os.path.join(out_dir, "weights.tsv"), "w")
    fid.writelines("feature_ids\tfeature_weights\n")
    for i in range(len(feature_ids)):
        fid.writelines("%s\t%.3e\n" %(feature_ids[i], W_all[i,-m2:].mean()))
    # fid.writelines("intercept\t%.3e\n" %W_all[-1,-m2:].mean())
    fid.writelines("#sigma\t%.3e\n" %sigma_)
    fid.close()

    # save psi
    fid = open(os.path.join(out_dir, "fractions.tsv"), "w")
    _line = "tran_id\tgene_id\ttransLen\tcounts\tFPKM\tPsi\tPsi_low\tPsi_high"
    fid.writelines(_line + "\n")
    for i in range(len(tran_ids)):
        psi_95 = get_CI(Psi_all[i,-m1:])[0,:]
        _line = "%s\t%s\t%d\t%.3e\t%.3e\t%.3f\t%.3f\t%.3f" %(tran_ids[i], 
            gene_ids[i], tran_len[i], Cnt_all[i,-m1:].mean(), 
            RPK_all[i,-m1:].mean(), Psi_all[i,-m1:].mean(), 
            psi_95[1], psi_95[0])
        fid.writelines(_line + "\n")
    fid.close()

    # save samples for all Psi
    if sample_num > 0:
        # import h5py
        # f = h5py.File(os.path.join(out_dir, "samples.h5"), "w")
        # f.create_dataset("gene_ids", data=gene_ids, compression="gzip")
        # f.create_dataset("tran_ids", data=tran_ids, compression="gzip")
        # f.create_dataset("features", data=feature_all, compression="gzip")
        # f.create_dataset("feature_ids", data=feature_ids, compression="gzip")
        # f.create_dataset("W_sample", data=W_all[:,-min(m2,sample_num):],
        #     compression="gzip", compression_opts=9)
        # f.create_dataset("Psi_sample", data=Psi_all[:,-min(m1,sample_num):],
        #     compression="gzip", compression_opts=9)
        # f.create_dataset("FPKM", data=RPK_all[:,-m1:].mean(axis=1),
        #     compression="gzip", compression_opts=9)
        # f.create_dataset("counts", data=Cnt_all[:,-m1:].mean(axis=1),
        #     compression="gzip", compression_opts=9)
        # f.create_dataset("sigma", data=np.array([sigma_]), compression="gzip")
        # f.close()

        W = W_all[:,-m2:].mean(axis=1)
        CNT = Cnt_all[:,-m1:].mean(axis=1)
        idx = np.arange(0, len(tran_ids), 2)
        priorY = np.zeros(len(tran_ids))
        priorY[idx] = np.dot(feature_all[idx,:], W)
        priorY[idx+1] = 0.0 - priorY[idx]
        
        samp_num = min(m1, sample_num)
        sample_file = os.path.join(out_dir, "samples.csv")
        fid = open(sample_file, "w")
        comment_line = "#tran_id,gene_id,count,prior_mean,prior_std,N_samples"
        fid.writelines(comment_line + "\n")
        for i in range(len(tran_ids)):
            name_part = "%s,%s" %(tran_ids[i], gene_ids[i])
            base_part = "%d,%.2e,%.2e" %(CNT[i], priorY[i], sigma_)
            data_part = ",".join(["%.2e" %(x) for x in Psi_all[i, -samp_num:]])
            fid.writelines(name_part + "," + base_part + "," + data_part + "\n")
        fid.close()

        bashCommand = "gzip -f %s" %(sample_file) 
        pro = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        output = pro.communicate()[0]
    