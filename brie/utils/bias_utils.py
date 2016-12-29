# This module is to estimate the sequence or position biases. It provides
# ways to access and save the bias parameters.

import pysam
import numpy as np
import pylab as pl

def norm_pdf(x, mu, sigma):
    return 1 / (sigma*np.sqrt(2*np.pi)) * np.exp(-1/2*((x-mu)/sigma)**2)

class FastaFile:
    """docstring for FastaFile"""
    def __init__(self, fasta_file):
        self.f = pysam.FastaFile(fasta_file)

    def get_seq(self, qref, start, stop):
        """get the sequence in a given region, the start is from 1.
        The start and stop index may still need double check."""
        return self.f.fetch(qref, start-1, stop)

class BiasFile:
    """docstring for BiasFile"""
    def __init__(self, bias_file=None):
        """get the bias parameters from the hdf5 file"""
        self.set_base_chain()
        self.pos5_bias = np.zeros((5, 20))
        self.pos3_bias = np.zeros((5, 20))
        self.pos5_unif = np.zeros((5, 20))
        self.pos3_unif = np.zeros((5, 20))
        self.pos5_prob = np.zeros((5, 20))
        self.pos3_prob = np.zeros((5, 20))
        self.percentile = np.zeros((5, 2))
        self.flen_mean, self.flen_std = 0, 0
        self.flen_sum1, self.flen_sum2 = 0, 0
        self.read_num = 0

        self.seq5_bias, self.seq3_bias = {}, {}
        self.seq5_unif, self.seq3_unif = {}, {}
        self.seq5_prob, self.seq3_prob = {}, {}
        for i in range(len(self.chain_len)):
            self.seq5_bias[str(i)] = np.zeros(4**self.chain_len[i])
            self.seq3_bias[str(i)] = np.zeros(4**self.chain_len[i])
            self.seq5_unif[str(i)] = np.zeros(4**self.chain_len[i])
            self.seq3_unif[str(i)] = np.zeros(4**self.chain_len[i])
            self.seq5_prob[str(i)] = np.zeros(4**self.chain_len[i])
            self.seq3_prob[str(i)] = np.zeros(4**self.chain_len[i])
        
        if bias_file is None: return
        
        fid = open(bias_file, "r")
        all_lines = fid.readlines()
        fid.close()
        
        self.flen_mean = float(all_lines[4].split("\t")[0])
        self.flen_std  = float(all_lines[4].split("\t")[1])
        self.flen_sum1 = float(all_lines[4].split("\t")[2])
        self.flen_sum2 = float(all_lines[4].split("\t")[3])
        self.read_num  = float(all_lines[4].split("\t")[4])
        for i in range(5,105):
            a, b = (i-5) // 20, (i-5) % 20
            if b == 0:
                self.percentile[a,:] = all_lines[i].split("|")[0].split("-")
            self.pos5_bias[a,b] = all_lines[i].split("\t")[1]
            self.pos3_bias[a,b] = all_lines[i].split("\t")[2]
            self.pos5_unif[a,b] = all_lines[i].split("\t")[3]
            self.pos3_unif[a,b] = all_lines[i].split("\t")[4]
            self.pos5_prob[a,b] = max(0, self.pos5_bias[a,b] / self.pos5_unif[a,b])
            self.pos3_prob[a,b] = max(0, self.pos3_bias[a,b] / self.pos3_unif[a,b])
            # self.pos5_prob[a,b] = self.pos5_bias[a,b] / self.pos5_unif[a,b]
            # self.pos3_prob[a,b] = self.pos3_bias[a,b] / self.pos3_unif[a,b]

        ii, cnt = all_lines[105].split("|")[0], -1
        for i in range(105,849):
            if ii == all_lines[i].split("|")[0]: 
                cnt += 1
            else:
                ii = all_lines[i].split("|")[0]
                cnt = 0
            self.seq5_bias[ii][cnt] = all_lines[i].split("\t")[1]
            self.seq3_bias[ii][cnt] = all_lines[i].split("\t")[2]
            self.seq5_unif[ii][cnt] = all_lines[i].split("\t")[3]
            self.seq3_unif[ii][cnt] = all_lines[i].split("\t")[4]
            self.seq5_prob[ii][cnt] = max(0, self.seq5_bias[ii][cnt] / self.seq5_unif[ii][cnt])
            self.seq3_prob[ii][cnt] = max(0, self.seq3_bias[ii][cnt] / self.seq3_unif[ii][cnt])
            self.base_chain[ii][cnt] = all_lines[i].split("\t")[0].split("|")[1]
            
    def add_bias_file(self, BF):
        self.pos5_bias += BF.pos5_bias
        self.pos3_bias += BF.pos3_bias
        self.pos5_unif += BF.pos5_unif
        self.pos3_unif += BF.pos3_unif
        for i in range(len(self.chain_len)):
            self.seq5_bias[str(i)] += BF.seq5_bias[str(i)]
            self.seq3_bias[str(i)] += BF.seq3_bias[str(i)]
            self.seq5_unif[str(i)] += BF.seq5_unif[str(i)]
            self.seq3_unif[str(i)] += BF.seq3_unif[str(i)]

        self.read_num  += BF.read_num
        self.flen_sum1 += BF.flen_sum1
        self.flen_sum2 += BF.flen_sum2
        if self.read_num > 0:
            self.flen_mean = self.flen_sum1 / (self.read_num+0.0)
            self.flen_std  = np.sqrt(self.flen_sum2*self.read_num - self.flen_sum1**2) / (self.read_num+0.0)

    def updata_prob(self):
        # self.pos5_prob = self.pos5_bias[a,b] / self.pos5_unif[a,b]
        # self.pos3_prob = self.pos3_bias[a,b] / self.pos3_unif[a,b]
        self.pos5_prob = self.pos5_bias / self.pos5_unif
        self.pos3_prob = self.pos3_bias / self.pos3_unif
        for i in range(len(self.chain_len)):
            self.seq5_prob[str(i)] = self.seq5_bias[str(i)] / self.seq5_unif[str(i)]
            self.seq3_prob[str(i)] = self.seq3_bias[str(i)] / self.seq3_unif[str(i)]
        self.flen_mean = self.flen_sum1 / (self.read_num+0.0)
        self.flen_std  = np.sqrt(self.flen_sum2*self.read_num - self.flen_sum1**2) / (self.read_num+0.0)

    def set_base_chain(self):
        """set the sub-base chain for the variable-length Markov model (VLMM),
        which was proposed by Reberts et al, Genome Biology, 2011: 
        Figure2 in supp 3. http://genomebiology.com/2011/12/3/r22/"""
        b1 = ["A","T","G","C"]
        b2, b3 = [], []
        for i in b1:
            for j in b1:
                b2.append(j+i)
                for k in b1:
                    b3.append(k+j+i)
        base_comb = [b1, b2, b3]

        self.chain_len = [1]*4 + [2]*3 + [3]*10 + [2]*2 + [1]*2
        self.base_chain = {}
        for i in range(21):
            self.base_chain[str(i)] = base_comb[self.chain_len[i]-1]

    def get_both_bias(self, seq, loc, ulen, end_num=5):
        """get the bias from the bias parameters"""
        prob = (self.get_seq_bias(seq, end_num) * 
                self.get_pos_bias(loc, ulen, end_num))
        return prob

    def get_seq_bias(self, seq, end_num):
        """get the sequence bias score"""
        if end_num == 5:
            parameters = self.seq5_prob
        elif end_num == 3:
            parameters = self.seq3_prob
        else:
            print("wrong end_num: %s" %str(end_num))
            return None

        prob = 1.0
        for j in range(len(seq)):
            _len = self.chain_len[j]
            _bas = seq[j-_len+1 : j+1]
            if self.base_chain[str(j)].count(_bas) == 0: continue
            _idx = self.base_chain[str(j)].index(_bas)
            prob = prob * parameters[str(j)][_idx]
        return prob

    def get_pos_bias(self, loc, ulen, end_num):
        """get the position bias score, the loc is base pair distance
        from the 5'end of the units"""
        if end_num == 5:
            parameters = self.pos5_prob
        elif end_num == 3:
            parameters = self.pos3_prob
        else:
            print("wrong end_num: %s" %str(end_num))
            return None

        bin1 = (ulen >= self.percentile[:,0]) * (ulen <= self.percentile[:,1])
        bin2 = 20.0 * loc / ulen
        prob = parameters[bin1, bin2]
        return prob

    def set_percentile(self, ulen, K=5):
        """set the percentiles by input the lengths of unitsets, i.e., ulen, 
        and number of percentiles, K."""
        perc_gap = np.linspace(0, 100, K+1)
        _percent = np.percentile(ulen, list(perc_gap))
        self.percentile = np.zeros((K, 2))
        for i in range(K):
            self.percentile[i, 0] = int(_percent[i])+1
            self.percentile[i, 1] = int(_percent[i+1])
            if i == 0:
                self.percentile[i,0] = 0
            elif i==4:
                self.percentile[i,1] = float("inf")

    def set_both_bias(self, seq, loc, ulen, weight, end_num=5, mode="bias"):
        """get the bias from the bias parameters"""
        self.set_seq_bias(seq, weight, end_num, mode)
        self.set_pos_bias(loc, ulen, weight, end_num, mode)

    def set_seq_bias(self, seq, weight, end_num=5, mode="bias"):
        """get the sequence bias score"""
        for j in range(len(seq)):
            _len = self.chain_len[j]
            _bas = seq[j-_len+1 : j+1]
            if self.base_chain[str(j)].count(_bas) == 0: continue
            _idx = self.base_chain[str(j)].index(_bas)
            if end_num == 5:
                if mode == "bias":
                    self.seq5_bias[str(j)][_idx] += weight
                elif mode == "unif":
                    self.seq5_unif[str(j)][_idx] += weight
            else:
                if mode == "bias":
                    self.seq3_bias[str(j)][_idx] += weight
                elif mode == "unif":
                    self.seq3_unif[str(j)][_idx] += weight

    def set_pos_bias(self, loc, ulen, weight, end_num=5, mode="bias"):
        """get the position bias score, the loc is base pair distance
        from the 5'end of the units"""
        bin1 = (ulen >= self.percentile[:,0]) * (ulen <= self.percentile[:,1])
        bin2 = int(20.0 * loc / (ulen + 0.0001))
        if end_num == 5:
            if mode == "bias":
                self.pos5_bias[bin1, bin2] += weight
            elif mode == "unif":
                self.pos5_unif[bin1, bin2] += weight
        else:
            if mode == "bias":
                self.pos3_bias[bin1, bin2] += weight
            elif mode == "unif":
                self.pos3_unif[bin1, bin2] += weight

    def save_file(self, out_file="out_file.bias"):
        """to save the bias file in BIAS FILE FORMAT"""
        fid = open(out_file, "w")
        fid.writelines("# BIAS PARAMETER FORMAT\n")
        fid.writelines("# fragment leng: 5 (mean, std, sum_fl, sum_fl^2, reads), line 5\n")
        fid.writelines("# position bias: 5*20*4 (name, b5, b3, u5, u3), line 6-105\n")
        fid.writelines("# sequence bias: 744*4 (name, b5, b3, u5, u3), line 106-849\n")
        fid.writelines("%.2f\t%.2f\t%.2e\t%.2e\t%.0f\n" %(self.flen_mean, self.flen_std,
            self.flen_sum1, self.flen_sum2, self.read_num))
        for i in range(self.pos5_bias.shape[0]):
            for j in range(self.pos5_bias.shape[1]):
                aLine = ("%.0f-%.0f|%d\t%.2e\t%.2e\t%.2e\t%.2e\n"
                         %(self.percentile[i,0], self.percentile[i,1], j, self.pos5_bias[i,j], 
                           self.pos3_bias[i,j], self.pos5_unif[i,j], self.pos3_unif[i,j]))
                fid.writelines(aLine)
        for i in sorted(self.base_chain.keys(), key=float):
            for j in range(len(self.base_chain[i])):
                aLine = ("%s|%s\t%.2e\t%.2e\t%.2e\t%.2e\n"
                         %(i, self.base_chain[i][j], self.seq5_bias[i][j], 
                           self.seq3_bias[i][j], self.seq5_unif[i][j], self.seq3_unif[i][j]))
                fid.writelines(aLine)
        fid.close()

    def plot_bias(self, mode=None):
        """plot of bias parameters: flen, pos5, pos3, seq5, seq3"""
        #fragment distribution
        if mode == "flen":
            xx = np.arange(0, 1000)
            yy = norm_pdf(xx, self.flen_mean, self.flen_std)
            pl.fill(xx, yy, 'k')#, linewidth=2.0)
            pl.xlabel("fragment length")
            pl.ylabel("$p(L)$")
            pl.xlim(0, 400)

        #position bias
        if mode == "pos5" or mode == "pos3":
            pl.plot(np.arange(20), np.ones(20), '--k')
            for i in range(5):
                _label="bin%d: %.0f-%.0f bp" %(i+1, self.percentile[i,0], self.percentile[i,1])
                if mode == "pos5":
                    pl.plot(np.arange(20)+0.5, self.pos5_prob[i,:], linewidth=2.0, label=_label)
                else:
                    pl.plot(np.arange(20)+0.5, self.pos3_prob[i,:], linewidth=2.0, label=_label)
            pl.legend(loc="best")
            pl.xlabel("fractional transcription position")
            pl.ylabel("bias weight")
            pl.ylim(0,2)

        #sequence bias
        if mode == "seq5" or mode == "seq3":
            base = ["A", "T", "G", "C"]
            _color = ["g", "r", "orange", "b"]
            if mode == "seq5":
                pl.plot(np.arange(21)-8, np.ones(21), '--k')
                pl.plot(np.zeros(2), np.array([0, 2.0]), '--k', linewidth=2.0)
                percent = np.zeros((4,21))
                for i in range(4):
                    for j in range(21):
                        _seq_bias = self.seq5_prob[str(j)]
                        percent[i,j] = np.sum(_seq_bias[i*4**(self.chain_len[j]-1): 
                            (i+1)*4**(self.chain_len[j]-1)]) / 4**(self.chain_len[j]-1)
                    pl.plot(np.arange(21)-8, percent[i,:], ":o", c=_color[i], label=base[i])
                pl.xlabel("offset from 3' fragment end")
                pl.xlim(-8,12)
            else:
                pl.plot(np.arange(21)-12, np.ones(21), '--k')
                pl.plot(np.zeros(2), np.array([0, 2.0]), '--k', linewidth=2.0)
                percent = np.zeros((4,21))
                for i in range(4):
                    for k in range(21):
                        j = 20 - k
                        _seq_bias = self.seq3_prob[str(j)]
                        percent[i,j] = np.sum(_seq_bias[i*4**(self.chain_len[j]-1): 
                            (i+1)*4**(self.chain_len[j]-1)]) / 4**(self.chain_len[j]-1)
                    pl.plot(np.arange(21)-12, percent[i,:], ":o", c=_color[i], label=base[i])
                pl.xlabel("offset from 3' fragment end")
                pl.xlim(-12,8)
            pl.legend(loc="best")
            pl.xlabel("offset from %s' fragment end" %mode[3])
            pl.ylabel("bias weight")
            pl.ylim(0.5,2)

        #legend only
        if mode == "legend":
            base = ["A", "T", "G", "C"]
            _color = ["g", "r", "orange", "b"]
            pl.axis('off')
            ax1 = pl.twinx()
            for i in range(len(base)):
                ax1.plot([], [], "o", c=_color[i], label=base[i])
            ax1.legend(numpoints=1, loc=4)
            ax1.axis('off')
            ax2 = pl.twinx()
            for i in range(5):
                _label="bin%d: %.0f-%.0f bp" %(i+1, self.percentile[i,0], self.percentile[i,1])
                ax2.plot([], [], linewidth=2.0, label=_label)
            ax2.legend(loc=3)
            ax2.axis('off')

    def plot_bias_full(self):
        pl.subplot(3,2,1)
        self.plot_bias(mode="flen")
        pl.subplot(3,2,2)
        self.plot_bias(mode="legend")
        pl.subplot(3,2,3)
        self.plot_bias(mode="pos5")
        pl.legend().set_visible(False)
        pl.subplot(3,2,4)
        self.plot_bias(mode="pos3")
        pl.legend().set_visible(False)
        pl.subplot(3,2,5)
        self.plot_bias(mode="seq5")
        pl.legend().set_visible(False)
        pl.subplot(3,2,6)
        self.plot_bias(mode="seq3")
        pl.legend().set_visible(False)
