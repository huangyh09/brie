# This module is to process different transcript state, and set information
# for the transcript states, and finally get ready to run MCMC sampling.

import sys
import numpy as np
from .sam_utils import load_samfile, fetch_reads
from .bias_utils import FastaFile, BiasFile
from .gtf_utils import Gene, Transcript

def normal_pdf(x, mu, sigma):
    RV = 1 / (sigma*np.sqrt(2*np.pi)) * np.exp(-1.0/2*((x-mu)/sigma)**2)
    return RV

class TranUnits:
    """docstring for TranUnits"""
    def __init__(self, transcript):
        self.chrom  = transcript.chrom
        self.strand = transcript.strand
        self.units  = transcript.exons
        self.loci   = np.array([],"int")
        for i in range(self.units.shape[0]):
            _loci = np.arange(self.units[i,0], self.units[i,1]+1)
            self.loci = np.append(self.loci, _loci)
        self.ulen = len(self.loci)

    def set_sequence(self, fastaFile=None):
        """set the sequence of the transcript units, with 20 bases longer in 
        each end"""
        if fastaFile is not None:
            self.seq = fastaFile.get_seq(self.chrom, self.units[0,0] - 20, 
                                                     self.units[0,0] - 1)
            for i in range(self.units.shape[0]):
                self.seq += fastaFile.get_seq(self.chrom, self.units[i,0], 
                                                          self.units[i,1])
            self.seq += fastaFile.get_seq(self.chrom, self.units[i,1] + 1, 
                                                      self.units[i,1] + 20)
        else:
            print("This is a Null sequence file.")

    def set_bias(self, biasFile=None, mode="both"):
        """set the bias for all loci, with the bias modes of sequence / 
        position / both; make sure setting quence before using sequence
        or both modes"""
        if biasFile is None: return
        self.bias_method = mode
        self.bias5 = np.ones(self.loci.shape[0])
        self.bias3 = np.ones(self.loci.shape[0])

        if ["seq", "sequence", "both"].count(mode) > 0:
            for i in range(len(self.bias5)):
                ipos = i + 20
                if self.strand == "+" or self.strand == "1":
                    _seq5 = self.seq[ipos-8  : ipos+13]
                    _seq3 = self.seq[ipos-12 : ipos+9 ][::-1]
                else:
                    _seq5 = self.seq[ipos-12 : ipos+9 ][::-1]
                    _seq3 = self.seq[ipos-8  : ipos+13]

                self.bias5[i] *= biasFile.get_seq_bias(_seq5, 5)
                self.bias3[i] *= biasFile.get_seq_bias(_seq3, 3)

        if ["pos", "position", "both"].count(mode) > 0:
            for i in range(len(self.bias5)):
                if self.strand == "+" or self.strand == "1":
                    _pos = i
                else:
                    _pos = len(self.loci) - 1 - i

                self.bias5[i] *= biasFile.get_pos_bias(_pos, self.ulen, 5)
                self.bias3[i] *= biasFile.get_pos_bias(_pos, self.ulen, 3)

    def get_index(self, loc):
        """get the location and the exon id of a loc on the transcript."""
        RV = [-1, -1]
        if   loc < self.units[0,0]:   RV[:] = [0, -2]
        elif loc > self.units[-1,-1]: RV[:] = [self.ulen-1, -3]
        else:
            cnt = 0
            for i in range(self.units.shape[0]):
                if loc >= self.units[i,0] and loc <= self.units[i,1]:
                    RV = [cnt + loc - self.units[i,0], i]
                    break
                else:
                    cnt += self.units[i,1] - self.units[i,0] + 1
        return RV

    def get_read_info(self, r1, r2, part_in=True):
        if r1 is None and r2 is None:
            return None
        idx5, idx3 = None, None
        mapq1, idx51, idx31 = 0.0, None, None
        mapq2, idx52, idx32 = 0.0, None, None
        if self.strand == "+" or self.strand == "1":
            if r1 is not None and r1.is_reverse:
                r1, r2 = r2, r1
            elif r2 is not None and r2.is_reverse == False:
                r1, r2 = r2, r1
        else:
            if r1 is not None and r1.is_reverse == False:
                r1, r2 = r2, r1
            elif r2 is not None and r2.is_reverse:
                r1, r2 = r2, r1

        if r1 is not None:
            mapq1 = 1.0 - 10 ** (0 - r1.mapq / 10.0)
            if self.strand == "+" or self.strand == "1":
                idx51 = self.get_index(r1.pos)
                idx31 = self.get_index(r1.aend - 1)
            else:
                idx31 = self.get_index(r1.pos)
                idx51 = self.get_index(r1.aend - 1)
            if idx51[1] <= -1 or idx31[1] <= -1: return None
            elif idx51[1] >= 0 and idx31[1] >= 0:
                if (abs(idx51[0] - idx31[0]) + 1 > r1.qlen + 3 or 
                    abs(idx51[0] - idx31[0]) + 1 < r1.qlen - 3): return None

        if r2 is not None:
            mapq2 = 1.0 - 10 ** (0 - r2.mapq / 10.0)
            if self.strand == "+" or self.strand == "1":
                idx52 = self.get_index(r2.pos)
                idx32 = self.get_index(r2.aend - 1)
            else:
                idx32 = self.get_index(r2.pos)
                idx52 = self.get_index(r2.aend - 1)
            if  idx52[1] <= -1 or idx32[1] <= -1: return None
            elif idx52[1] >= 0 and idx32[1] >= 0:
                if (abs(idx52[0] - idx32[0]) + 1 > r2.qlen + 3 or 
                    abs(idx52[0] - idx32[0]) + 1 < r2.qlen - 3): return None

        if r1 is None: 
            flen = abs(idx52[0] - idx32[0]) + 1
            if idx32[1] >= 0: idx3 = idx32[0]
        elif r2 is None:
            flen = abs(idx51[0] - idx31[0]) + 1
            if idx51[1] >= 0: idx5 = idx51[0]
        else:
            flen = abs(idx32[0] - idx51[0]) + 1
            if idx51[1] >= 0: idx5 = idx51[0]
            if idx32[1] >= 0: idx3 = idx32[0]
            
        RV = {}
        RV["idx5"] = idx5
        RV["idx3"] = idx3
        RV["flen"] = flen
        RV["prob"] = max(mapq1, mapq2)
        return RV

    def set_reads(self, reads1=[], reads2=[], bias_mode="unif",  
        flen_mean=None, flen_std=None):
        """identify whether a read (pair) or is in this units, and return the 
        identity of the reads in this units, the units specific fragment  
        length bias scores. The bias score can be based on both (mode=both) 
        ends or single end5 (mode=end5) or single end3 (mode=end5), uniform 
        (mode=unif). Make sure the loc of read1 is smaller than read2."""
        if len(reads1) == 0 and len(reads2) == 0:
            print("Please input paired-end reads or singled end reads!")
            sys.exit(1)
        elif (len(reads1) * len(reads2)) != 0 and len(reads2) != len(reads1):
            print("Please input the same number of both mates of the reads!")
            sys.exit(1)

        self.rcnt = max(len(reads1), len(reads2))
        self.idx5 = np.ones(self.rcnt, "float")
        self.idx3 = np.ones(self.rcnt, "float")
        self.Rmat = np.ones(self.rcnt, "bool")
        self.flen = np.ones(self.rcnt, "float")
        self.proB = np.ones(self.rcnt, "float")
        self.proU = np.ones(self.rcnt, "float")
        self.bias_mode = bias_mode
        self.efflen_bias = 0
        self.efflen_unif = 0

        for i in range(self.rcnt):
            if len(reads1) == 0: r1 = None
            else: r1 = reads1[i]
            if len(reads2) == 0: r2 = None
            else: r2 = reads2[i]

            rinfo = self.get_read_info(r1, r2)
            if rinfo is None: 
                self.Rmat[i] = False
                self.flen[i] = None
                self.proB[i] = None
                self.proU[i] = None
                self.idx5[i] = None
                self.idx3[i] = None
            else:
                self.Rmat[i] = True
                self.flen[i] = rinfo["flen"]
                self.proB[i] = rinfo["prob"]
                self.proU[i] = rinfo["prob"]
                self.idx5[i] = rinfo["idx5"]
                self.idx3[i] = rinfo["idx3"]
                if self.bias_mode == "unif": continue
                elif self.bias_mode != "end3" and rinfo["idx5"] is not None:
                    self.proB[i] *= self.bias5[rinfo["idx5"]]
                elif self.bias_mode != "end5" and rinfo["idx3"] is not None:
                    self.proB[i] *= self.bias3[rinfo["idx3"]]

        # fragement distribution
        flen = self.flen[self.Rmat]
        self.probs = np.zeros(self.ulen)
        if sum(self.Rmat) == 0: self.probs[0] = 1.0
        elif np.unique(flen).shape[0] >= 10:
            if flen_std is None: flen_std = np.std(flen)
            if flen_mean is None: flen_mean = np.mean(flen)
            x = np.arange(self.ulen) + 1
            self.probs[:] = normal_pdf(x, flen_mean, flen_std)
            if sum(self.probs) > 0: self.probs /= sum(self.probs)
        else:
            for i in np.unique(flen): 
                self.probs[int(i)-1] = np.mean(flen==i) #be careful here.
             
        # effective length
        self.biasLen = np.zeros(self.ulen)
        for i in range(1, self.ulen+1):
            # self.efflen_unif += self.probs[i-1] * (self.ulen-i+1)
            self.efflen_unif += 1
            if self.bias_mode == "unif" or self.probs[i-1] == 0: continue
            for j in range(self.ulen - i + 1):
                if self.strand == "+" or self.strand == "1":
                    pos5, pos3 = j, j+i-1
                else:
                    pos3, pos5 = j, j+i-1
                if   self.bias_mode == "end5": _bias = self.bias5[pos5]
                elif self.bias_mode == "end3": _bias = self.bias3[pos3]
                else : _bias = self.bias5[pos5] * self.bias3[pos3]
                self.biasLen[i-1] += _bias
            self.efflen_bias += self.probs[i-1] * self.biasLen[i-1]

        # reads probability
        for i in range(self.rcnt):
            if self.Rmat[i] == False: continue
            fL = int(self.flen[i])
            self.proU[i] *= self.probs[fL-1] / (self.ulen - fL + 1)
            if self.bias_mode != "unif":
                self.proB[i] *= (self.probs[fL-1] / self.biasLen[fL-1])

class TranSplice:
    def __init__(self, Gene):
        self.gene = Gene
        self.stop = Gene.stop
        self.chrom = Gene.chrom
        self.start = Gene.start
        self.unitSet = []
        for i in range(len(self.gene.trans)):
            self.unitSet.append(TranUnits(self.gene.trans[i]))
        self.read1p = []
        self.read2p = []
        self.read1u = []
        self.read2u = []

    def set_sequence(self, fastafile):
        """get the sequence from the genome sequence in FastaFile object"""
        for i in range(len(self.unitSet)):
            self.unitSet[i].set_sequence(fastafile)

    def set_bias(self, biasfile, mode="both"):
        """get the bias parameters from the BiasFile object
        mode: both, seq, pos"""
        for i in range(len(self.unitSet)):
            self.unitSet[i].set_bias(biasfile, mode)
        self.bias_in = True

    def set_reads(self, samfile, rm_duplicate=True, inner_only=True,
                  mapq_min=0, mismatch_max=5, rlen_min=1, is_mated=True):
        """fetch the reads on this transcript from sam file"""
        reads = fetch_reads(samfile, self.chrom, self.start, self.stop,  
                            rm_duplicate, inner_only, mapq_min, 
                            mismatch_max, rlen_min,   is_mated)

        self.read1p += reads["reads1"]
        self.read2p += reads["reads2"]
        self.read1u += reads["reads1u"]
        self.read2u += reads["reads2u"]

    def get_ready(self, bias_mode="unif", flen_mean=None, flen_std=None, 
                  mate_mode="pair", auto_min=200):
        """get the location index of the transcript, need implimentation
        in future. Then, we could remove the Rmat, and flen from ReadSet
        object, and set the Rmat and flen by the info of states."""
        if mate_mode == "single" or (mate_mode == "auto" and 
            len(self.read1p) < auto_min):
            self.read1u += self.read1p
            self.read2u += self.read2p
            self.read1p = []
            self.read2p = []
            
        rcnt = [len(self.read1p), len(self.read1u), len(self.read2u)]
        unit_cnt = len(self.unitSet)
        self.Rmat = np.ones((sum(rcnt), unit_cnt), "bool")
        self.flen = np.ones((sum(rcnt), unit_cnt), "float")
        self.proB = np.ones((sum(rcnt), unit_cnt), "float")
        self.proU = np.ones((sum(rcnt), unit_cnt), "float")
        self.efflen_bias = np.zeros(unit_cnt, "float")
        self.efflen_unif = np.zeros(unit_cnt, "float")

        for i in range(unit_cnt):
            _units = self.unitSet[i]
            for j in range(3):
                _idx = np.arange(sum(rcnt[:j]), sum(rcnt[:j+1]))
                _reads1, _reads2 = [], []
                if   j == 0: _reads1, _reads2 = self.read1p, self.read2p
                elif j == 1: _reads1 = self.read1u
                elif j == 2: _reads2 = self.read2u
                if (len(_reads1) + len(_reads2)) == 0: continue

                _units.set_reads(_reads1, _reads2, bias_mode, flen_mean, flen_std)
                self.Rmat[_idx, i] = _units.Rmat
                self.flen[_idx, i] = _units.flen
                self.proB[_idx, i] = _units.proB
                self.proU[_idx, i] = _units.proU

                self.efflen_unif[i] += sum(_units.Rmat) * _units.efflen_unif
                self.efflen_bias[i] += sum(_units.Rmat) * _units.efflen_bias

            if sum(self.Rmat[:, i]) > 0:
                self.efflen_unif[i] /= sum(self.Rmat[:, i])
                self.efflen_bias[i] /= sum(self.Rmat[:, i])
            else:
                self.efflen_unif[i] = _units.ulen
                self.efflen_bias[i] = _units.ulen

