# This module is to process fasta files and some utils.
# download bigWigSummary binary file:
# http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigSummary

import pysam
import itertools
import subprocess
import numpy as np


class FastaFile:
    """docstring for FastaFile"""
    def __init__(self, fasta_file):
        self.f = pysam.FastaFile(fasta_file)

    def get_seq(self, qref, start, stop):
        """get the sequence in a given region, the start is from 1.
        The start and stop index may still need double check."""
        return self.f.fetch(qref, start-1, stop)


def rev_seq(seq):
    _tmp = []
    _tmp[:] = seq
    for j in range(len(_tmp)):
        if _tmp[j] == "A": _tmp[j] = "T"
        elif _tmp[j] == "T": _tmp[j] = "A"
        elif _tmp[j] == "G": _tmp[j] = "C"
        elif _tmp[j] == "C": _tmp[j] = "G"
    RV = "".join(_tmp[::-1])
    return RV


def get_motif(seq_full, motif, mode="counts"):
    """get the counts of motif in a sequence"""
    cnt = 0
    for i in range(len(seq_full)-len(motif)+1):
        if seq_full[i:i+len(motif)] == motif:
            cnt += 1
    if mode == "counts":
        return cnt
    elif mode == "frequency":
        return cnt / (len(seq_full)-len(motif)+1.0)
    elif mode == "normalized":
        return cnt / (len(seq_full)-len(motif)+1.0) / (0.25**len(motif))
    else:
        return None


def get_kmer_all(kmax=5, kmin=1, seqs="ATGC"):
    """generate kmers"""
    RV = []
    for i in range(kmin, kmax+1):
        for _seq in itertools.product(seqs, repeat=i): 
            RV.append("".join(_seq))
    return RV


def get_factorID(phast_in=False):
    """generate the factor names"""
    RV = ["SS5.I1", "SS3.I1", "SS5.I2", "SS3.I2", "logLen.C1", 
          "logLen.I1", "logLen.A", "logLen.I2", "logLen.C2", 
          "logLen.A_I1", "logLen.A_I2", "logLen.I1_I2"]

    regID = ["C1", "I1_5p", "I1_3p", "A", "I2_5p", "I2_3p", "C2"]
    if phast_in is True:
        RV += ["phastCons.%s" %x for x in regID]

    regK = [3, 2, 3, 4, 3, 2, 3]
    for i in range(len(regK)):
        kmers = get_kmer_all(kmax=regK[i], kmin=1, seqs="ATGC")
        RV += ["%s.%s" %(x, regID[i]) for x in kmers]

    return RV


def get_factor(tran, ref_file, phast_file):
    """Get sequence factors. Need bigWigSummary in $Path for Phast file.
    """
    # RV = np.zeros(735)
    RV = {}
    if tran.exonNum != 3:
        print("This is not a triplet of exons. Please check.")
        return RV
    else:
        exons = tran.exons
        chrom = tran.chrom
    fastaFile = FastaFile(ref_file)

    # logLength #8
    logLen = [np.log(exons[0,1]-exons[0,0]+1),  #C1
              np.log(exons[1,0]-exons[0,1]-1),  #I1
              np.log(exons[1,1]-exons[1,0]+1),  #A
              np.log(exons[2,0]-exons[1,1]-1),  #I2
              np.log(exons[2,1]-exons[2,0]+1)]  #C2
    if tran.strand != "+" and tran.strand != "1":
        logLen = logLen[::-1]
    logLen += [logLen[2]/logLen[1], logLen[2]/logLen[3], logLen[1]/logLen[3]]
        
    # Splice site score #4
    if tran.strand == "+" or tran.strand == "1":
        SS_seq = [fastaFile.get_seq(chrom, exons[0,1]-3,  exons[0,1]+8), #I1'5
                  fastaFile.get_seq(chrom, exons[1,0]-17, exons[1,0]+3), #I1'3
                  fastaFile.get_seq(chrom, exons[1,1]-3,  exons[1,1]+8), #I2'5
                  fastaFile.get_seq(chrom, exons[2,0]-17, exons[2,0]+3)] #I2'3
    else:
        SS_seq = [fastaFile.get_seq(chrom, exons[2,0]-8, exons[2,0]+3 ), #I1'5
                  fastaFile.get_seq(chrom, exons[1,1]-3, exons[1,1]+17), #I1'3
                  fastaFile.get_seq(chrom, exons[1,0]-8, exons[1,0]+3 ), #I2'5
                  fastaFile.get_seq(chrom, exons[0,1]-3, exons[0,1]+17)] #I2'3
        SS_seq = [rev_seq(x) for x in SS_seq]

    # PhastCons score #7
    regID = ["C1", "I1_5p", "I1_3p", "A", "I2_5p", "I2_3p", "C2"]
    regions = [[exons[0,0], exons[0,1]],        #C1
               [exons[0,1]+1, exons[0,1]+300],  #I1_5p
               [exons[1,0]-300, exons[1,0]-1],  #I1_3p
               [exons[1,0], exons[1,1]],        #A
               [exons[1,1]+1, exons[1,1]+300],  #I2_5p
               [exons[2,0]-300, exons[2,0]-1],  #I2_3p
               [exons[2,0], exons[2,1]]]        #C2
    if tran.strand != "+" and tran.strand != "1":
        regions = regions[::-1]

    cons_val = []
    if phast_file is not None:
        for i in range(len(regions)):
            bashCommand = "bigWigSummary %s %s %d %d 1" %(phast_file, chrom, 
                regions[i][0], regions[i][1]) 
            pro = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
            output = pro.communicate()[0]
            try:
                cons_val.append(float(output))
            except ValueError:
                cons_val.append(0.0) #np.nan
                print "No PhastCons data for %s. Treated as Zero." %tran.tranID
                # print output

    # 1 to 4mers #676
    regK = [3, 2, 3, 4, 3, 2, 3]
    kmer_frq = []
    for i in range(len(regions)):
        kmers = get_kmer_all(kmax=regK[i], kmin=1, seqs="ATGC")
        seq = fastaFile.get_seq(chrom, regions[i][0], regions[i][1])
        if tran.strand != "+" and tran.strand != "1":
            seq = rev_seq(seq)
        for j in range(len(kmers)):
            kmer_frq.append(get_motif(seq, kmers[j], mode="frequency"))

    RV["SS_seq"] = SS_seq
    RV["factor_val"] = np.array(logLen + cons_val + kmer_frq)
    return RV
    

def motif_score(msa, pwm_msa=None):
    """calculate motif scores

    Parameters
    ----------
    msa: list or array, element is string
    pwm_msa: list or array or None, element is string
        the msa for calculating pwm, with smooth
        if None, use msa, without smooth
        
    Return:
    -------
    score: array
        normalized motif score for each sequence
        100 means the best score according to pwm
        0 means the score for null pwm (by random)
        negetive score can happen when is poorer than null pwm
    """
    motif_len = len(msa[0])
    data = np.zeros((len(msa), motif_len), dtype="S1")
    for i in range(len(msa)):
        tmp = []
        tmp[:] = msa[i].upper()
        data[i,:] = tmp
    
    if pwm_msa is None: 
        pwmS = data
        pwm_add = 0.0
    else:
        pwm_add = 0.01 # for smooth the pwm
        pwmS = np.zeros((len(pwm_msa), motif_len), dtype="S1")
        for i in range(len(pwm_msa)):
            tmp = []
            tmp[:] = pwm_msa[i].upper()
            pwmS[i,:] = tmp
        
    pwm = np.zeros((4, motif_len))
    for i in range(motif_len):
        pwm[0,i] = (sum(pwmS[:,i]=="A")+pwm_add) / (pwmS.shape[0] + pwm_add*4)
        pwm[1,i] = (sum(pwmS[:,i]=="T")+pwm_add) / (pwmS.shape[0] + pwm_add*4)
        pwm[2,i] = (sum(pwmS[:,i]=="G")+pwm_add) / (pwmS.shape[0] + pwm_add*4)
        pwm[3,i] = (sum(pwmS[:,i]=="C")+pwm_add) / (pwmS.shape[0] + pwm_add*4)
        
    score = np.zeros(len(msa))
    s_max = np.sum(np.log2(pwm.max(axis=0)))
    #s_min = np.sum(np.log2(pwm.min(axis=0)))
    s_min = pwm.shape[1] * np.log2(1.0/pwm.shape[0]) #random is prefered as zero
    for i in range(data.shape[0]):
        for j in range(motif_len):
            if   data[i,j] == "A": score[i] += np.log2(pwm[0, j])
            elif data[i,j] == "T": score[i] += np.log2(pwm[1, j])
            elif data[i,j] == "G": score[i] += np.log2(pwm[2, j])
            elif data[i,j] == "C": score[i] += np.log2(pwm[3, j])
    score = (score - s_min) / (s_max-s_min) * 100
    
    return score

