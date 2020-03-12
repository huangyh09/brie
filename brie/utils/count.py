import sys
import numpy as np
from .sam_utils import load_samfile, fetch_reads

def _check_SE_event(gene):
    """Check SE event"""
    if (len(gene.trans) != 2 or 
        gene.trans[0].exons.shape[0] != 3 or
        gene.trans[1].exons.shape[0] != 2 or
        np.mean(gene.trans[0].exons[[0, 2], :] == 
                gene.trans[1].exons) != 1):
        return False
    else:
        return True


def _get_segment(exons, read):
    """Get the length of segments by devidinig a read into exons.
    The segments include one for each exon and two edges.
    """
    if read is None:
        return None
        
    _seglens = [0] * (exons.shape[0] + 2)
    _seglens[0] = np.sum(read.positions < exons[0, 0])
    _seglens[-1] = np.sum(read.positions > exons[-1, -1])
    for i in range(exons.shape[0]):
        _seglens[i + 1] = np.sum(
            read.positions >= exons[i, 0] and read.positions <= exons[i, 1])
    return _seglens


def check_reads_compatible(transcript, reads, edge_hang=10, junc_hang=2):
    """Check if reads are compatible with a transcript
    """
    is_compatible = [True] * len(reads)
    for i in range(len(reads)):
        _segs = _get_segment(transcript.exons, reads[i])
        
        # check mismatch to regions not in this transcript
        if len(reads[i].positions) - sum(_segs) >= junc_hang:
            is_compatible[i] = False
            continue

        # check if edge hang is too short
        if (_segs[0] > 0 or _segs[-1] > 0) and sum(_segs[1:-1]) < edge_hang:
            is_compatible[i] = False

    return np.array(is_compatible)


def SE_reads_count(gene, samfile, edge_hang=10, junc_hang=2, **kwargs):
    """Count the categorical reads mapped to a splicing event

    rm_duplicate=True, inner_only=True,
    mapq_min=0, mismatch_max=5, rlen_min=1, is_mated=True
    """
    # Check SE event
    if _check_SE_event(gene) == False:
        print("This is not exon-skipping event!")
        exit()

    # Fetch reads (TODO: customise fetch_reads function, e.g., FLAG)
    reads = fetch_reads(samfile, gene.chrom, gene.start, gene.stop, **kwargs)

    # Check reads compatible
    is_isoform1 = check_reads_compatible(gene.trans[0], reads["reads1"])
    is_isoform2 = check_reads_compatible(gene.trans[1], reads["reads1"])
    if len(reads["reads2"]) > 0:
        is_isoform1 *= check_reads_compatible(gene.trans[0], reads["reads2"])
        is_isoform2 *= check_reads_compatible(gene.trans[1], reads["reads2"])

    is_isoform1.append(check_reads_compatible(gene.trans[0], reads["reads1u"]))
    is_isoform2.append(check_reads_compatible(gene.trans[1], reads["reads1u"]))

    is_isoform1.append(check_reads_compatible(gene.trans[0], reads["reads2u"]))
    is_isoform2.append(check_reads_compatible(gene.trans[1], reads["reads2u"]))

    # return Reads matrix
    Rmat = np.zeros((len(is_isoform1), 2), dtype=bool)
    Rmat[:, 0] = is_isoform1
    Rmat[:, 1] = is_isoform2
    return Rmat


def get_count_matrix(g, g_idx, sam_files, edge_hang=10, junc_hang=2):
    RV = []
    for s in range(len(sam_files)):
        _RV = SE_reads_count(g, sam_files[s], edge_hang=10, junc_hang=2, 
            rm_duplicate=True, inner_only=False, mapq_min=0, mismatch_max=5, 
            rlen_min=1, is_mated=True)

        if _RV["Rmat"].shape[0] == 0:
            continue

        M = _RV["Rmat"]
        K = 2**(np.arange(M.shape[1]))
        code_id, code_cnt = np.unique(np.dot(M, K), return_counts=True)
        
        count_dict = {}
        for i in range(len(code_id)):
            count_dict["%d" %(code_id[i])] = code_cnt[i]
            
        RV.append("%d\t%d\t%s" %(g_idx + 1, s + 1, str(count_dict)))
    
    RV_line = ""
    if len(RV) > 0:
        RV_line = "\n".join(RV) + "\n"
        
    return RV_line



def SE_probability(gene, rlen=75, edge_hang=10, junc_hang=2):
    """Get read categorical probability of each isoform.
    In exon-skipping (SE) event, there are two isoform:
    isoform1 for exon inclusion and isoform2 for exon exclusion.
    
    Here, we only treat single-end reads. For paired-end reads,
    we treat it as the single-end by only using the most informative
    mate, namely the mate mapped to least number of isoform(s).
    
    isoform1: l1 + l2 + l3 + rlen - 2 * edge_hang
        p1: l2 + rlen - 2 * junc_hang
        p3: l1 + l3 - 2 * edge_hang + 2 * junc_hang
    isoform2: l1 + l3 + rlen - 2 * edge_hang
        p1: rlen - 2 * junc_hang
        p3: l1 + l3 - 2 * edge_hang + 2 * junc_hang
    """
    # check SE event
    if _check_SE_event(gene) == False:
        print("This is not exon-skipping event!")
        exit()
    
    l1, l2, l3 = gene.trans[0].exons[:, 1] - gene.trans[0].exons[:, 0]
    prob_mat = np.zeros((2, 3))
    
    # Isoform 1
    len_isoform1 = l1 + l2 + l3 + rlen - 2 * edge_hang
    prob_mat[0, 0] = (l2 + rlen - 2 * junc_hang) / len_isoform1
    prob_mat[0, 2] = (l1 + l3 - 2 * edge_hang + 2 * junc_hang) / len_isoform1
    
    # Isoform 2
    len_isoform2 = l1 + l3 + rlen - 2 * edge_hang
    prob_mat[1, 0] = (rlen - 2 * junc_hang) / len_isoform2
    prob_mat[1, 1] = (l1 + l3 - 2 * edge_hang + 2 * junc_hang) / len_isoform2
    
    return prob_mat

