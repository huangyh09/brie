"""Count reads in transcript for well-based scRNA-seq or bulk RNA-Seq"""

import sys
import time
import numpy as np
import multiprocessing
from .sam_utils import load_samfile, fetch_reads, check_pysam_chrom

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
        print("This is not exon-skipping event: %s!" %(gene.geneID))
        exit()
    
    l1, l2, l3 = gene.trans[0].exons[:, 1] - gene.trans[0].exons[:, 0]
    prob_mat = np.zeros((2, 3))
    
    # Isoform 1
    len_isoform1 = l1 + l2 + l3 + rlen - 2 * edge_hang
    prob_mat[0, 0] = (l2 + rlen - 2 * junc_hang) / len_isoform1
    prob_mat[0, 2] = (l1 + l3 - 2 * edge_hang + 2 * junc_hang) / len_isoform1
    
    # Isoform 2
    len_isoform2 = l1 + l3 + rlen - 2 * edge_hang
    prob_mat[1, 1] = (rlen - 2 * junc_hang) / len_isoform2
    prob_mat[1, 2] = (l1 + l3 - 2 * edge_hang + 2 * junc_hang) / len_isoform2
    
    return prob_mat


def SE_effLen(gene, rlen=75, edge_hang=10, junc_hang=2):
    """Get effective length matrix for three read categories from two isoforms.
    
    In exon-skipping (SE) event, there are two isoform:
    isoform1 for exon inclusion and isoform2 for exon exclusion.
    and three read groups:
    group1: uniquely from isoform1
    group2: uniquely from isoform2
    group3: ambiguous identity
    
    Here, we only treat single-end reads. For paired-end reads,
    we treat it as the single-end by only using the most informative
    mate, namely the mate mapped to least number of isoform(s).
    
    isoform1: l1 + l2 + l3 + rlen - 2 * edge_hang
        read group1: l2 + rlen - 2 * junc_hang
        read group3: l1 + l3 - 2 * edge_hang + 2 * junc_hang
    isoform2: l1 + l3 + rlen - 2 * edge_hang
        read group2: rlen - 2 * junc_hang
        read group3: l1 + l3 - 2 * edge_hang + 2 * junc_hang
    """
    # check SE event
    if _check_SE_event(gene) == False:
        print("This is not exon-skipping event: %s!" %(gene.geneID))
        exit()
    
    l1, l2, l3 = gene.trans[0].exons[:, 1] - gene.trans[0].exons[:, 0]
    isoLen_mat = np.zeros((2, 3))
    
    # isoform length
    len_isoform1 = l1 + l2 + l3 + rlen - 2 * edge_hang
    len_isoform2 = l1 + l3 + rlen - 2 * edge_hang
    
    # segments
    isoLen_mat[0, 0] = l2 + rlen - 2 * junc_hang
    isoLen_mat[1, 1] = rlen - 2 * junc_hang
    isoLen_mat[0, 2] = l1 + l3 - 2 * edge_hang + 2 * junc_hang
    isoLen_mat[1, 2] = l1 + l3 - 2 * edge_hang + 2 * junc_hang
    
    # prob_mat = isoLen_mat / isoLen_mat.sum(1, keepdims=True)
    
    return isoLen_mat


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
            (read.positions >= exons[i, 0]) * (read.positions <= exons[i, 1]))
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
            continue
            
        # check if exon has been skipped
        if len(_segs) > 4:
            for j in range(2, len(_segs) - 2):
                if (_segs[j-1] >= junc_hang and _segs[j+1] >= junc_hang and
                    transcript.exons[j-1, 1] - transcript.exons[j-1, 0] - 
                    _segs[j] >= junc_hang):
                    is_compatible[i] = False
                    break

    return np.array(is_compatible)


def fetch_reads_count(gene, sam_file, event_type="SE", RNA_type="spliced",
                      edge_hang=10, junc_hang=2, **kwargs):
    """Count the categorical reads mapped to a splicing event

    rm_duplicate=True, inner_only=True,
    mapq_min=0, mismatch_max=5, rlen_min=1, is_mated=True
    """
    # Check SE event
    if event_type == "SE" and _check_SE_event(gene) == False:
        print("This is not exon-skipping event!")
        exit()

    samFile, _chrom = check_pysam_chrom(sam_file, gene.chrom)

    # Fetch reads (TODO: customise fetch_reads function, e.g., FLAG)
    reads = fetch_reads(samFile, _chrom, gene.start, gene.stop, **kwargs)

    # Check reads compatible
    n_readsPE = len(reads["reads1"])
    n_readsU1 = len(reads["reads1u"])
    n_readsU2 = len(reads["reads2u"])
    
    n_reads = n_readsPE + n_readsU1 + n_readsU2
    n_trans = len(gene.trans)
    
    Rmat = np.zeros((n_reads , n_trans), dtype=bool)
    for i in range(n_trans):
        if RNA_type.upper() == "UNSPLICED":
            _tran = gene.trans[i].make_premRNA()
        else:
            _tran = gene.trans[i]
            
        idx_PE = np.arange(0, n_readsPE)
        idx_U1 = np.arange(n_readsPE, n_readsPE + n_readsU1)
        idx_U2 = np.arange(n_readsPE + n_readsU1, n_reads)
        
        Rmat[idx_PE, i] = check_reads_compatible(_tran, reads["reads1"])
        if len(reads["reads2"]) > 0:
            Rmat[idx_PE, i] *= check_reads_compatible(_tran, reads["reads2"])
            
        Rmat[idx_U1, i] = check_reads_compatible(_tran, reads["reads1u"])
        Rmat[idx_U2, i] = check_reads_compatible(_tran, reads["reads2u"])
    
    return Rmat

def get_count_matrix(**kwargs):
    print('Warning: get_count_matrix() is now renamed by _count_one_bam()')
    return _count_one_bam(**kwargs)


def _count_one_bam(genes, sam_file, sam_num, event_type="SE", 
                   edge_hang=10, junc_hang=2):
    """Counting reads for all genes in one bam file
    """
    RV = []
    for g in range(len(genes)):
        _Rmat = fetch_reads_count(
            genes[g], sam_file, event_type, edge_hang=edge_hang, 
            junc_hang=junc_hang, 
            rm_duplicate=True, inner_only=False, mapq_min=0, trimLen_max=5, 
            rlen_min=1, is_mated=True
        )

        if _Rmat.shape[0] == 0:
            continue

        K = 2**(np.arange(_Rmat.shape[1]))
        code_id, code_cnt = np.unique(np.dot(_Rmat, K), return_counts=True)
        
        count_dict = {}
        for i in range(len(code_id)):
            count_dict["%d" %(code_id[i])] = code_cnt[i]
            
        RV.append("%d\t%d\t%s\n" %(sam_num + 1, g + 1, str(count_dict)))
    return RV


### Parallel computing
def get_smartseq_matrix(genes, sam_table, out_dir, event_type="SE", 
    edge_hang=10, junc_hang=2, nproc=1, verbose=False):
    """
    Wrap function to run get_count_matrix in multiple processing
    """
    global TOTAL_BAMs, START_TIME, FID, PROCESSED
    FID = None
    PROCESSED = 0
    START_TIME = time.time()
    TOTAL_BAMs = sam_table.shape[0]

    def _show_progress(RV=None):    
        global PROCESSED, TOTAL_BAMs, START_TIME, FID
        if RV is not None: 
            FID.writelines(RV)

            PROCESSED += 1
            bar_len = 20
            run_time = time.time() - START_TIME
            percents = 100.0 * PROCESSED / TOTAL_BAMs
            filled_len = int(bar_len * percents / 100)
            bar = '=' * filled_len + '-' * (bar_len - filled_len)
            
            sys.stdout.write('\r[BRIE2] [%s] %.1f%% cells done in %.1f sec.' 
                % (bar, percents, run_time))
            sys.stdout.flush()
        return RV

    # write results into file
    FID = open(out_dir + "/read_count.mtx", "w")
    FID.writelines("%" + "%MatrixMarket matrix coordinate integer general\n")
    FID.writelines("%d\t%d\t%d\n" %(sam_table.shape[0], len(genes), 0))
    
    if nproc <= 1:
        for s in range(len(sam_table[:, 0])):
            sam_file = sam_table[s, 0]
            res = _count_one_bam(genes, sam_file, s, event_type, edge_hang, 
                                 junc_hang)
            _show_progress(res)
    else:
        pool = multiprocessing.Pool(processes=nproc)
        result = []
        for s in range(len(sam_table[:, 0])):
            sam_file = sam_table[s, 0]
            result.append(pool.apply_async(_count_one_bam, 
                (genes, sam_file, s, event_type, edge_hang, junc_hang), 
                callback=_show_progress))
        pool.close()
        pool.join()
    
    FID.close()

    print("")
    print("[BRIE2] %d genes have been processed." %(len(genes)))
    return None
