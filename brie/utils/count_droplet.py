"""Count reads in transcript from droplet based data"""

import sys
import time
import numpy as np
import multiprocessing
from .base_utils import match
from .sam_utils import load_samfile, fetch_reads, check_pysam_chrom
from .count import check_reads_compatible, _check_SE_event


def get_droplet_UMIcount(gene, samFile, event_type="SE", edge_hang=10, 
        junc_hang=2, CB_tag='CB', UMI_tag='UR', verbose=False, **kwargs):
    """Count the categorical reads mapped to a splicing event

    rm_duplicate=True, inner_only=True,
    mapq_min=0, mismatch_max=5, rlen_min=1, is_mated=True
    """
    # Check SE event
    if event_type == "SE" and _check_SE_event(gene) == False:
        print("This is not exon-skipping event!")
        exit()

    # Fetch reads (TODO: customise fetch_reads function, e.g., FLAG)
    reads = fetch_reads(samFile, gene.chrom, gene.start, gene.stop, **kwargs)
    _n_reads = (len(reads["reads1u"]) + len(reads["reads2u"]) + 
                len(reads["reads1"]) + len(reads["reads2"]))

    if verbose == True or verbose >= 1:
        print("%d reads fetched on %s" %(_n_reads, gene.geneID))
        print("R1, R2, paired_R1, paired_R2:", len(reads["reads1u"]), 
            len(reads["reads2u"]), len(reads["reads1"]), len(reads["reads2"]))

    # Check reads have CB_tag and UMI_tag
    reads["reads1u"] = [r for r in reads["reads1u"] if r.has_tag(CB_tag)]
    reads["reads2u"] = [r for r in reads["reads2u"] if r.has_tag(CB_tag)]
    reads["reads1u"] = [r for r in reads["reads1u"] if r.has_tag(UMI_tag)]
    reads["reads2u"] = [r for r in reads["reads2u"] if r.has_tag(UMI_tag)]

    reads["reads1"] = [r for r in reads["reads1"] if r.has_tag(CB_tag)]
    reads["reads2"] = [r for r in reads["reads2"] if r.has_tag(CB_tag)]
    reads["reads1"] = [r for r in reads["reads1"] if r.has_tag(UMI_tag)]
    reads["reads2"] = [r for r in reads["reads2"] if r.has_tag(UMI_tag)]
    
    _n_reads = (len(reads["reads1u"]) + len(reads["reads2u"]) + 
                len(reads["reads1"]) + len(reads["reads2"]))

    if verbose == True or verbose >= 1:
        print("%d reads fetched on %s with CB & UMI" %(_n_reads, gene.geneID))
        print("R1, R2, paired_R1, paired_R2:", len(reads["reads1u"]), 
            len(reads["reads2u"]), len(reads["reads1"]), len(reads["reads2"]))

    # Check reads compatible
    n_readsPE = len(reads["reads1"])
    n_readsU1 = len(reads["reads1u"])
    n_readsU2 = len(reads["reads2u"])
    
    n_reads = n_readsPE + n_readsU1 + n_readsU2
    n_trans = len(gene.trans)
    
    R_CB = []
    R_UR = []
    if n_readsPE > 0:
        print('Warning: here assumes mate1 & ' +
              'mate2 have the same cell & UMI barocdes.')

    R_UR += [x.get_tag(UMI_tag) for x in reads["reads1"]]
    R_UR += [x.get_tag(UMI_tag) for x in reads["reads1u"]]
    R_UR += [x.get_tag(UMI_tag) for x in reads["reads2u"]]

    R_CB += [x.get_tag(CB_tag) for x in reads["reads1"]]
    R_CB += [x.get_tag(CB_tag) for x in reads["reads1u"]]
    R_CB += [x.get_tag(CB_tag) for x in reads["reads2u"]]

    Rmat = np.zeros((n_reads , n_trans), dtype=bool)
    for i in range(n_trans):
        idx_PE = np.arange(0, n_readsPE)
        idx_U1 = np.arange(n_readsPE, n_readsPE + n_readsU1)
        idx_U2 = np.arange(n_readsPE + n_readsU1, n_reads)
        
        Rmat[idx_PE, i] = (
            check_reads_compatible(gene.trans[i], reads["reads1"], edge_hang, junc_hang) * 
            check_reads_compatible(gene.trans[i], reads["reads2"], edge_hang, junc_hang)
        )
        Rmat[idx_U1, i] = check_reads_compatible(
            gene.trans[i], reads["reads1u"], edge_hang, junc_hang)
        Rmat[idx_U2, i] = check_reads_compatible(
            gene.trans[i], reads["reads2u"], edge_hang, junc_hang)

    return Rmat, R_CB, R_UR


def encode_reads(Rmat, R_CB, R_UR, cell_list, g_idx, merge_UMIs=True,
        matched_reads_only=False, verbose=False):
    """Encode reads
    """
    ## merge UMIs (sort, merge)
    if merge_UMIs and len(R_UR) > 0:
        # sort
        CB_UMI = [R_CB[i] + R_UR[i] for i in range(len(R_CB))]
        sort_idx = np.argsort(CB_UMI) # Note, categorical may be faster
        
        Rmat = Rmat[sort_idx, :]
        R_CB = [R_CB[x] for x in sort_idx]
        R_UR = [R_UR[x] for x in sort_idx]
        CB_UMI = [CB_UMI[x] for x in sort_idx]

        # merge
        _uniq_idx = []
        _curr_bar = 'NA'
        for i in range(len(CB_UMI)):
            if CB_UMI[i] != _curr_bar:
                _curr_idx = i
                _curr_bar = CB_UMI[i]
                _uniq_idx.append(_curr_idx)
                # print(_curr_bar, "New")
            else:
                Rmat[_curr_idx, :] *= Rmat[i, :]
                # print(CB_UMI[i])

        Rmat = Rmat[_uniq_idx, :]
        R_CB = [R_CB[x] for x in _uniq_idx]
        R_UR = [R_UR[x] for x in _uniq_idx]

        if verbose == True or verbose >= 1:
            print("Merged %d reads into %d UMIs" %(len(CB_UMI), len(_uniq_idx)))

    ## remove reads unmatched to any transcript; to reduce candidate reads
    if matched_reads_only:
        idx = np.where(Rmat.sum(axis = 1) > 0)[0]
        Rmat = Rmat[idx, :]
        R_CB = [R_CB[i] for i in idx]
        R_UR = [R_UR[i] for i in idx]

    ## Encode read compatibility & cell barcodes
    K = 2**(np.arange(Rmat.shape[1]))
    _R_code = np.dot(Rmat, K)
    _CB_ids = match(R_CB, cell_list, uniq_ref_only=False)

    ## remove unlisted barcodes
    # print(np.mean(_CB_ids == None))
    _R_code = _R_code[_CB_ids != None]
    _CB_ids = _CB_ids[_CB_ids != None].astype(float).astype(int)
    
    ## sort _CB_ids
    _idx_sort = np.argsort(_CB_ids)
    _CB_ids = _CB_ids[_idx_sort]
    _R_code = _R_code[_idx_sort]
    _CB_uniq, _CB_uniq_idx = np.unique(_CB_ids, return_index=True)

    RV = []
    for c in range(len(_CB_uniq)):
        i1 = _CB_uniq_idx[c]
        if c < len(_CB_uniq) - 1:
            i2 = _CB_uniq_idx[c+ 1]
        else:
            i2 = len(_CB_ids)

        code_id, code_cnt = np.unique(_R_code[i1:i2], return_counts=True)
        count_dict = {}
        for i in range(len(code_id)):
            count_dict["%d" %(code_id[i])] = code_cnt[i]

        RV_line = "%d\t%d\t%s\n" %(_CB_uniq[c] + 1, g_idx + 1, str(count_dict))
        RV.append(RV_line)

    return RV


def _count_one_gene(sam_file, genes, g, cell_list, event_type="SE", 
    edge_hang=10, junc_hang=2, CB_tag='CB', UMI_tag='UR', verbose=False):
    """Counting UMIs for all cells on one gene
    """
    # Load bam file
    samFile, _chrom = check_pysam_chrom(sam_file, genes[g].chrom)
    genes[g].chrom = _chrom

    if verbose == True or verbose >= 1:
        print("")
        print("[BRIE2] parsing gene %d: %s, %s" 
            %(g + 1, genes[g].geneName, genes[g].geneID))
        print("[BRIE2] transcript lengths:", [x.tranL for x in genes[g].trans])

    _Rmat, _R_CB, _R_UR = get_droplet_UMIcount(genes[g], samFile, 
        event_type, edge_hang, junc_hang, CB_tag, UMI_tag, 
        rm_duplicate=True, inner_only=False, mapq_min=0, trimLen_max=15, 
        rlen_min=1, is_mated=True, verbose=verbose)

    if _Rmat.shape[0] == 0:
        RV = None
    else:
        RV = encode_reads(_Rmat, _R_CB, _R_UR, cell_list, g, verbose)
    
    return RV


def get_droplet_matrix(genes, sam_file, cell_list, out_dir, event_type="SE", 
        edge_hang=10, junc_hang=2, CB_tag='CB', UMI_tag='UR', nproc=1, 
        verbose=False):
    """Fetch UMI count matrix for droplet based scRNA-seq data
    Note, trimLen_max is 15 here; different from get_count_matrix with 5.
    """
    global START_TIME, PROCESSED, TOTAL_GENE, FID
    FID = None
    PROCESSED = 0
    TOTAL_GENE = len(genes)
    START_TIME = time.time()
        
    def _show_progress(RV=None):    
        global PROCESSED, TOTAL_GENE, START_TIME, FID
        if RV is not None: 
            FID.writelines(RV)
            
            PROCESSED += 1
            bar_len = 20
            run_time = time.time() - START_TIME
            percents = 100.0 * PROCESSED / TOTAL_GENE
            filled_len = int(bar_len * percents / 100)
            bar = '=' * filled_len + '-' * (bar_len - filled_len)
            
            sys.stdout.write('\r[BRIE2] [%s] %.1f%% genes done in %.1f sec.' 
                % (bar, percents, run_time))
            sys.stdout.flush()
        return RV


    FID = open(out_dir + "/read_count.mtx", "w")
    FID.writelines("%" + "%MatrixMarket matrix coordinate integer general\n")
    FID.writelines("%d\t%d\t%d\n" %(cell_list.shape[0], len(genes), 0))

    # processing each gene with multiple processors
    if nproc <= 1:
        for g in range(len(genes)):
            res = _count_one_gene(sam_file, genes, g, cell_list, event_type, 
                edge_hang, junc_hang, CB_tag, UMI_tag, verbose)
            _show_progress(res)
    else:
        pool = multiprocessing.Pool(processes=nproc)
        result = []
        for g in range(len(genes)):
            result.append(pool.apply_async(_count_one_gene, 
                (sam_file, genes, g, cell_list, event_type, 
                 edge_hang, junc_hang, CB_tag, UMI_tag, verbose), 
                callback=_show_progress))
        pool.close()
        pool.join()
    
    FID.close()

    print("")
    print("[BRIE2] %d genes have been processed." %(len(genes)))
    return None
