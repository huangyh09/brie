# This module is to parse the sam/bam file with the pysam tool, 
# which can be found at: http://pysam.readthedocs.org

import sys
import pysam

def load_samfile(sam_file):
    """To load an indexed bam file"""
    ftype = sam_file.split(".")[-1]
    if ftype != "bam" and ftype != "sam":
        print("Error: file type need suffix of bam or sam.")
        sys.exit(1)
    # print("Loading a %s" %ftype + " with file name: %s" %sam_file)
    if ftype == "bam":
        samfile = pysam.Samfile(sam_file, "rb")
    else:
        samfile = pysam.Samfile(sam_file, "r")
    return samfile


def fetch_reads(samfile, chrom, start, end, rm_duplicate=True, inner_only=True,
                mapq_min=0, mismatch_max=10**5, rlen_min=1,  is_mated=True):
    """To fetch the reads in a given region from a pysam AlignmentFile.

    Args:
        samfile: A Samfile object in pysam.
        chrom: A string of chromosome, e.g., "IV", "chr10".
        start: An integer of the start position for mapped reads.
        end: An integer of the end position for mapped reads.
        rm_duplicate: A bool for only keeping the first one of duplicates.
        inner_only: A bool for only keeping fully region matched reads.
        mapq_min: An integer of the minimum of map quality.
        mismatch_min: An integer of the minimum of mismatched bases.
        rlen_min: An integer of the minimum of read length.
        is_paired: A bool for mating paired-end reads.

    Returns:
        A dict containing lists of mated reads1 and reads2, and unmated reads1u
        and reads2u, i.e.,
        {'reads1': [r11, r21, ...]
         'reads2': [r12, r22, ...]
         'reads1u': [r*1, r*1, ...]
         'reads2u': [r*2, r*2, ...]}
        reads1 is the 5-end of the fragment, and reads2 is the 3-end of the 
        fragment.

    Raises:
        ValueError: An error occurred when fetching reads.
        AssertionError: An error occurred when fetching reads.
    """
    #part 1. check the input and fetch the reads
    chrom  = str(chrom)
    if chrom in samfile.references:
        pass
    else:
        chrom_parts = chrom.split("chr")
        if len(chrom_parts) <= 1:
            chrom = chrom_parts[0]
        else:
            chrom = chrom_parts[1]
    try:
        reads = samfile.fetch(chrom, start, end)

    except ValueError:
        print("Cannot fetch reads in region: %s:%d-%d" %(chrom, start, end))
    except AssertionError:
        print("AssertionError in region: %s:%d-%d" %(chrom, start, end))
        print(" - Check that your BAM file is indexed!")

    #part 2. get reads and filter some of them
    qname1, qname2 = [], []
    reads1, reads2 = [], []
    r_prev = None
    for r in reads:
        # filter 4: only keep the first one of duplicates
        if (rm_duplicate and r_prev is not None and r_prev.qname == r.qname and 
            r_prev.positions == r.positions): r_prev = r; continue
        r_prev = r
        # filter 1: only particially mapped to the regions
        if inner_only == True and (r.pos < start or r.aend > end): continue
        # filter 2: too low map quality
        if r.mapq < mapq_min: continue
        # filter 3: too long mismatch
        if r.rlen - len(r.positions) > mismatch_max: continue
        # filter 5: too short reads length
        if r.rlen < rlen_min: continue
        
        if r.is_read2:
            reads2.append(r)
            qname2.append(r.qname)
        else:
            reads1.append(r)
            qname1.append(r.qname)

    #part 2.1 chech the mate reads' query
    FLAG = True
    if len(qname1) > 0:
        for i in range(len(qname1)-1):
            if qname1[i][-1] != qname1[i+1][-1]:
                FLAG = False
                break
    if FLAG and len(qname2) > 0:
        for i in range(len(qname1)-1):
            if qname1[i][-1] != qname1[i+1][-1]:
                FLAG = False
                break

    if FLAG:
        for i in range(len(qname1)):
            qname1[i] = qname1[i][:-1]
        for i in range(len(qname2)):
            qname2[i] = qname2[i][:-1]


    # part 3. mate the reads
    rv_reads1, rv_reads2 = [], []
    rv_reads1u, rv_reads2u = [], []
    if is_mated == True:
        idx1 = sorted(range(len(qname1)), key=qname1.__getitem__)
        idx2 = sorted(range(len(qname2)), key=qname2.__getitem__)

        i1, i2 = 0, 0
        while i1 < len(idx1) and i2 < len(idx2):
            if qname1[idx1[i1]] == qname2[idx2[i2]]:
                rv_reads1.append(reads1[idx1[i1]])
                rv_reads2.append(reads2[idx2[i2]])
                i1, i2 = i1 + 1, i2 + 1
            elif qname1[idx1[i1]] < qname2[idx2[i2]]:
                rv_reads1u.append(reads1[idx1[i1]])
                i1 += 1
            elif qname1[idx1[i1]] > qname2[idx2[i2]]:
                rv_reads2u.append(reads2[idx2[i2]])
                i2 += 1
        for i in range(i1, len(idx1)):
            rv_reads1u.append(reads1[idx1[i]])
        for i in range(i2, len(idx2)):
            rv_reads2u.append(reads2[idx2[i]])
    else:
        rv_reads1u, rv_reads2u = reads1, reads2

    # part 4. return reads
    RV = {}
    RV["reads1"] = rv_reads1
    RV["reads2"] = rv_reads2
    RV["reads1u"] = rv_reads1u
    RV["reads2u"] = rv_reads2u
    return RV
    
