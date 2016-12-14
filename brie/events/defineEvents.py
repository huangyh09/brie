# Define alternative splicing events
#
# The 5 core functions in this file is almost identical
# as Yarden's original files, which may need some check,
# particularly for RI, as the original file used another
# function.

import os, operator, string

LETTERS = string.uppercase


def A3SS(DtoA_F, AtoD_F, DtoA_R, AtoD_R, gff3_f,
         flanking='commonshortest', multi_iso=False):
    """
    Define alt. 3' splice sites
    A3SS are events where a donor is spliced to >1 acceptor, and the 
    acceptors share a downstream donor site in the same exon.
    """
    print "Generating alternative 3\' splice sites (A3SS)"
    if os.path.isfile(gff3_f):
        print("  - Found file, skipping...")
        return
    out = open(gff3_f, 'w')
 
    for donor in DtoA_F:
        # Check to see whether this donor splice site has more than one acceptor
        # site.
        if len(DtoA_F[donor]) > 1:  # if it does, then see if any of them share a downstream 5'ss.
            acceptors = list(DtoA_F[donor])
            nextDonorCounts = {}
            nextDonorToAcceptors = {}   # downstream donor -> list of preious acceptors
            for acceptor in acceptors:
                for x in list(AtoD_F[acceptor].elements()):
                    try:
                        nextDonorCounts[x] += 1
                    except:
                        nextDonorCounts[x] = 1
                for x in list(AtoD_F[acceptor]):
                    try:
                        nextDonorToAcceptors[x].append(acceptor)
                    except:
                        nextDonorToAcceptors[x] = [acceptor]

            acceptorlistnameToDonor = {}

            # Find all possible downstream donor sites for this set of alt. acceptor sites
            for nextDonor in nextDonorToAcceptors:
                acceptorlist = nextDonorToAcceptors[nextDonor]
                if len(acceptorlist) > 1:
                    acceptorlist = [[x, int(x.split(":")[1])] for x in acceptorlist]
                    acceptorlist.sort(key=operator.itemgetter(1))
                    acceptorlistname = ";".join([x[0] for x in acceptorlist])
                    try:
                        acceptorlistnameToDonor[acceptorlistname].append(nextDonor)
                    except:
                        acceptorlistnameToDonor[acceptorlistname] = [nextDonor]

            if len(acceptorlistnameToDonor) > 0:
                chrom, coord, strand = donor.split(":")
                donorCoord = donor.split(":")[1]
              
                for acceptorlistname in acceptorlistnameToDonor:
                    acceptorCoords = [x.split(":")[1] for x in acceptorlistname.split(";")]

                    # Get the previous acceptor sites and their frequencies
                    prevAcceptorCounts = DtoA_R[donor]
                    prevAcceptors = [int(x.split(":")[1]) for x in list(DtoA_R[donor])]
                    nextDonorDict = {}
                    for nextDonor in acceptorlistnameToDonor[acceptorlistname]:
                        nextDonorDict[int(nextDonor.split(":")[1])] = \
                            nextDonorCounts[nextDonor]
                    nextDonors = nextDonorDict.keys()

                    if flanking == 'shortest': 
                        if strand == '+': 
                            prevAcceptorCoord = str(sorted(prevAcceptors)[-1])
                            nextDonorCoord = str(sorted(nextDonors)[0])
                        else:
                            prevAcceptorCoord = str(sorted(prevAcceptors)[0])
                            nextDonorCoord = str(sorted(nextDonors)[-1])
                    elif flanking == 'longest':
                        if strand == '+': 
                            prevAcceptorCoord = str(sorted(prevAcceptors)[0])
                            nextDonorCoord = str(sorted(nextDonors)[-1])
                        else:
                            prevAcceptorCoord = str(sorted(prevAcceptors)[-1])
                            nextDonorCoord = str(sorted(nextDonors)[0])
                    elif flanking in ['commonshortest', 'commonlongest']:
                        prevAcceptorOptions = [int(x.split(":")[1]) for x in prevAcceptorCounts \
                            if prevAcceptorCounts[x] == max(prevAcceptorCounts.values())]
                        nextDonorOptions = [x for x in nextDonorDict \
                            if nextDonorDict[x] == max(nextDonorDict.values())]
                        if flanking == 'commonshortest':
                            if strand == '+':
                                prevAcceptorCoord = str(sorted(prevAcceptorOptions)[0])
                                nextDonorCoord = str(sorted(nextDonorOptions)[-1])
                            else: 
                                prevAcceptorCoord = str(sorted(prevAcceptorOptions)[-1])
                                nextDonorCoord = str(sorted(nextDonorOptions)[0])
                        else:
                            if strand == '+':
                                prevAcceptorCoord = str(sorted(prevAcceptorOptions)[-1])
                                nextDonorCoord = str(sorted(nextDonorOptions)[0])
                            else: 
                                prevAcceptorCoord = str(sorted(prevAcceptorOptions)[0])
                                nextDonorCoord = str(sorted(nextDonorOptions)[-1])

                    # Output output 2-iso or multiple isoforms
                    if multi_iso:
                        if strand == '+':
                            upexon = ":".join([chrom, prevAcceptorCoord, donorCoord, strand])
                            altexon = ":".join([chrom, "|".join(acceptorCoords),\
                                nextDonorCoord, strand])
                            name = "@".join([upexon, altexon]) 
                            out.write("\t".join([chrom, 'A3SS', 'gene', prevAcceptorCoord,\
                                nextDonorCoord, '.', strand, '.',\
                                "ID=" + name + ";Name=" + name]) + "\n")
                            for i in range(len(acceptorCoords)):
                                out.write("\t".join([chrom, 'A3SS', 'mRNA', prevAcceptorCoord,\
                                    nextDonorCoord, '.', strand, '.',\
                                    "ID=" + name + "." + LETTERS[i] + ";Parent=" + name]) + "\n")
                            for i in range(len(acceptorCoords)): 
                                out.write("\t".join([chrom, 'A3SS', 'exon', prevAcceptorCoord,\
                                    donorCoord, '.', strand, '.',\
                                    "ID=" + name + "." + LETTERS[i] + ".up;Parent=" + name + "." +\
                                    LETTERS[i]]) + "\n")
                                out.write("\t".join([chrom, 'A3SS', 'exon', acceptorCoords[i],\
                                    nextDonorCoord, '.', strand, '.',\
                                    "ID=" + name + "." + LETTERS[i] + ".altss;Parent=" + name + "." +\
                                    LETTERS[i]]) + "\n")
                        else:
                            upexon = ":".join([chrom, donorCoord, prevAcceptorCoord, strand])
                            altexon = ":".join([chrom, "|".join(acceptorCoords),\
                                nextDonorCoord, strand])
                            name = "@".join([upexon, altexon]) 
                            out.write("\t".join([chrom, 'A3SS', 'gene', nextDonorCoord,\
                                prevAcceptorCoord, '.', strand, '.',\
                                "ID=" + name + ";Name=" + name]) + "\n")
                            for i in range(len(acceptorCoords)):
                                out.write("\t".join([chrom, 'A3SS', 'mRNA', nextDonorCoord,\
                                    prevAcceptorCoord, '.', strand, '.',\
                                    "ID=" + name + "." + LETTERS[i] + ";Parent=" + name]) + "\n")
                            for i in range(len(acceptorCoords)): 
                                out.write("\t".join([chrom, 'A3SS', 'exon', donorCoord,\
                                    prevAcceptorCoord, '.', strand, '.',\
                                    "ID=" + name + "." + LETTERS[i] + ".up;Parent=" + name + "." +\
                                    LETTERS[i]]) + "\n")
                                out.write("\t".join([chrom, 'A3SS', 'exon', nextDonorCoord,\
                                    acceptorCoords[i], '.', strand, '.',\
                                    "ID=" + name + "." + LETTERS[i] + ".altss;Parent=" + name + "." +\
                                    LETTERS[i]]) + "\n")
    
                    # Otherwise do 2-isoform
                    else:
                        for iso1 in range(len(acceptorCoords)):
                            for iso2 in range(iso1 + 1, len(acceptorCoords)):
                                if strand == '+':
                                    upexon = ":".join([chrom, prevAcceptorCoord, donorCoord, strand])
                                    altexon = ":".join([chrom, "|".join([acceptorCoords[iso1],\
                                        acceptorCoords[iso2]]), nextDonorCoord, strand])
                                    name = "@".join([upexon, altexon]) 
                                    out.write("\t".join([chrom, 'A3SS', 'gene', prevAcceptorCoord,\
                                        nextDonorCoord, '.', strand, '.',\
                                        "ID=" + name + ";Name=" + name]) + "\n")
                                    out.write("\t".join([chrom, 'A3SS', 'mRNA', prevAcceptorCoord,\
                                        nextDonorCoord, '.', strand, '.',\
                                        "ID=" + name + ".A;Parent=" + name]) + "\n")
                                    out.write("\t".join([chrom, 'A3SS', 'exon', prevAcceptorCoord,\
                                        donorCoord, '.', strand, '.',\
                                        "ID=" + name + ".A.up;Parent=" + name + ".A"]) + "\n")
                                    out.write("\t".join([chrom, 'A3SS', 'exon', acceptorCoords[iso2],\
                                        nextDonorCoord, '.', strand, '.',\
                                        "ID=" + name + ".A.coreAndExt;Parent=" + name + ".A"]) + "\n")
                                    out.write("\t".join([chrom, 'A3SS', 'mRNA', prevAcceptorCoord,\
                                        nextDonorCoord, '.', strand, '.',\
                                        "ID=" + name + ".B;Parent=" + name]) + "\n")
                                    out.write("\t".join([chrom, 'A3SS', 'exon', prevAcceptorCoord,\
                                        donorCoord, '.', strand, '.',\
                                        "ID=" + name + ".B.up;Parent=" + name + ".B"]) + "\n")
                                    out.write("\t".join([chrom, 'A3SS', 'exon', acceptorCoords[iso1],\
                                        nextDonorCoord, '.', strand, '.',\
                                        "ID=" + name + ".B.core;Parent=" + name + ".B"]) + "\n")
                                else:
                                    upexon = ":".join([chrom, donorCoord, prevAcceptorCoord, strand])
                                    altexon = ":".join([chrom, "|".join([acceptorCoords[iso1],\
                                        acceptorCoords[iso2]]), nextDonorCoord, strand])
                                    name = "@".join([upexon, altexon]) 
                                    out.write("\t".join([chrom, 'A3SS', 'gene', nextDonorCoord,\
                                        prevAcceptorCoord, '.', strand, '.',\
                                        "ID=" + name + ";Name=" + name]) + "\n")
                                    out.write("\t".join([chrom, 'A3SS', 'mRNA', nextDonorCoord,\
                                        prevAcceptorCoord, '.', strand, '.',\
                                        "ID=" + name + ".A;Parent=" + name]) + "\n")
                                    out.write("\t".join([chrom, 'A3SS', 'exon', donorCoord,\
                                        prevAcceptorCoord, '.', strand, '.',\
                                        "ID=" + name + ".A.up;Parent=" + name + ".A"]) + "\n")
                                    out.write("\t".join([chrom, 'A3SS', 'exon', nextDonorCoord,\
                                        acceptorCoords[iso2], '.', strand, '.',\
                                        "ID=" + name + ".A.coreAndExt;Parent=" + name + ".A"]) + "\n")
                                    out.write("\t".join([chrom, 'A3SS', 'mRNA', nextDonorCoord,\
                                        prevAcceptorCoord, '.', strand, '.',\
                                        "ID=" + name + ".B;Parent=" + name]) + "\n")
                                    out.write("\t".join([chrom, 'A3SS', 'exon', donorCoord,\
                                        prevAcceptorCoord, '.', strand, '.',\
                                        "ID=" + name + ".B.up;Parent=" + name + ".B"]) + "\n")
                                    out.write("\t".join([chrom, 'A3SS', 'exon', nextDonorCoord,\
                                        acceptorCoords[iso1], '.', strand, '.',\
                                        "ID=" + name + ".B.core;Parent=" + name + ".B"]) + "\n")

    out.close()



def A5SS(DtoA_F, AtoD_F, DtoA_R, AtoD_R, gff3_f,
         flanking='commonshortest', multi_iso=False):
    """
    Define alt. 5' splice sites
    A3SS are events where >1 donor is spliced to an acceptor, and the 
    donors share an upstream acceptor site in the same exon.
    """
    print "Generating alternative 5\' splice sites (A5SS)"
    if os.path.isfile(gff3_f):
        print("  - Found file, skipping...")
        return
    out = open(gff3_f, 'w')
 
    for acceptor in AtoD_R:
        # Check to see whether this acceptor splice site has more than one donor 
        # site.
        if len(AtoD_R[acceptor]) > 1:   # if it does, then see if any of them share an upstream 3'ss.
            donors = list(AtoD_R[acceptor])     # these are the alt 5'ss 
            prevAcceptorCounts = {}             # coordinate -> counts
            prevAcceptorToDonors = {}           # splice site -> list of donors

            # Iterate through each alt donor and get their upstream acceptors
            for donor in donors:
                for x in list(DtoA_R[donor].elements()):
                    try:
                        prevAcceptorCounts[x] += 1
                    except:
                        prevAcceptorCounts[x] = 1
                for x in list(DtoA_R[donor]):
                    try:
                        prevAcceptorToDonors[x].append(donor)
                    except:
                        prevAcceptorToDonors[x] = [donor]

            donorlistnameToAcceptor = {}

            # Find all possible acceptor sites for this list of alt. donor sites
            for prevAcceptor in prevAcceptorToDonors:
                donorlist = prevAcceptorToDonors[prevAcceptor]
                if len(donorlist) > 1:                  # this is an alt. 5'ss event
                    donorlist = [[x, int(x.split(":")[1])] for x in donorlist]
                    donorlist.sort(key=operator.itemgetter(1))
                    donorlistname = ";".join([x[0] for x in donorlist])
                    try:
                        donorlistnameToAcceptor[donorlistname].append(prevAcceptor)
                    except:
                        donorlistnameToAcceptor[donorlistname] = [prevAcceptor]
            
            if len(donorlistnameToAcceptor) > 0: 
                chrom, coord, strand = acceptor.split(":")
                acceptorCoord = acceptor.split(":")[1]

                for donorlistname in donorlistnameToAcceptor:
                    donorCoords = [x.split(":")[1] for x in donorlistname.split(";")] 

                    # Get the next donor sites and their frequencies
                    nextDonorCounts = AtoD_F[acceptor]
                    nextDonors = [int(x.split(":")[1]) for x in list(AtoD_F[acceptor])]
                    prevAcceptorDict = {}
                    for prevAcceptor in donorlistnameToAcceptor[donorlistname]:
                        prevAcceptorDict[int(prevAcceptor.split(":")[1])] = \
                            prevAcceptorCounts[prevAcceptor]
                    prevAcceptors = prevAcceptorDict.keys()

                    if flanking == 'shortest': 
                        if strand == '+': 
                            prevAcceptorCoord = str(sorted(prevAcceptors)[-1])
                            nextDonorCoord = str(sorted(nextDonors)[0])
                        else:
                            prevAcceptorCoord = str(sorted(prevAcceptors)[0])
                            nextDonorCoord = str(sorted(nextDonors)[-1])
                    elif flanking == 'longest':
                        if strand == '+': 
                            prevAcceptorCoord = str(sorted(prevAcceptors)[0])
                            nextDonorCoord = str(sorted(nextDonors)[-1])
                        else:
                            prevAcceptorCoord = str(sorted(prevAcceptors)[-1])
                            nextDonorCoord = str(sorted(nextDonors)[0])
                    elif flanking in ['commonshortest', 'commonlongest']:
                        prevAcceptorOptions = [x for x in prevAcceptorDict \
                            if prevAcceptorDict[x] == max(prevAcceptorDict.values())]
                        nextDonorOptions = [int(x.split(":")[1]) for x in nextDonorCounts \
                            if nextDonorCounts[x] == max(nextDonorCounts.values())]
                        if flanking == 'commonshortest':
                            if strand == '+':
                                prevAcceptorCoord = str(sorted(prevAcceptorOptions)[0])
                                nextDonorCoord = str(sorted(nextDonorOptions)[-1])
                            else: 
                                prevAcceptorCoord = str(sorted(prevAcceptorOptions)[-1])
                                nextDonorCoord = str(sorted(nextDonorOptions)[0])
                        else:
                            if strand == '+':
                                prevAcceptorCoord = str(sorted(prevAcceptorOptions)[-1])
                                nextDonorCoord = str(sorted(nextDonorOptions)[0])
                            else: 
                                prevAcceptorCoord = str(sorted(prevAcceptorOptions)[0])
                                nextDonorCoord = str(sorted(nextDonorOptions)[-1])

                    # Output multiple isoforms if necessary
                    if multi_iso:
                        if strand == '+':
                            altexon = ":".join([chrom, prevAcceptorCoord,\
                                "|".join(donorCoords), strand])
                            dnexon = ":".join([chrom, acceptorCoord, nextDonorCoord, strand])
                            name = "@".join([altexon, dnexon]) 
                            out.write("\t".join([chrom, 'A5SS', 'gene', prevAcceptorCoord,\
                                nextDonorCoord, '.', strand, '.',\
                                "ID=" + name + ";Name=" + name]) + "\n")
                            for i in range(len(donorCoords)):
                                out.write("\t".join([chrom, 'A5SS', 'mRNA', prevAcceptorCoord,\
                                    nextDonorCoord, '.', strand, '.',\
                                    "ID=" + name + "." + LETTERS[i] + ";Parent=" + name]) + "\n")
                            for i in range(len(donorCoords)): 
                                out.write("\t".join([chrom, 'A5SS', 'exon', prevAcceptorCoord,\
                                    donorCoords[i], '.', strand, '.',\
                                    "ID=" + name + "." + LETTERS[i] + ".altss;Parent=" + name + "." +\
                                    LETTERS[i]]) + "\n")
                                out.write("\t".join([chrom, 'A5SS', 'exon', acceptorCoord,\
                                    nextDonorCoord, '.', strand, '.',\
                                    "ID=" + name + "." + LETTERS[i] + ".dn;Parent=" + name + "." +\
                                    LETTERS[i]]) + "\n")
                        else:
                            altexon = ":".join([chrom, prevAcceptorCoord,\
                                "|".join(donorCoords), strand])
                            dnexon = ":".join([chrom, nextDonorCoord, acceptorCoord, strand])
                            name = "@".join([altexon, dnexon]) 
                            out.write("\t".join([chrom, 'A5SS', 'gene', nextDonorCoord,\
                                prevAcceptorCoord, '.', strand, '.',\
                                "ID=" + name + ";Name=" + name]) + "\n")
                            for i in range(len(donorCoords)):
                                out.write("\t".join([chrom, 'A5SS', 'mRNA', nextDonorCoord,\
                                    prevAcceptorCoord, '.', strand, '.',\
                                    "ID=" + name + "." + LETTERS[i] + ";Parent=" + name]) + "\n")
                            for i in range(len(donorCoords)): 
                                out.write("\t".join([chrom, 'A5SS', 'exon', donorCoords[i],\
                                    prevAcceptorCoord, '.', strand, '.',\
                                    "ID=" + name + "." + LETTERS[i] + ".altss;Parent=" + name + "." +\
                                    LETTERS[i]]) + "\n")
                                out.write("\t".join([chrom, 'A5SS', 'exon', nextDonorCoord,\
                                    acceptorCoord, '.', strand, '.',\
                                    "ID=" + name + "." + LETTERS[i] + ".dn;Parent=" + name + "." +\
                                    LETTERS[i]]) + "\n")
                    
                    # Otherwise do 2 isoform.
                    else:
                        for iso1 in range(len(donorCoords)):
                            for iso2 in range(iso1 + 1, len(donorCoords)):
                                if strand == '+':
                                    altexon = ":".join([chrom, prevAcceptorCoord,\
                                        "|".join([donorCoords[iso1], donorCoords[iso2]]),\
                                        strand])
                                    dnexon = ":".join([chrom, acceptorCoord, nextDonorCoord, strand])
                                    name = "@".join([altexon, dnexon]) 
                                    out.write("\t".join([chrom, 'A5SS', 'gene', prevAcceptorCoord,\
                                        nextDonorCoord, '.', strand, '.',\
                                        "ID=" + name + ";Name=" + name]) + "\n")
                                    out.write("\t".join([chrom, 'A5SS', 'mRNA', prevAcceptorCoord,\
                                        nextDonorCoord, '.', strand, '.',\
                                        "ID=" + name + ".A;Parent=" + name]) + "\n")
                                    out.write("\t".join([chrom, 'A5SS', 'exon', prevAcceptorCoord,\
                                        donorCoords[iso2], '.', strand, '.',\
                                        "ID=" + name + ".A.coreAndExt;Parent=" + name + ".A"]) + "\n")
                                    out.write("\t".join([chrom, 'A5SS', 'exon', acceptorCoord,\
                                        nextDonorCoord, '.', strand, '.',\
                                        "ID=" + name + ".A.dn;Parent=" + name + ".A"]) + "\n")
                                    out.write("\t".join([chrom, 'A5SS', 'mRNA', prevAcceptorCoord,\
                                        nextDonorCoord, '.', strand, '.',\
                                        "ID=" + name + ".B;Parent=" + name]) + "\n")
                                    out.write("\t".join([chrom, 'A5SS', 'exon', prevAcceptorCoord,\
                                        donorCoords[iso1], '.', strand, '.',\
                                        "ID=" + name + ".B.core;Parent=" + name + ".B"]) + "\n")
                                    out.write("\t".join([chrom, 'A5SS', 'exon', acceptorCoord,\
                                        nextDonorCoord, '.', strand, '.',\
                                        "ID=" + name + ".B.dn;Parent=" + name + ".B"]) + "\n")
                                else:
                                    altexon = ":".join([chrom, prevAcceptorCoord,\
                                        "|".join([donorCoords[iso1], donorCoords[iso2]]),\
                                        strand])
                                    dnexon = ":".join([chrom, nextDonorCoord, acceptorCoord, strand])
                                    name = "@".join([altexon, dnexon]) 
                                    out.write("\t".join([chrom, 'A5SS', 'gene', nextDonorCoord,\
                                        prevAcceptorCoord, '.', strand, '.',\
                                        "ID=" + name + ";Name=" + name]) + "\n")
                                    out.write("\t".join([chrom, 'A5SS', 'mRNA', nextDonorCoord,\
                                        prevAcceptorCoord, '.', strand, '.',\
                                        "ID=" + name + ".A;Parent=" + name]) + "\n")
                                    out.write("\t".join([chrom, 'A5SS', 'exon', donorCoords[iso1],\
                                        prevAcceptorCoord, '.', strand, '.',\
                                        "ID=" + name + ".A.coreAndExt;Parent=" + name + ".A"]) + "\n")
                                    out.write("\t".join([chrom, 'A5SS', 'exon', nextDonorCoord,\
                                        acceptorCoord, '.', strand, '.',\
                                        "ID=" + name + ".A.dn;Parent=" + name + ".A"]) + "\n")
                                    out.write("\t".join([chrom, 'A5SS', 'mRNA', nextDonorCoord,\
                                        prevAcceptorCoord, '.', strand, '.',\
                                        "ID=" + name + ".B;Parent=" + name]) + "\n")
                                    out.write("\t".join([chrom, 'A5SS', 'exon', donorCoords[iso2],\
                                        prevAcceptorCoord, '.', strand, '.',\
                                        "ID=" + name + ".B.core;Parent=" + name + ".B"]) + "\n")
                                    out.write("\t".join([chrom, 'A5SS', 'exon', nextDonorCoord,\
                                        acceptorCoord, '.', strand, '.',\
                                        "ID=" + name + ".B.dn;Parent=" + name + ".B"]) + "\n")
    out.close()


def SE(DtoA_F, AtoD_F, DtoA_R, AtoD_R, gff3_f,
       flanking='commonshortest'):
    """
    Define skipped exons.

    Arguments are:
    splice site dictionaries, followed by output file and a keyword to 
    denote method of selecting flanking exon coordinates:
      - shortest
      - longest
      - commonshortest
      - commonlongest
    """
    print "Generating skipped exons (SE)"
    if os.path.isfile(gff3_f):
        print("  - Found file, skipping...")
        return
    out = open(gff3_f, 'w')
 
    for acceptor in AtoD_R:                             # this acceptor is the SE acceptor
        chrom, coord, strand = acceptor.split(":") 
        donors = list(AtoD_F[acceptor])                 # these are possible 5'ss of the SE
        for donor in donors:                            # consider each SE
            upDonors = list(AtoD_R[acceptor])           # these are 5'ss that splice to SE 
            if donor in DtoA_F:
                dnAcceptors = list(DtoA_F[donor])       # these are 3'ss that SE splices to
                for upDonor in upDonors:
                    for dnAcceptor in dnAcceptors:
                        if dnAcceptor in list(DtoA_F[upDonor]):   # this is an SE
                            ss1list = map(int, [x.split(":")[1] for x in list(DtoA_R[upDonor].elements())])
                            ss2 = upDonor.split(":")[1]
                            ss3 = acceptor.split(":")[1]
                            ss4 = donor.split(":")[1]
                            ss5 = dnAcceptor.split(":")[1]
                            ss6list = map(int, [x.split(":")[1] for x in list(AtoD_F[dnAcceptor].elements())])
  
                            if flanking == 'shortest': 
                                if strand == '+': 
                                    ss1 = str(sorted(ss1list)[-1])
                                    ss6 = str(sorted(ss6list)[0])
                                else:
                                    ss1 = str(sorted(ss1list)[0])
                                    ss6 = str(sorted(ss6list)[-1])
                            elif flanking == 'longest':
                                if strand == '+': 
                                    ss1 = str(sorted(ss1list)[0])
                                    ss6 = str(sorted(ss6list)[-1])
                                else:
                                    ss1 = str(sorted(ss1list)[-1])
                                    ss6 = str(sorted(ss6list)[0])
                            elif flanking in ['commonshortest', 'commonlongest']:
                                ss1cts = map(ss1list.count, ss1list)
                                ss1options = [ss1list[x] for x in range(len(ss1cts)) if ss1cts[x] == max(ss1cts)]
                                ss6cts = map(ss6list.count, ss6list)
                                ss6options = [ss6list[x] for x in range(len(ss6cts)) if ss6cts[x] == max(ss6cts)]
                                if flanking == 'commonshortest':
                                    if strand == '+':
                                        ss1 = str(sorted(ss1options)[0])
                                        ss6 = str(sorted(ss6options)[-1])
                                    else: 
                                        ss1 = str(sorted(ss1options)[-1])
                                        ss6 = str(sorted(ss6options)[0])
                                else:
                                    if strand == '+':
                                        ss1 = str(sorted(ss1options)[-1])
                                        ss6 = str(sorted(ss6options)[0])
                                    else: 
                                        ss1 = str(sorted(ss1options)[0])
                                        ss6 = str(sorted(ss6options)[-1])

                            # Iterate through each possible set of flanking exons
                            # for i in range(len(ss1list)):
                            #    for j in range(len(ss6list)):
                            if strand == '+': 
                                upexon = ":".join([chrom, ss1, ss2, strand])
                                seexon = ":".join([chrom, ss3, ss4, strand])
                                dnexon = ":".join([chrom, ss5, ss6, strand])
                                name = "@".join([upexon, seexon, dnexon])
                                out.write("\t".join([chrom, 'SE', 'gene', ss1, ss6,\
                                    '.', strand, '.', "ID=" + name + ";Name=" + name]) + "\n")
                                out.write("\t".join([chrom, 'SE', 'mRNA', ss1, ss6,\
                                    '.', strand, '.', "ID=" + name + ".A;Parent=" + name]) + "\n")
                                out.write("\t".join([chrom, 'SE', 'exon', ss1, ss2,\
                                    '.', strand, '.', "ID=" + name + ".A.up;Parent=" + name + ".A"]) + "\n")
                                out.write("\t".join([chrom, 'SE', 'exon', ss3, ss4,\
                                    '.', strand, '.', "ID=" + name + ".A.se;Parent=" + name + ".A"]) + "\n")
                                out.write("\t".join([chrom, 'SE', 'exon', ss5, ss6,\
                                    '.', strand, '.', "ID=" + name + ".A.dn;Parent=" + name + ".A"]) + "\n")
                                out.write("\t".join([chrom, 'SE', 'mRNA', ss1, ss6,\
                                    '.', strand, '.', "ID=" + name + ".B;Parent=" + name]) + "\n")
                                out.write("\t".join([chrom, 'SE', 'exon', ss1, ss2,\
                                    '.', strand, '.', "ID=" + name + ".B.up;Parent=" + name + ".B"]) + "\n")
                                out.write("\t".join([chrom, 'SE', 'exon', ss5, ss6,\
                                    '.', strand, '.', "ID=" + name + ".B.dn;Parent=" + name + ".B"]) + "\n")
                                
                            else:
                                upexon = ":".join([chrom, ss2, ss1, strand])
                                seexon = ":".join([chrom, ss4, ss3, strand])
                                dnexon = ":".join([chrom, ss6, ss5, strand])
                                name = "@".join([upexon, seexon, dnexon])
                                out.write("\t".join([chrom, 'SE', 'gene', ss6, ss1,\
                                    '.', strand, '.', "ID=" + name + ";Name=" + name]) + "\n")
                                out.write("\t".join([chrom, 'SE', 'mRNA', ss6, ss1,\
                                    '.', strand, '.', "ID=" + name + ".A;Parent=" + name]) + "\n")
                                out.write("\t".join([chrom, 'SE', 'exon', ss2, ss1,\
                                    '.', strand, '.', "ID=" + name + ".A.up;Parent=" + name + ".A"]) + "\n")
                                out.write("\t".join([chrom, 'SE', 'exon', ss4, ss3,\
                                    '.', strand, '.', "ID=" + name + ".A.se;Parent=" + name + ".A"]) + "\n")
                                out.write("\t".join([chrom, 'SE', 'exon', ss6, ss5,\
                                    '.', strand, '.', "ID=" + name + ".A.dn;Parent=" + name + ".A"]) + "\n")
                                out.write("\t".join([chrom, 'SE', 'mRNA', ss6, ss1,\
                                    '.', strand, '.', "ID=" + name + ".B;Parent=" + name]) + "\n")
                                out.write("\t".join([chrom, 'SE', 'exon', ss2, ss1,\
                                    '.', strand, '.', "ID=" + name + ".B.up;Parent=" + name + ".B"]) + "\n")
                                out.write("\t".join([chrom, 'SE', 'exon', ss6, ss5,\
                                    '.', strand, '.', "ID=" + name + ".B.dn;Parent=" + name + ".B"]) + "\n")

    out.close()


def MXE(DtoA_F, AtoD_F, DtoA_R, AtoD_R, gff3_f,
       flanking='commonshortest'):
    """
    Define mutually exclusive exons.

    Arguments are:
    splice site dictionaries, followed by output file and a keyword to 
    denote method of selecting flanking exon coordinates:
      - shortest
      - longest
      - commonshortest
      - commonlongest
    """
    print "Generating mutually exclusive exons (MXE)"
    if os.path.isfile(gff3_f):
        print("  - Found file, skipping...")
        return
    out = open(gff3_f, 'w')

    visited = {} 
    for mxe1acceptor in AtoD_R:                             # this acceptor is the acceptor for the MXE1
        visited[mxe1acceptor] = 0
        chrom, coord, strand = mxe1acceptor.split(":") 
        mxe1donors = [x for x in list(AtoD_F[mxe1acceptor]) if x in DtoA_F]
                                                            # these are donors for MXE1
        for mxe1donor in mxe1donors:                        
            upDonors = list(AtoD_R[mxe1acceptor])           # these are the upstream 5'ss that splice to MXE1 
            dnAcceptorsFrom1 = list(DtoA_F[mxe1donor])      # these are possible downstream acceptors from MXE1 
            for upDonor in upDonors:
                mxe2acceptors = [x for x in list(DtoA_F[upDonor]) if x not in visited] 
                                                            # these are potential 3'ss for MXE2
                for mxe2acceptor in mxe2acceptors: 
                    mxe2donors = [x for x in list(AtoD_F[mxe2acceptor]) if x in DtoA_F]
                                                            # these are potential 5'ss for MXE2
                    # proper MXEs will share a common downstream acceptor, from MXE1 donor and MXE2 donor
                    # and also be non-overlapping
                    for mxe2donor in mxe2donors:

                        dnAcceptorsFrom2 = list(DtoA_F[mxe2donor])
                                                            # these are possible downstream acceptors from MXE2
                        dnAcceptors = set.intersection(set(dnAcceptorsFrom1), set(dnAcceptorsFrom2))
                        if len(dnAcceptors) > 0:
                            for dnAcceptor in dnAcceptors:

                                # Put MXEs in strand order
                                acceptorcoords = map(int, [mxe1acceptor.split(":")[1], mxe2acceptor.split(":")[1]])
                                donorcoords = map(int, [mxe1donor.split(":")[1], mxe2donor.split(":")[1]])
                                if strand == '+':
                                    if acceptorcoords[1] < acceptorcoords[0]:
                                        acceptorcoords = acceptorcoords[::-1]     # sort by coordinate
                                        donorcoords = donorcoords[::-1]
                                else:
                                    if acceptorcoords[0] < acceptorcoords[1]:     # sort by coordinate
                                        acceptorcoords = acceptorcoords[::-1]
                                        donorcoords = donorcoords[::-1]

                                # Make sure MXEs are non-overlapping
                                if (strand == '+' and donorcoords[0] < acceptorcoords[1]) or \
                                    (strand == '-' and acceptorcoords[1] < donorcoords[0]):

                                    # Output these MXEs
                                    ss1list = map(int, [x.split(":")[1] for x in list(DtoA_R[upDonor].elements())])
                                    ss2 = upDonor.split(":")[1]
                                    ss3 = str(acceptorcoords[0])
                                    ss4 = str(donorcoords[0])
                                    ss5 = str(acceptorcoords[1])
                                    ss6 = str(donorcoords[1])
                                    ss7 = dnAcceptor.split(":")[1]
                                    ss8list = map(int, [x.split(":")[1] for x in list(AtoD_F[dnAcceptor].elements())])

                                    if flanking == 'shortest': 
                                        if strand == '+': 
                                            ss1 = str(sorted(ss1list)[-1])
                                            ss8 = str(sorted(ss8list)[0])
                                        else:
                                            ss1 = str(sorted(ss1list)[0])
                                            ss8 = str(sorted(ss8list)[-1])
                                    elif flanking == 'longest':
                                        if strand == '+': 
                                            ss1 = str(sorted(ss1list)[0])
                                            ss8 = str(sorted(ss8list)[-1])
                                        else:
                                            ss1 = str(sorted(ss1list)[-1])
                                            ss8 = str(sorted(ss8list)[0])
                                    elif flanking in ['commonshortest', 'commonlongest']:
                                        ss1cts = map(ss1list.count, ss1list)
                                        ss1options = [ss1list[x] for x in range(len(ss1cts)) \
                                            if ss1cts[x] == max(ss1cts)]
                                        ss8cts = map(ss8list.count, ss8list)
                                        ss8options = [ss8list[x] for x in range(len(ss8cts)) \
                                            if ss8cts[x] == max(ss8cts)]
                                        if flanking == 'commonshortest':
                                            if strand == '+':
                                                ss1 = str(sorted(ss1options)[0])
                                                ss8 = str(sorted(ss8options)[-1])
                                            else: 
                                                ss1 = str(sorted(ss1options)[-1])
                                                ss8 = str(sorted(ss8options)[0])
                                        else:
                                            if strand == '+':
                                                ss1 = str(sorted(ss1options)[-1])
                                                ss8 = str(sorted(ss8options)[0])
                                            else: 
                                                ss1 = str(sorted(ss1options)[0])
                                                ss8 = str(sorted(ss8options)[-1])

                                        # Iterate through each possible set of flanking exons
                                        # for i in range(len(ss1list)):
                                        #    for j in range(len(ss6list)):
                                        if strand == '+': 
                                            upexon = ":".join([chrom, ss1, ss2, strand])
                                            mxe1 = ":".join([chrom, ss3, ss4, strand])
                                            mxe2 = ":".join([chrom, ss5, ss6, strand])
                                            dnexon = ":".join([chrom, ss7, ss8, strand])
                                            name = "@".join([upexon, mxe1, mxe2, dnexon])
                                            out.write("\t".join([chrom, 'MXE', 'gene', ss1, ss8,\
                                                '.', strand, '.', "ID=" + name + ";Name=" + name]) + "\n")
                                            out.write("\t".join([chrom, 'MXE', 'mRNA', ss1, ss8,\
                                                '.', strand, '.', "ID=" + name + ".A;Parent=" + name]) + "\n")
                                            out.write("\t".join([chrom, 'MXE', 'exon', ss1, ss2,\
                                                '.', strand, '.', "ID=" + name + ".A.up;Parent=" + name + ".A"]) + "\n")
                                            out.write("\t".join([chrom, 'MXE', 'exon', ss3, ss4,\
                                                '.', strand, '.', "ID=" + name + ".A.mxe1;Parent=" + name + ".A"]) + "\n")
                                            out.write("\t".join([chrom, 'MXE', 'exon', ss7, ss8,\
                                                '.', strand, '.', "ID=" + name + ".A.dn;Parent=" + name + ".A"]) + "\n")
                                            out.write("\t".join([chrom, 'MXE', 'mRNA', ss1, ss8,\
                                                '.', strand, '.', "ID=" + name + ".B;Parent=" + name]) + "\n")
                                            out.write("\t".join([chrom, 'MXE', 'exon', ss1, ss2,\
                                                '.', strand, '.', "ID=" + name + ".B.up;Parent=" + name + ".B"]) + "\n")
                                            out.write("\t".join([chrom, 'MXE', 'exon', ss5, ss6,\
                                                '.', strand, '.', "ID=" + name + ".B.mxe2;Parent=" + name + ".B"]) + "\n")
                                            out.write("\t".join([chrom, 'MXE', 'exon', ss7, ss8,\
                                                '.', strand, '.', "ID=" + name + ".B.dn;Parent=" + name + ".B"]) + "\n")
                                            
                                        else:
                                            upexon = ":".join([chrom, ss2, ss1, strand])
                                            mxe1 = ":".join([chrom, ss4, ss3, strand])
                                            mxe2 = ":".join([chrom, ss6, ss5, strand])
                                            dnexon = ":".join([chrom, ss8, ss7, strand])
                                            name = "@".join([upexon, mxe1, mxe2, dnexon])
                                            out.write("\t".join([chrom, 'MXE', 'gene', ss8, ss1,\
                                                '.', strand, '.', "ID=" + name + ";Name=" + name]) + "\n")
                                            out.write("\t".join([chrom, 'MXE', 'mRNA', ss8, ss1,\
                                                '.', strand, '.', "ID=" + name + ".A;Parent=" + name]) + "\n")
                                            out.write("\t".join([chrom, 'MXE', 'exon', ss2, ss1,\
                                                '.', strand, '.', "ID=" + name + ".A.up;Parent=" + name + ".A"]) + "\n")
                                            out.write("\t".join([chrom, 'MXE', 'exon', ss4, ss3,\
                                                '.', strand, '.', "ID=" + name + ".A.mxe1;Parent=" + name + ".A"]) + "\n")
                                            out.write("\t".join([chrom, 'MXE', 'exon', ss8, ss7,\
                                                '.', strand, '.', "ID=" + name + ".A.dn;Parent=" + name + ".A"]) + "\n")
                                            out.write("\t".join([chrom, 'MXE', 'mRNA', ss8, ss1,\
                                                '.', strand, '.', "ID=" + name + ".B;Parent=" + name]) + "\n")
                                            out.write("\t".join([chrom, 'MXE', 'exon', ss2, ss1,\
                                                '.', strand, '.', "ID=" + name + ".B.up;Parent=" + name + ".B"]) + "\n")
                                            out.write("\t".join([chrom, 'MXE', 'exon', ss6, ss5,\
                                                '.', strand, '.', "ID=" + name + ".B.mxe2;Parent=" + name + ".B"]) + "\n")
                                            out.write("\t".join([chrom, 'MXE', 'exon', ss8, ss7,\
                                                '.', strand, '.', "ID=" + name + ".B.dn;Parent=" + name + ".B"]) + "\n")
    out.close()



# TODO: need more test on RI
def RI(DtoA_F, AtoD_F, DtoA_R, AtoD_R, gff3_f,
       multi_iso=False):
    """
    Define retained introns.

    Arguments are:
    splice site dictionaries, followed by output file and a keyword to 
    denote method of selecting flanking exon coordinates:
      - shortest
      - longest
      - commonshortest
      - commonlongest 
    """
    # return

    print "Generating retained introns (RI)"
    if os.path.isfile(gff3_f):
        print("  - Found file, skipping...")
        return
    out = open(gff3_f, 'w')

    # print "ATOD_F: ", AtoD_F

    for acceptor in AtoD_F:                                 # iterate through acceptors
        chrom, acceptorcoord, strand = acceptor.split(":")
        donors = [x for x in list(AtoD_F[acceptor]) if x in DtoA_R]
                                                           # get the 5'ss in same exon for these acceptors

        if len(donors) > 1:
            for donor in donors:                           # iterate through these 5'ss
                rilist = []
                riAcceptors = [x for x in list(DtoA_R[donor]) if x in AtoD_R]
                # print "riAcceptors: ", riAcceptors
                                                           # get the upstream 3'ss for this 5'ss
                for riAcceptor in riAcceptors:             # iterate through these 3'ss and get their 5'ss
                    riDonors = [x for x in list(AtoD_R[riAcceptor]) if x in DtoA_R]
                    for riDonor in riDonors:               # if the 3'ss upstream of these 5'ss is the upstream
                        upAcceptors = list(DtoA_R[riDonor])# acceptor, this is a retained intron 
                        if acceptor in upAcceptors:
                            rilist.append([riDonor, riAcceptor])
                if len(rilist) > 0:
                    ss1 = acceptor.split(":")[1]
                    donorlist = list(set([x[0].split(":")[1] for x in rilist]))
                    acceptorlist = list(set([x[1].split(":")[1] for x in rilist]))
                    ss4 = donor.split(":")[1]

                    if multi_iso:
                        if strand == '+':
                            upexon = ":".join([chrom, ss1, "|".join(donorlist), strand])
                            dnexon = ":".join([chrom, "|".join(acceptorlist), ss4, strand])
                            name = "@".join([upexon, dnexon])
                            out.write("\t".join([chrom, 'RI', 'gene', ss1, ss4,\
                                '.', strand, '.', "ID=" + name + ";Name=" + name]) + "\n")
                            out.write("\t".join([chrom, 'RI', 'mRNA', ss1, ss4,\
                                '.', strand, '.', "ID=" + name + ".A;Parent=" + name]) + "\n")
                            for i in range(len(rilist)):
                                out.write("\t".join([chrom, 'RI', 'mRNA', ss1, ss4,\
                                    '.', strand, '.', "ID=" + name + "." + LETTERS[i + 1] + ";Parent=" + name]) + "\n")
                            out.write("\t".join([chrom, 'RI', 'exon', ss1, ss4,\
                                '.', strand, '.', "ID=" + name + ".ri;Parent=" + name + ".A"]) + "\n")
                            for i in range(len(rilist)):
                                out.write("\t".join([chrom, 'RI', 'exon', ss1, donorlist[i],\
                                    '.', strand, '.', "ID=" + name + ".up;Parent=" + name + "." + LETTERS[i + 1]]) + "\n")
                                out.write("\t".join([chrom, 'RI', 'exon', acceptorlist[i], ss4,\
                                    '.', strand, '.', "ID=" + name + ".dn;Parent=" + name + "." + LETTERS[i + 1]]) + "\n")

                        else:
                            upexon = ":".join([chrom, ss1, "|".join(donorlist), strand])
                            dnexon = ":".join([chrom, "|".join(acceptorlist), ss4, strand])
                            name = "@".join([upexon, dnexon])
                            out.write("\t".join([chrom, 'RI', 'gene', ss4, ss1,\
                                '.', strand, '.', "ID=" + name + ";Name=" + name]) + "\n")
                            out.write("\t".join([chrom, 'RI', 'mRNA', ss4, ss1,\
                                '.', strand, '.', "ID=" + name + ".A;Parent=" + name]) + "\n")
                            for i in range(len(rilist)):
                                out.write("\t".join([chrom, 'RI', 'mRNA', ss4, ss1,\
                                    '.', strand, '.', "ID=" + name + "." + LETTERS[i + 1] + ";Parent=" + name]) + "\n")
                            out.write("\t".join([chrom, 'RI', 'exon', ss4, ss1,\
                                '.', strand, '.', "ID=" + name + ".up;Parent=" + name + ".A"]) + "\n")
                            for i in range(len(rilist)):
                                out.write("\t".join([chrom, 'RI', 'exon', donorlist[i], ss1,\
                                    '.', strand, '.', "ID=" + name + ".up;Parent=" + name + "." + LETTERS[i + 1]]) + "\n")
                                out.write("\t".join([chrom, 'RI', 'exon', ss4, acceptorlist[i],\
                                    '.', strand, '.', "ID=" + name + ".dn;Parent=" + name + "." + LETTERS[i + 1]]) + "\n")

                    # 2-isoform case
                    else:
                        for iso1 in range(len(donorlist)):
                            for iso2 in range(len(acceptorlist)):
                                if strand == '+':
                                    upexon = ":".join([chrom, ss1, donorlist[iso1], strand])
                                    dnexon = ":".join([chrom, acceptorlist[iso2], ss4, strand])
                                    name = "@".join([upexon, dnexon])

                                    out.write("\t".join([chrom, 'RI', 'gene', ss1, ss4,\
                                        '.', strand, '.', "ID=" + name + ";Name=" + name]) + "\n")
                                    out.write("\t".join([chrom, 'RI', 'mRNA', ss1, ss4,\
                                        '.', strand, '.', "ID=" + name + ".A;Parent=" + name]) + "\n")
                                    out.write("\t".join([chrom, 'RI', 'exon', ss1, ss4,\
                                        '.', strand, '.', "ID=" + name + ".ri;Parent=" + name + ".A"]) + "\n")
                                    out.write("\t".join([chrom, 'RI', 'mRNA', ss1, ss4,\
                                        '.', strand, '.', "ID=" + name + ".B;Parent=" + name]) + "\n")
                                    out.write("\t".join([chrom, 'RI', 'exon', ss1, donorlist[iso1],\
                                        '.', strand, '.', "ID=" + name + ".up;Parent=" + name + ".B"]) + "\n")
                                    out.write("\t".join([chrom, 'RI', 'exon', acceptorlist[iso2], ss4,\
                                        '.', strand, '.', "ID=" + name + ".dn;Parent=" + name + ".B"]) + "\n")

                                else:
                                    upexon = ":".join([chrom, ss1, donorlist[iso1], strand])
                                    dnexon = ":".join([chrom, acceptorlist[iso2], ss4, strand])
                                    name = "@".join([upexon, dnexon])
                                    
                                    out.write("\t".join([chrom, 'RI', 'gene', ss4, ss1,\
                                        '.', strand, '.', "ID=" + name + ";Name=" + name]) + "\n")
                                    out.write("\t".join([chrom, 'RI', 'mRNA', ss4, ss1,\
                                        '.', strand, '.', "ID=" + name + ".A;Parent=" + name]) + "\n")
                                    out.write("\t".join([chrom, 'RI', 'exon', ss4, ss1,\
                                        '.', strand, '.', "ID=" + name + ".ri;Parent=" + name + ".A"]) + "\n")
                                    out.write("\t".join([chrom, 'RI', 'mRNA', ss4, ss1,\
                                        '.', strand, '.', "ID=" + name + ".B;Parent=" + name]) + "\n")
                                    out.write("\t".join([chrom, 'RI', 'exon', donorlist[iso1], ss1,\
                                        '.', strand, '.', "ID=" + name + ".up;Parent=" + name + ".B"]) + "\n")
                                    out.write("\t".join([chrom, 'RI', 'exon', ss4, acceptorlist[iso2],\
                                        '.', strand, '.', "ID=" + name + ".dn;Parent=" + name + ".B"]) + "\n")
                                # Record seen RI
                                #seen_RIs[name] = True
    out.close()

