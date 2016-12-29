##
## Make a GFF annotation of alternative splicing
## events from a set of UCSC tables.
##
import os
import sys
import time
from optparse import OptionParser, OptionGroup

import parseTables
from defineEvents import SE, MXE, RI, A3SS, A5SS

# import gffutils
# import gffutils.helpers as helpers


def prepareSplicegraph(anno_file, ftype):
    """
    Prepare splicegraph for use in defining events.
    anno_file could be multiple files with delimiter of comma;
    ftype should have the same numbers as anno_file, or just one.
    Returns the splice site dictionaries.
    """
    DtoA_F = {}
    AtoD_F = {}
    DtoA_R = {}
    AtoD_R = {}

    ftype = ftype.split(",")
    anno_fnames = anno_file.split(",")
    if len(ftype) == 1 and len(anno_fnames) > 1:
        ftype = ftype * len(anno_fnames)

    for i in range(len(anno_fnames)): 
        print "Reading table", anno_fnames[i]
        DtoA_F, AtoD_F, DtoA_R, AtoD_R = \
            parseTables.populateSplicegraph(anno_fnames[i], ftype[i],
                                            DtoA_F, AtoD_F, DtoA_R, AtoD_R)
    DtoA_F, AtoD_F, DtoA_R, AtoD_R = \
        parseTables.cleanSplicegraph(DtoA_F, AtoD_F, DtoA_R, AtoD_R)

    return DtoA_F, AtoD_F, DtoA_R, AtoD_R



def defineAllSplicing(anno_file, ftype, gff3dir,
                      flanking='commonshortest',
                      multi_iso=False,
                      genome_label=None,
                      sanitize=False,
                      event_types=["SE", "RI", "MXE", "A3SS", "A5SS"]):
    """
    A wrapper to define all splicing events: SE, RI, MXE, A3SS, A5SS
    RI does not use the "flanking criteria".
    """
    if isinstance(multi_iso, str):
        multi_iso = eval(multi_iso)

    DtoA_F, AtoD_F, DtoA_R, AtoD_R = prepareSplicegraph(anno_file, ftype)

    # Encode the flanking exons rule in output directory
    # gff3dir = os.path.join(gff3dir, flanking)
    if os.path.isfile(gff3dir):
        print "Error: %s is a file!" %(gff3dir)
        sys.exit(1)
    try:
        os.makedirs(gff3dir)
    except OSError:
        pass

    if genome_label is not None:
        genome_label = ".%s" %(genome_label)
    else:
        genome_label = ""

    fname_all = []

    if "SE" in event_types:
        out_fname = os.path.join(gff3dir, "SE%s.gff3" %(genome_label))
        fname_all.append(out_fname)
        SE(DtoA_F, AtoD_F, DtoA_R, AtoD_R, out_fname, flanking)

    if "RI" in event_types:
        out_fname = os.path.join(gff3dir, "RI%s.gff3" %(genome_label))
        fname_all.append(out_fname)
        RI(DtoA_F, AtoD_F, DtoA_R, AtoD_R, out_fname, multi_iso)

    if "MXE" in event_types:
        out_fname = os.path.join(gff3dir, "MXE%s.gff3" %(genome_label))
        fname_all.append(out_fname)
        MXE(DtoA_F, AtoD_F, DtoA_R, AtoD_R, out_fname, flanking)

    if "A3SS" in event_types:
        out_fname = os.path.join(gff3dir, "A3SS%s.gff3" %(genome_label))
        fname_all.append(out_fname)
        A3SS(DtoA_F, AtoD_F, DtoA_R, AtoD_R, out_fname, flanking, multi_iso)

    if "A5SS" in event_types:
        out_fname = os.path.join(gff3dir, "A5SS%s.gff3" %(genome_label))
        fname_all.append(out_fname)
        A5SS(DtoA_F, AtoD_F, DtoA_R, AtoD_R, out_fname, flanking, multi_iso)
        
    # If asked, sanitize the annotation in place
    if sanitize:
        for _fname in fname_all:
            print "Sanitizing %s" %(_fname)
            helpers.sanitize_gff_file(_fname, in_place=True)


def main():
    """
    Make GFF annotation. Takes GFF tables directory
    and an output directory.
    """
    # load parse command line options
    parser = OptionParser()
    parser.add_option("--anno_file", "-a", default=None,
                        help="The annotation files used in making the "
                        "annotation. You could input multiple files; use comma"
                        "',' as delimiter.")
    parser.add_option("--anno_type", default="gtf",
                        help="The type of each annotation file. Use one for "
                        "all files or set for each file. Use comma ',' as "
                        "delimiter. You could choose 'ucsc', 'gtf', 'gff3'. "
                        "[default: %default]")
    parser.add_option("--output_dir", "-o", help="Output directory.")
    parser.add_option("--flanking-rule", default="commonshortest",
                        help="Rule to use when defining exon trios. "
                        "E.g. \'commonshortest\' to use the most common "
                        "and shortest regions are flanking exons to an "
                        "alternative trio. [default: %default]")
    parser.add_option("--multi-iso", action="store_true",
                        default=False, help="If passed, generates "
                        "multi-isoform annotations. Off by default.")
    parser.add_option("--genome-label", help="If given, used as label for "
                        "genome in output files.")
    parser.add_option("--sanitize", default=False, action="store_true",
                        help="If passed, sanitize the annotation. "
                        "Off by default.")
    (options, args) = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        print("Welcome to Splicing-event-maker!\n")
        print("use -h or --help for help on argument.")
        sys.exit(1)
    
    # process
    if options.anno_file is None:
        print("Error: need --anno_file for annotation.")
        sys.exit(1)

    if options.output_dir is None:
        output_dir = os.path.join(os.path.split(options.anno_file)[0], 
            "AS_events")
    else:
        output_dir = options.output_dir

    print "Making GFF alternative events annotation..."
    print "  - Input annotation files: %s" %(options.anno_file)
    print "  - Output dir: %s" %(output_dir)

    t1 = time.time()
    defineAllSplicing(options.anno_file, options.anno_type, output_dir,
                         flanking=options.flanking_rule,
                         multi_iso=options.multi_iso,
                         genome_label=options.genome_label,
                         sanitize=options.sanitize)
    t2 = time.time()
    print "Took %.2f minutes to make the annotation." \
          %((t2 - t1)/60.)
          

if __name__ == "__main__":
    main()
