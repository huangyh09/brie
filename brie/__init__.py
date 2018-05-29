
__version__ = "0.1.4"

# import pyximport; pyximport.install()
# from .utils.fasta import *


from .utils.bias_utils import FastaFile, BiasFile
from .utils.tran_utils import TranUnits, TranSplice
from .utils.sam_utils import load_samfile, fetch_reads
from .utils.gtf_utils import Gene, Transcript, loadgene, savegene
