"""
BRIE package
============
Bayesian Regression for Isoform Estimate - a python toolbox for splicing 
analysis in single cells.
The official documentation together with examples and tutorials can be found
at https://brie.readthedocs.io/.
"""

from .version import __version__
from .models.model_TFProb import BRIE2
from .utils.io_utils import load_brie_count
from .utils.bias_utils import FastaFile, BiasFile
from .utils.tran_utils import TranUnits, TranSplice
from .utils.sam_utils import load_samfile, fetch_reads
from .utils.gtf_utils import Gene, Transcript, loadgene, savegene


__all__ = [
    "__version__",
    "utils",
    "models"
]
