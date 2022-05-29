"""
BRIE package
============
Bayesian Regression for Isoform Estimate - a python toolbox for splicing 
analysis in single cells.
The official documentation together with examples and tutorials can be found
at https://brie.readthedocs.io/.
"""

# from ._cli import cli
from .version import __version__

# direct classes or functions
# from .models.model_TFProb import BRIE2
from .utils.io_utils import read_brieMM, read_h5ad, read_gff, read_npz
from .utils import io_utils as io
from .utils.base_utils import match

# set simplified alias
from . import plot as pl
# from .models import tools as tl
from .utils import preprocessing as pp


# __all__ = [
#     "__version__",
#     "utils",
#     "models"
# ]
