"""
Brie - Bayesian regression for isoform estimate
See: http://brie.sourceforge.net
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path
import brie

here = path.abspath(path.dirname(__file__))

# Get the long description from the relevant file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()
    
reqs = ['numpy', 'pysam', 'h5py']

setup(
    name='brie',

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version=brie.__version__,

    description='BRIE: Bayesian regression for isoform estimate',
    long_description=long_description,

    # The project's main homepage.
    url='http://brie-rna.sourceforge.net',

    # Author details
    author='Yuanhua Huang',
    author_email='Y.Huang@ed.ac.uk',

    # Choose your license
    license='MIT',

    # What does your project relate to?
    keywords=['splicing isoform estimate', 'Bayesian regression', 
              'single cell RNA-seq', 'Markov chain Monte Carlo'],

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=find_packages(),

    entry_points={
          'console_scripts': [
              'brie = brie.brie:main',
              'brie-diff = brie.brie_diff:main',
              'brie-factor = brie.brie_factor:main',
              'brie-event = brie.events.event_maker:main',
              'brie-event-filter = brie.events.event_filter:main',
              ],
          }, 

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    
    install_requires=reqs,

    py_modules = ['brie']

    # buid the distribution: python setup.py sdist
    # upload to pypi: twine upload dist/...

)