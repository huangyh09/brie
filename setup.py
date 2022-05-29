"""
Brie - Bayesian regression for isoform estimate
See: https://brie.readthedocs.io
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from pathlib import Path

# Set __version__ for the project.
exec(open("./brie/version.py").read())

setup(
    name='brie',

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version=__version__,

    description='BRIE: Bayesian regression for isoform estimate',
    long_description=Path("README.rst").read_text("utf-8"),

    # The project's main homepage.
    url='https://brie.readthedocs.io',

    # Author details
    author='Yuanhua Huang',
    author_email='yuanhua@hku.hk',

    # Choose your license
    license='Apache-2.0',

    # What does your project relate to?
    keywords=['RNA splicing', 'Bayesian regression',
              'single cell RNA-seq', 'variantional inference'],

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=find_packages(),

    entry_points={
        'console_scripts': [
            # 'brie = brie:_cli.cli',
            'brie = brie.bin.brie_main:main',
            'brie-quant = brie.bin.quant:main',
            'brie-count = brie.bin.count:main',
            'brie1 = brie.version1.brie:main',
            'brie1-diff = brie.version1.brie_diff:main',
        ],
    },
    
    python_requires='>=3.5',

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html

    # install_requires=reqs,
    install_requires=[
        l.strip() for l in Path('requirements.txt').read_text('utf-8').splitlines()
    ],

    extras_require={
        'docs': [
            # 'sphinx >= 1.4',
            'sphinx_bootstrap_theme']},

    py_modules=['brie']

    # buid the distribution: python setup.py sdist
    # upload to pypi: twine upload dist/...

)
