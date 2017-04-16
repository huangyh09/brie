============
Installation
============

Easy install
============

* Required packages in python: `numpy`, `matplotlib`, `pysam`, `h5py`

  * we suggest using Anaconda_ distribute, which includes most packages (except 
    `pysam` here), and provides you a user specific environment, i.e., all your 
    python packages go to a single folder. Thus, you don't need the root to 
    install packages.

  * you could install `pysam` by pypi in terminal, or download_ and install it 
    in the same way as `BRIE`.

  .. _Anaconda: http://continuum.io/downloads
  .. _download: https://github.com/pysam-developers/pysam

* You can install `BRIE` simply via pypi in terminal (suggested), or upgrade 
  by add ``--upgrade`` as follows:

  ::

    pip install brie

    pip install --upgrade --no-deps brie


Source code
===========

* Alternatively, you also could download the source code via GitHub (latest 
  version, suggested) or Sourceforge (any version) and run python setup in 
  terminal:

  * GitHub: https://github.com/huangyh09/brie

  * Sourceforge: http://sourceforge.net/projects/brie-rna/

  ::

    python setup.py install

* In any case, if had the permission error for installation as you are not 
  root, add ``--user``.


Test
====

In order to test the installation, you could type ``brie``. If successful, you
will see the following output.

.. code-block:: html

  Welcome to Brie!

  use -h or --help for help on argument.

If installation is sucessful, but can't run it, then check whether the directory 
which contains the executable binary file is added to PATH environment. 

.. code-block:: html

  brie: command not found

Usually the directory is ``~/.local/bin`` if you don't use Anaconda. You could add 
the path into PATH environment variable, by write the following line into ``.profile`` 
or ``.bashrc`` file.

:: 
  
  export PATH="~/.local/bin:$PATH"