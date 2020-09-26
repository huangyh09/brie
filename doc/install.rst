============
Installation
============

Environment setting
===================
We recommend to create a separated `conda environment`_ for running BRIE2, which
heavily depends on TensorFlow and TensorFlow-probability.

.. code-block:: bash
  
  conda create -n TFProb python=3.7

replace ``-n TFProb`` with ``-p ANY_PATH/TFProb`` to specify the path for conda 
environment. Then activate the environment by ``conda activate TFProb`` or the 
full path, before install more packages.

.. _conda environment: https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html


Easy install
============

BRIE2 is available through `pypi`_. To install, type the following command 
line, and add ``-U`` for upgrading:

.. code-block:: bash

  pip install -U brie

Alternatively, you can install from this GitHub repository for latest (often 
development) version by following command line

.. code-block:: bash

  pip install -U git+https://github.com/huangyh09/brie

In either case, if you don't have write permission for your current Python 
environment, add ``--user``, but check the previous section on create your own
conda environment.

.. _pypi: https://pypi.org/project/brie


GPU usage
=========
With TensorFlow backend, BRIE2 can benefit from using GPUs. Here is one way to 
set up GPU configurations with NVIDIA GPU on Ubuntu:

.. code-block:: bash

  pip install -U tensorflow-gpu
  conda install -c anaconda cupti 
  conda install -c anaconda cudnn
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda/extras/CUPTI/lib64

For more information on GPU configuration, please refer to the 
`Tensorflow documentation`_, or `anaconda GPU`_.

.. _Tensorflow documentation: https://www.tensorflow.org/guide/gpu
.. _anaconda GPU: https://docs.anaconda.com/anaconda/user-guide/tasks/gpu-packages/


.. note::
   At the moment, TensorFlow calls all available GPUs, which is not nessary. 
   You can specify the card you want to use by add the following variable before
   you command line ``CUDA_VISIBLE_DEVICES=3 brie-quant -i my_count.h5ad``
   


Test
====

In order to test the installation, you could type ``brie-quant``. If successful,
you will see the following output.

.. code-block:: html

  Welcome to brie-quant in BRIE v2.0.2!

  use -h or --help for help on argument.

If installation is sucessful, but can't run it, then check whether the directory 
which contains the executable binary file is added to PATH environment. 

.. code-block:: html

  brie-quant: command not found

Usually the directory is ``~/.local/bin`` if you don't use Anaconda. You could add 
the path into PATH environment variable, by write the following line into ``.profile`` 
or ``.bashrc`` file.

.. code-block:: html
  
  export PATH="~/.local/bin:$PATH"


If you have any issue, please report it to the issue on `brie issues`_.

.. _brie issues: https://github.com/huangyh09/brie/issues

