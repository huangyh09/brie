============
Installation
============

Environment setting (optional)
==============================
For using BRIE2, we recommend creating a separated `conda environment`_ for 
easier and cleaner management of dependent packages, particularly TensorFlow and 
TensorFlow-probability, which are under active development.

The following command line in terminal (Linux or MacOS) will create a conda 
environment with name ``TFProb``, probably in the path ``~/.conda/envs/TFProb``:

.. code-block:: bash
  
  conda create -n TFProb python=3.7

Alternatively, you can create the environment in another path by replacing 
``-n TFProb`` with ``-p ANY_PATH/TFProb`` to specify the path for conda 
environment. Then you can check your environment by ``conda env list`` and 
activate the environment by ``conda activate TFProb`` or the full path, before 
start installing brie and other packages.

.. _conda environment: https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html


Easy install
============

BRIE2 is available on `PYPI`_. To install, type the following command 
line, and add ``-U`` for upgrading:

.. code-block:: bash

  pip install -U brie

Alternatively, you can install from this GitHub repository for the latest (often 
development) version by the following command line

.. code-block:: bash

  pip install -U git+https://github.com/huangyh09/brie

In either case, if you don't have write permission for your current Python 
environment, add ``--user``, but check the previous section above for creating 
your own conda environment.

.. _PYPI: https://pypi.org/project/brie


GPU usage
=========
One of the key benefits of using TensorFlow backend is its direct support of 
GPU for substantial speedups (generally >10x). Here is one way to set up GPU 
configurations with NVIDIA GPU on Ubuntu:

.. code-block:: bash

  pip install -U tensorflow-gpu
  conda install -c anaconda cudatoolkit

Make sure that you have compatible versions between tensorflow and NVIDIA CUDA. 
You can check TF's test `here <https://www.tensorflow.org/install/source#gpu>`_.
For more information on GPU configuration, please refer to the 
`Tensorflow documentation`_ or `anaconda GPU`_.

.. _Tensorflow documentation: https://www.tensorflow.org/guide/gpu
.. _anaconda GPU: https://docs.anaconda.com/anaconda/user-guide/tasks/gpu-packages/

Once successfully configured, you will see log info like the following when 
using ``brie-quant``:

.. code-block:: html

  $ brie-quant
    I tensorflow/stream_executor/platform/default/dso_loader.cc:53] 
    Successfully opened dynamic library libcudart.so.11.0
    Welcome to brie-quant in BRIE v2.0.5!


.. note::
   At the moment, TensorFlow calls all available GPUs, which is not necessary. 
   You can specify the card (e.g., card 3) by adding the below variable before 
   your command line ``CUDA_VISIBLE_DEVICES=3 brie-quant -i my_count.h5ad``
   


Test
====

In order to test the installation, you could type ``brie-quant``. If successful,
you will see the following output.

.. code-block:: html

  Welcome to brie-quant in BRIE v2.0.2!

  use -h or --help for help on argument.

If you install BRIE successfully, but can't run it and see the error below, 
then check whether its directory is added to ``PATH`` environment. 

.. code-block:: html

  brie-quant: command not found

Usually, the directory is ``~/.local/bin`` if you don't use Anaconda. You could add 
the path into ``PATH`` environment variable, by writing the following line into 
``.profile`` or ``.bashrc`` file.

.. code-block:: html
  
  export PATH="~/.local/bin:$PATH"


If you have any issues, please report them to the issue on `brie issues`_.

.. _brie issues: https://github.com/huangyh09/brie/issues

