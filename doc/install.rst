============
Installation
============

Environment setting (optional)
==============================
For using BRIE2, we recommend creating a separated `conda environment`_ for 
easier and cleaner management of dependent packages, particularly TensorFlow and 
TensorFlow-probability, which are under active development.

The following command line in terminal (Linux or MacOS) will create a conda 
environment with name ``TFProb``, probably in the path ``~/.conda/envs/TFProb``
(tested on Python 3.7 to 3.11):

.. code-block:: bash
  
  conda create -n TFProb python=3.11

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
GPU for substantial speedups (generally >10x). In general, the above 
installation should directly support GPU use as default by using the newest 
tensorflow and tensorflow-probability. 

To check if the installation is compatible with GPUs, you can print out the 
detectable GPU cards, as below (it gives `[ ]` if failing to setup properly):

.. code-block:: html

  $ python3 -c "import tensorflow as tf; print(tf.config.list_physical_devices('GPU'))"
    [PhysicalDevice(name='/physical_device:GPU:0', device_type='GPU'), ...]

However, occasionally, the newest TensorFlow may not be stable or compatible
widely, for example the `tensorflow[and-cuda]==2.16.1` is not compatible with GPUs,
see discussion `here <https://github.com/tensorflow/tensorflow/issues/63362#issuecomment-2053849484>`_.

Here is one way to use a lower version for GPU configurations with NVIDIA GPU 
on Ubuntu (tested on 21/04/2024):

.. code-block:: bash

  pip install tensorflow-probability==0.23.0
  pip install tensorflow[and-cuda]==2.15.1

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

