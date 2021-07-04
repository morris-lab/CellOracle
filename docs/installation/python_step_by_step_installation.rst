.. _python_step_by_step_installation:

Python dependent library installation troubleshooting
=====================================================

0. (Optional) Make a new conda environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This step is optional, but we recommend installing CellOracle in an independent conda environment to avoid dependent software conflicts.
Please make a new python environment for celloracle and install dependent libraries in it.

::

    conda create -n celloracle_env python=3.6
    conda activate celloracle_env



1. Add conda channels
^^^^^^^^^^^^^^^^^^^^^
Installation of some libraries requires non-default anaconda channels. Please add the channels below. Instead, you can explicitly enter the channel when you install a library.

::

    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge


2. Install `velocyto <http://velocyto.org/velocyto.py/install/index.html>`_
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Please install velocyto with the following commands or `the author's instruction <http://velocyto.org/velocyto.py/install/index.html>`_ .

::

    conda install numpy scipy cython numba matplotlib scikit-learn h5py>=3.1.0 click pysam llvm louvain

Then

::

    pip install velocyto

It was reported that some compile errors might occur during the installation of velocyto on MacOS.
Various errors were reported, and you need to find the best solution depending on your error.
You may find the solution with these links below.

- `Solution 1: Install Xcode <https://developer.apple.com/xcode/>`_. Please try this first.
- `Solution 2: Install macOS_SDK_headers <https://stackoverflow.com/a/53057706/10641716>`_. This solution is needed in addition to Solution-1 if your OS is macOS Mojave.
- `Solution 3 <https://github.com/morris-lab/CellOracle/issues/3>`_. This is the solution reported by a CellOracle user. Thank you very much!
- `Other solutions on Velocyto GitHub issue page <https://github.com/velocyto-team/velocyto.py/issues?q=>`_

3. Install `scanpy <https://scanpy.readthedocs.io/en/stable/installation.html>`_
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Please install scanpy with the following commands or `the author's instruction <https://scanpy.readthedocs.io/en/stable/installation.html>`_ .

::

    conda install scanpy



4. Install other python libraries
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Please install other python libraries below with the following commands.

::

    conda install goatools pyarrow tqdm joblib jupyter gimmemotifs==0.14.4 genomepy==0.8.4


5. install celloracle
^^^^^^^^^^^^^^^^^^^^^
::

    pip install git+https://github.com/morris-lab/CellOracle.git
