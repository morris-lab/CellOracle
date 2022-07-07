.. _python_step_by_step_installation:

Python dependent library installation troubleshooting
=====================================================


`velocyto <http://velocyto.org/velocyto.py/install/index.html>`_
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you failed CellOracle installation because of velocyto installation error, please try to install velocyto using the following commands or please follow `the author's instruction <http://velocyto.org/velocyto.py/install/index.html>`_ .

::

    conda install numpy scipy cython numba matplotlib scikit-learn h5py click pysam llvm louvain

Then

::

    pip install velocyto

Some MacOS users have encountered compiler errors during the installation of velocyto.
Several different errors have been reported.
You may find the solution for each of these in the links below.

- `Solution 1: Install Xcode <https://developer.apple.com/xcode/>`_. Please try this first.
- `Solution 2: Install macOS_SDK_headers <https://stackoverflow.com/a/53057706/10641716>`_. This solution is needed in addition to Solution-1 if your OS is macOS Mojave.
- `Solution 3 <https://github.com/morris-lab/CellOracle/issues/3>`_. This is the solution reported by a CellOracle user. Thank you very much!
- `Other solutions on Velocyto GitHub issue page <https://github.com/velocyto-team/velocyto.py/issues?q=>`_

`scanpy <https://scanpy.readthedocs.io/en/stable/installation.html>`_
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you failed CellOracle installation because of Scanpy installation error, please try to install Scanpy using the following commands or `the author's instruction <https://scanpy.readthedocs.io/en/stable/installation.html>`_ .

::

    conda install scanpy


`pybedtools <https://daler.github.io/pybedtools/>`_
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you have an error related to the pybedtools installation, please install pybedtools using anaconda.

::

    conda install -c bioconda pybedtools

.. note::
   Pybedtools is python wrapper of bedtools. If you get "NotImplementedError" during ATAC-seq peak processing step, please try to install bedtools.

   ::

       sudo apt-get install bedtools



Install gimmemotifs with conda
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you failed installation in the gimmemotifs installation step, please install gimmemotifs solely first.
We recommend installing gimmemotifs using conda prior to celloracle installation.

::

    conda install -c bioconda gimmemotifs


.. note::
   If your jupyter notebook kernel keep dying during TF motif scan, there are mainly two possibilities below.

   1. Your kernel might be killed because of memory shortage. Please make sure you have enough memory.

   2. If you have enough memory, there is most likely a problem with the gimmemotifs installation in your environment.
   Please uninstall gimmemotifs and re-install gimmemotifs using conda.

   ::

       pip uninstall gimmemotifs -y
       conda install -c bioconda gimmemotifs



.. warning::
   We found gimmemotifs might have installation issue with anaconda bioconda channel (at 7/7/2022).
   If you fail to install gimmemotifs with anaconda, we recommend installing it using anaconda's offline mode.
   More information can be found `here <https://github.com/vanheeringen-lab/gimmemotifs/issues/271>`_.


If you encounter an error related to "certifi"
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you encounter the following error, it means the error is caused by versiom mismatch of "certifi".
See `this page  <https://stackoverflow.com/questions/50129762/graphlab-create-2-1-installation-fails-to-uninstall-certifi-a-distutils-insta>`_. for more information.

::

    ERROR: Cannot uninstall 'certifi'. It is a distutils installed project and thus we cannot accurately determine which files belong to it which would lead to only a partial uninstall.

In this case, please add "--ignore-installed certifi " to the installation command.

::

    pip install celloracle --ignore-installed certifi
