.. _install:

Installation
============

``celloracle`` uses several python libraries and R library. Please follow this guide below to install the dependent software of celloracle.

.. _require:

Python Requirements
-------------------

- ``celloracle`` was developed with python 3.6. We do not support python 2.7x or python <=3.5.
- ``celloracle`` was developed in Linux and macOS. We do not guarantee that ``celloracle`` works in Windows OS.
- If you want to use ``celloracle`` in Windows OS, we recommend using `Windows Subsystem for Linux <https://docs.microsoft.com/en-us/windows/wsl/install-win10>`_ to built environment. In that case, please select Ubuntu18.04 LTS.
- We highly recommend using `anaconda <https://www.continuum.io/downloads>`_ to setup python environment.
- Please install all dependent libraries before installing ``celloracle`` according to the instructions below.
- ``celloracle``  is still beta version and it is not available through PyPI or anaconda distribution yet. Please install ``celloracle`` from GitHub repository according to the instruction below.


0. (Optional) Make a new environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This step is optional. Please make a new python environment for celloracle and install dependent libraries in it if you get some software conflicts.

::

    conda create -n celloracle_env python=3.6
    conda activate celloracle_env



1. Add conda channels
^^^^^^^^^^^^^^^^^^^^^
Installation of some libraries below requires non-default anaconda channels. Please add the channels below. Instead, you can explicitly enter the channel when you install a library.

::

    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge


2. Install `velocyto <http://velocyto.org/velocyto.py/install/index.html>`_
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Please install velocyto with the following commands or `the author's instruction <http://velocyto.org/velocyto.py/install/index.html>`_ .
On Mac OS, you may have a compile error during velocyto installation. I recommend you installing `Xcode <https://developer.apple.com/xcode/>`_ in that case.


::

    conda install numpy scipy cython numba matplotlib scikit-learn h5py click

Then

::

    pip install velocyto

3. Install `scanpy <https://scanpy.readthedocs.io/en/stable/installation.html>`_
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Please install scanpy with the following commands or `the author's instruction <https://scanpy.readthedocs.io/en/stable/installation.html>`_ .

::

    conda install seaborn scikit-learn statsmodels numba pytables python-igraph louvain

Then

::

    pip install scanpy

4. Install `gimmemotifs <https://gimmemotifs.readthedocs.io/en/master/installation.html>`_
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Please install gimmemotifs with the following commands or `the author's instruction <https://gimmemotifs.readthedocs.io/en/master/installation.html>`_ .


::

    conda install gimmemotifs genomepy=0.5.5


5. Install other python libraries
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Please install other python libraries below with the following commands.

::

    conda install goatools pyarrow tqdm joblib jupyter


6. install celloracle from github
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
::

    pip install git+https://github.com/morris-lab/CellOracle.git



R requirements
--------------

``celloracle`` use R library for the network analysis and scATAC-seq analysis.
Please install `R <https://www.r-project.org>`_ (>=3.5) and R libraries below according to the author's instruction.

`Seurat <https://satijalab.org/seurat/install.html>`_
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Please install ``Seurat`` with the following r-script or `the author's instruction <https://satijalab.org/seurat/install.html>`_ .
``celloracle`` is compatible with both Seurat V2 and V3.
If you use only ``scanpy`` for the scRNA-seq preprocessing and do not use ``Seurat`` , you can skip installation of ``Seurat``.

in R console,

.. code-block:: r

   install.packages('Seurat')

`Cicero <https://cole-trapnell-lab.github.io/cicero-release/docs/#installing-cicero>`_
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Please install ``Cicero`` with the following r-script or `the author's instruction <https://cole-trapnell-lab.github.io/cicero-release/docs/#installing-cicero>`_ .
If you have no plan for scATAC-seq analysis and just want to use  ``celloracle`` with a default TF information which was supplied with celloracle, you can skip installation of ``Cicero``.

in R console,

.. code-block:: r

   if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
   BiocManager::install("cicero", version = "3.8")

`igraph <https://igraph.org/r/>`_
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Please install ``igraph`` with the following r-script or `the author's instruction <https://igraph.org/r/>`_ .

in R console,

.. code-block:: r

   install.packages("igraph")


`linkcomm <https://cran.r-project.org/web/packages/linkcomm/index.html>`_
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Please install ``linkcomm`` with the following r-script or `the author's instruction <https://cran.r-project.org/web/packages/linkcomm/index.html>`_ .

in R console,

.. code-block:: r

   install.packages("linkcomm")

`rnetcarto <https://github.com/cran/rnetcarto/blob/master/src/rgraph/README.md>`_
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
``rnetcarto`` installation has to be done with several steps. Please install rnetcarto with `the author's instruction <https://github.com/cran/rnetcarto/blob/master/src/rgraph/README.md>`_ .
You need to install `the GNU Scientific Libraries <https://www.gnu.org/software/gsl/>`_ before installing rnetcarto. Detailed instruction can be found `here <https://github.com/cran/rnetcarto/blob/master/src/rgraph/README.md>`_ .


Check installation
^^^^^^^^^^^^^^^^^^
These R libraries above are necessary for the network analysis in celloracle. You can check installation using celloracle's function.

in python console,

.. code-block:: Python

   import celloracle as co
   co.network_analysis.test_R_libraries_installation()

Please make sure that all R libraries are installed. The following message will be shown when all R libraries are appropriately installed.

| checking R library installation: igraph -> OK
| checking R library installation: linkcomm -> OK
| checking R library installation: rnetcarto -> OK
