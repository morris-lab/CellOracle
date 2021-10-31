.. _install:

Installation
============

CellOracle uses several python libraries and R libraries. Please follow this guide below to install CellOracle and its dependent software.

.. _require:

Docker image
------------

- Pre-built docker image is available through `Docker Hub <https://hub.docker.com/repository/docker/kenjikamimoto126/celloracle_ubuntu>`_ .

::

    docker push kenjikamimoto126/celloracle_ubuntu:latest


- This docker image was built based on Ubuntu 20.04.
- Python dependent packages and celloracle are installed under an anaconda environment, celloracle_env. This environment will be activated automatically when you log in.
- R dependent libraries for network analysis are installed. Also, Seurat V3, Monocle3, and Cicero are installed.
- After logging in, the user switches from the root user to the following user. Username: user. Password: pass.


Install CellOracle
------------------

System Requirements
^^^^^^^^^^^^^^^^^^^

- Operating system: macOS or Linux are highly recommended. CellOraclewas developed and tested in Linux and macOS.
- We found that the celloracle calculation may be EXTREMELY SLOW under an environment of Windows Subsystem for Linux (WSL). We do not recommend using WSL.
- While you can install CellOracle in Windows OS, please do so at your own risk and responsibility. We DO NOT provide any support for the use in the Windows OS.

- Memory: 16 G byte or more.  Memory usage also depends on your data. Especially in silico perturbation requires large amount of memory.
- CPU: Core i5 or better processor. GRN inference supports multicore calculation. Higher number of CPU cores enables fast calculation.


Python Requirements
^^^^^^^^^^^^^^^^^^^

- CellOracle was developed with python 3.6. We do not support python 2.7x or python <=3.5.
- Please install all dependent libraries before installing CellOracle according to the instructions below.
- CellOracle is still a beta version and it is not available through PyPI or anaconda distribution yet. Please install CellOracle from our GitHub repository according to the instruction below.

CellOracle installation using conda and pip
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**1. Make a conda environment**
  We recommend installing CellOracle in an independent conda environment to avoid dependent software conflicts.
  Please make a new python environment for celloracle and install dependent libraries in it.

  ::

      conda create -n celloracle_env python=3.6
      conda activate celloracle_env

  Installation of some libraries requires non-default anaconda channels. Please add the channels below. Instead, you can explicitly enter the channel when you install a library.

  ::

      conda config --add channels defaults
      conda config --add channels bioconda
      conda config --add channels conda-forge

**2. Install dependencies using conda**

  Run the following command to install some dependencies prior to celloracle installation.

  ::

      conda install pybedtools numba

**3. Install CellOracle and other dependencies**

  ::

      pip install git+https://github.com/morris-lab/CellOracle.git


You may have an error in the installation process of CellOracle dependent libraries.
If you have an error, please look at the troubleshooting page.


.. toctree::
   :maxdepth: 1

   python_step_by_step_installation


R requirements
^^^^^^^^^^^^^^

CellOracle uses R libraries to calculate network graph score.
Please install `R <https://www.r-project.org>`_ (>=3.5) and R libraries below.

.. note::
   These R libraries are needed for network analysis.
   CellOracle gene perturbation simulation does not require the R libraries.
   **You can skip R library installation if you do not perform network analysis.**

.. code-block:: r

   install.packages("igraph")
   install.packages("rnetcarto")
   install.packages("linkcomm")

If you have an error when installing these R libraries above, please look at the troubleshooting tips below.


.. toctree::
   :maxdepth: 1

   R_step_by_step_installation



Check installation
^^^^^^^^^^^^^^^^^^

.. toctree::
   check_installation



Optional R libraries for input data preparation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We provide many working examples for input data preparation.
These R packages below are not in the part of the CellOracle library itself and not necessary.
However you can use them in the input data preparation step if you want.
Please install them on demand.
If you want to try CellOracle main tutorials, :doc:`networkanalysis` and :doc:`simulation`, you DO NOT need to install the libraries below.

- `Seurat <https://satijalab.org/seurat/install.html>`_
- `Cicero <https://cole-trapnell-lab.github.io/cicero-release/docs/#installing-cicero>`_
