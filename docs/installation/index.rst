.. _install:

Installation
============

Please follow the guide below to install CellOracle and its dependent software.

.. _require:

Docker image
------------

- Pre-built docker image is available through `Docker Hub <https://hub.docker.com/repository/docker/kenjikamimoto126/celloracle_ubuntu>`_ .

::

    docker pull kenjikamimoto126/celloracle_ubuntu:0.18.0


- This docker image was built based on Ubuntu 20.04 and python 3.10.
- Python dependent packages and celloracle are installed in an anaconda environment, celloracle_env. This environment will be activated automatically when you log in.

.. toctree::
   :maxdepth: 1

   docker_additional_information



Install CellOracle
------------------

System Requirements
^^^^^^^^^^^^^^^^^^^

- Operating system: macOS or Linux are highly recommended. CellOracle was developed and tested in Linux and macOS.
- CellOracle cannot be installed in Windows OS because some of the dependencies cannot be built on Windows OS.
- It is possible to install CellOracle on a Windows Subsystem for Linux (WSL), but we found that the CellOracle calculations may be extremely slow on a WSL. We do not recommend using WSL.
- While you can install CellOracle using the WSL, please do so at your own risk and responsibility. We DO NOT provide any support for the use with the Windows OS.

- Memory: 16 G byte or more.  Memory usage also depends on the size of your scRNA-seq dataset. Please note that in silico perturbation, may require large amount of memory.
- CPU: Core i5 or better processor. CellOracle's GRN inference function supports multicore calculation. Utilizing more CPU cores enables faster calculations.


Python Requirements
^^^^^^^^^^^^^^^^^^^

- CellOracle was developed using python 3.8. It is also tested with python 3.10. We do not support python 2.7x or python <=3.5.
- CellOracle is available through `PyPI <https://pypi.org/project/celloracle/>`_. You can also download and install celloracle from our `CellOracle  GitHub repository <https://github.com/morris-lab/CellOracle>`_ .

CellOracle installation using conda and pip
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  We recommend installing CellOracle in an independent conda environment to avoid dependent software conflicts.
  Please make a new python environment for celloracle and install dependent libraries in it.

  ::

      conda create -n celloracle_env python=3.10
      conda activate celloracle_env


  You can install CellOracle and all dependencies with the following command.

  ::

      pip install celloracle


If you get an error, please look at the troubleshooting page below.
We also recommend using our pre-built docker image if you have program conflicts between CellOracle's dependent libraries and your pre-existing software.

.. toctree::
   :maxdepth: 1

   python_step_by_step_installation


.. warning::
  Previous version of CellOracle required some R packages. But after celloracle 0.10.0, the R functions were replaced with python codes, and **celloracle does NOT require R packages after version 0.10.0.**


Check installation
^^^^^^^^^^^^^^^^^^

You can check the installed library version as follows. Please make sure that all dependencies are appropriately installed.

In python console,

.. code-block:: Python

   import celloracle as co
   co.check_python_requirements()



Other optional softwares for input data preparation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We provide many working examples for input data preparation.
These softwares below are not in the part of the CellOracle library itself and not necessary.
However, some users may like to use these R libraries to preprocess the input data before beginning a CellOracle analysis.
If this applies to you, please install these R libraries yourself and follow the related tutorials.
If you just want to try CellOracle main tutorials, :doc:`networkanalysis` and :doc:`simulation`, you DO NOT need to install the libraries below.

- `Cicero <https://cole-trapnell-lab.github.io/cicero-release/docs/#installing-cicero>`_
