.. _scrnaprocess:

scRNA-seq data preparation
==========================

Overview
--------

Before starting a CellOracle analysis, the scRNA-seq data must be preprocessed.
Please prepare the scRNA-seq data as an `anndata object` using Scanpy.

.. note::
   Scanpy is a python toolkit for scRNA-seq data analysis. If you are new to Scanpy, please read the documentation to learn it in advance.

    - scanpy documentation: https://scanpy.readthedocs.io/en/stable/
    - anndata documentation: https://anndata.readthedocs.io/en/latest/

.. warning::
   In this section, we will introduce an example of how to prepare the **input data** for CellOracle analysis.
   **This is NOT the CellOracle analysis itself.** We do NOT use CellOracle in this scRNA-seq data preprocessing steps.


A. scRNA-seq data preprocessing with Scanpy
-------------------------------------------
Please download the notebook from `here <https://github.com/morris-lab/CellOracle/blob/master/docs/notebooks/03_scRNA-seq_data_preprocessing/scanpy_preprocessing_with_Paul_etal_2015_data.ipynb>`_ .
Or please click below to view the content.

.. toctree::

   ../notebooks/03_scRNA-seq_data_preprocessing/scanpy_preprocessing_with_Paul_etal_2015_data


B. scRNA-seq data preprocessing with Seurat
-------------------------------------------

R notebook ... comming in the future update.


.. note::

   If you used ``Seurat`` for preprocessing, you will need to convert your Seurat object into an anndata object.
   CellOracle has a python API and command-line API to help users with this data conversion.
   Please go to the documentation of CellOracle's API documentation for more information.
