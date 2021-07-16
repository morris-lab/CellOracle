.. _scrnaprocess:

scRNA-seq data preparation
==========================

Overview
--------

In advance to CellOrale analysis, scRNA-seq data should be processed.
Please prepare scRNA-seq data as an `anndata` using `scanpy`.

.. note::
   `scanpy` is a python toolkit for scRNA-seq data analysis. If you are new to scanpy, pelase read the documentation to learn it in advance.

    - scanpy documentation: https://scanpy.readthedocs.io/en/stable/
    - anndata documentation: https://anndata.readthedocs.io/en/latest/

.. warning::
   In this section, we intend to introduce an example of how to prepare the **input data** for CellOracle analysis.
   **This is NOT the CellOracle analysis itself.** We do NOT use celloracle in this notebook.


A. scRNA-seq data preprocessing with scanpy
-------------------------------------------
Please download notebooks from `here <https://github.com/morris-lab/CellOracle/blob/master/docs/notebooks/03_scRNA-seq_data_preprocessing/scanpy_preprocessing_with_Paul_etal_2015_data.ipynb>`_ .
Or please click below to view the content.

.. toctree::

   ../notebooks/03_scRNA-seq_data_preprocessing/scanpy_preprocessing_with_Paul_etal_2015_data


B. scRNA-seq data preprocessing with Seurat
-------------------------------------------

R notebook ... comming in the future update.


.. note::

   If you use ``Seurat`` for preprocessing, you need to convert the scRNA-seq data (Seurat object) into anndata to analyze the data with ``celloracle``.
   ``celloracle`` has a python API and command-line API to convert a Seurat object into an anndata.
   Please go to the documentation of celloracle's API documentation for more information.
