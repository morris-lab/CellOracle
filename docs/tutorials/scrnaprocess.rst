.. _scrnaprocess:

Single-cell RNA-seq data preprocessing
======================================

Network analysis and simulation in celloracle require processed scRNA-seq data.
We recommend processing scRNA-seq data using either scanpy or Seurat.
If you want to know general idea of scRNA-seq data processing, please go to `the documentation of scanpy <https://scanpy.readthedocs.io/en/stable/>`_ and `the documentation of Seurat <https://satijalab.org/seurat/vignettes.html>`_ .
This step needs to include data quality check, filtering, normalization, dimensional reduction, clustering and cluster annotation.
If you have preprocessed scRNA-seq data, such as anndata or seurat object, you can skip this part.

.. note::

   several tips


A. scRNA-seq data preprocessing with scanpy
-------------------------------------------
scanpy is a python library for the analysis of scRNA-seq data.

In this tutorial, we use the almost same code as  `one of scanpy tutorials <https://scanpy-tutorials.readthedocs.io/en/latest/paga-paul15.html>`_ with some modification.

Python notebook

.. toctree::

   ../notebooks/03_scRNA-seq_data_preprocessing/scanpy_preprocessing_with_Paul_etal_2015_data


B. scRNA-seq data preprocessing with Seurat
-------------------------------------------

R notebook ... comming soon

.. toctree::

  ../notebooks/03_scRNA-seq_data_preprocessing/
