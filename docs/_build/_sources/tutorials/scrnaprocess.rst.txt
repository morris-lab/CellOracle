.. _scrnaprocess:

Single-cell RNA-seq data preprocessing
======================================

Network analysis and simulation in celloracle will be performed using scRNA-seq data. The scRNA-seq data should include the compontnes below.
 - Gene expression matrix; mRNA counts before scaling and transformation.
 - Clustering results.
 - Dimensional reduction results.

In addition to these minimum requirements, we highly recommend doing these analyses below in the preprocessing step.
 - Data quality check and Cell/Gene filtering.
 - Normalization
 - Identification of highly variable genes

We recommend processing scRNA-seq data using either scanpy or Seurat.
If you are not familiar with the general workflow of scRNA-seq data processing, please go to `the documentation of scanpy <https://scanpy.readthedocs.io/en/stable/>`_ and `the documentation of Seurat <https://satijalab.org/seurat/vignettes.html>`_ .

If you already have preprocessed scRNA-seq data, which includes necessary information above, you can skip this part.



A. scRNA-seq data preprocessing with scanpy
-------------------------------------------
``scanpy`` is a python library for the analysis of scRNA-seq data.

In this tutorial, we introduce an example of scRNA-seq preprocessing for celloracle with ``scanpy``.
We wrote the notebook based on `one of scanpy's tutorials <https://scanpy-tutorials.readthedocs.io/en/latest/paga-paul15.html>`_ with some modification.


Python notebook

.. toctree::

   ../notebooks/03_scRNA-seq_data_preprocessing/scanpy_preprocessing_with_Paul_etal_2015_data


B. scRNA-seq data preprocessing with Seurat
-------------------------------------------

R notebook ... comming soon


.. note::

   If you use ``Seurat`` for preprocessing, you need to convert the scRNA-seq data (Seurat object) into anndata to analyze the data with ``celloracle``.
   ``celloracle`` have python API and command-line API to convert Seurat object into anndata.
   Please go to the documentation of celloracle API for more information.
