.. _tutorial:

Tutorial
========

This tutorial demonstrates how to use CellOracle functions with a demo dataset.
Once you are familiar with CellOracle’s workflow, please replace the demo data with your own data to begin your analysis.


What the tutorial covers
------------------------

1. Main CellOracle analysis
^^^^^^^^^^^^^^^^^^^^^^^^

- :doc:`networkanalysis`: This notebook introduces how to construct sample-specific GRN models. It also contains examples of network analyses that use graph theory.
- :doc:`simulation` : This notebook performs in silico gene perturbation analysis using GRN models.


.. note::
   A demo dataset is available for each of the tutorial notebooks above. These datasets allow you to begin exploring CellOracle even if you do not have any data at any step in the analysis pipeline.

2. How to prepare input data
^^^^^^^^^^^^^^^^^^^^^^^^^
We recommend getting started with CellOracle using the provided demo dataset.
When you want to apply CellOracle to your own scRNA-seq or scATAC-seq dataset, please refer to the following tutorials to learn how to prepare input data.

- :doc:`scrnaprocess`: This notebook explains the preprocessing steps for scRNA-seq data.
- :doc:`base_grn`: This tutorial explains how to prepare the input data for TF motif scanning.
- :doc:`motifscan`: This tutorial describes the TF motif scan pipeline for base-GRN construction.

.. note::
   In the input data preparation steps, we introduce single-cell preprocessing workflows using other well-known libraries.
   However, please note these precursor notebooks are not considered part of CellOracle program itself.
   They are intended to help users leverage pre-existing tools in a way that is consistent with CellOracle’s workflow.
   **CellOracle is NOT a pipeline for scRNA-seq or scATAC-seq data preprocessing**.


Prerequisites
-------------

- This tutorial assumes that you have adequate **Python programming experience**. In particular, we assume you are familiar with the following python data science libraries: **jupyter, pandas, and matplotlib**.

- Also, this tutorial assume that you are familiar with basic scRNA-seq data analysis. In particular, we assume you basic understanding of scRNA-seq analysis using `Scanpy and Anndata <https://scanpy.readthedocs.io/en/stable/>`_ , which is a python toolkit for single-cell analysis.

- CellOracle provides pre-build base-GRNs, and it is not necessary to construct a custom base-GRN. However, if you would like to use a custom base-GRN constructed using your scATAC-seq data, we recommend doing basic scATAC-seq data analysis using `Cicero <https://cole-trapnell-lab.github.io/cicero-release/>`_.　In this case, please learn the basic usage Cicero prior to start constructing your base-GRN.

- CellOracle is a python package. For the installation of CellOracle, we recommend using `Anaconda <https://anaconda.org>`_ , If you are not familiar with Anaconda or python environment management, please use `our pre-built docker image <https://morris-lab.github.io/CellOracle.documentation/installation/index.html#docker-image>`_.


Code and data availability
--------------------------


- We provide links for the notebook in each section.

- You can download the demo input data using CellOracle data loading function in the notebooks.

- We provide intermediate files. You may begin the demo at any point in the pipeline.

Getting started
----------------

If you are running CellOracle for the first time, we recommend getting started with the :doc:`networkanalysis` and then, proceed to :doc:`simulation`.
We provide demo scRNA-seq dataset and base-GRN data as follows.
The demo scRNA-seq dataset and base GRN data can be loaded within the notebook.


- scRNA-seq data: Hematopoiesis dataset published by `Paul et al. (2015) <https://www.sciencedirect.com/science/article/pii/S0092867415014932?via%3Dihub>`_ .
- Base-GRN: Base-GRN generated from `the Mouse sci-ATAC-seq Atlas <https://atlas.gs.washington.edu/mouse-atac/>`_ .

Using these tutorials, you can easily start CellOracle analysis with this dataset.
You can reproduce hematopoiesis network analysis and perturbation simulation results that are shown in `our bioRxiv preprint <https://www.biorxiv.org/content/10.1101/2020.02.17.947416v3>`_ .


Index
-----

GRN model construction and network analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2
   :numbered: 1


   networkanalysis

In silico gene perturbation with GRNs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2
   :numbered: 1

   simulation


Prepare input data
^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2
   :numbered: 1

   scrnaprocess
   pseudotime
   base_grn

TF motif scan for base-GRN construction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2
   :numbered: 1

   motifscan
