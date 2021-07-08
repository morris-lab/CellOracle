.. _tutorial:

Tutorial
========

This tutorial aims to introduce how to use CellOracle functions using the demo dataset.
Once you get used to CellOracle codes, please replace demo data with your data to investigate it.

What the tutorial covers
------------------------

- :doc:`networkanalysis`: This notebook introduces how to construct sample-specific GRN models. It also contains examples of network analysis with graph theory.
- :doc:`simulation` : This notebook performs in silico gene perturbation analysis using GRN models.

Demo dataset is available in the tutorial notebooks above. You can try CellOracle quickly.
But If you want to apply CellOracle to your scRNA-seq or scATAC dataset, please refer to the following tutorials to know how to prepare input data.

- :doc:`scrnaprocess`: This notebook explains preprocessing steps for scRNA-seq data.
- :doc:`base_grn`: This tutorial explains how to prepare input data for TF motif scan.
- :doc:`motifscan`: This tutorial describes the TF motif scan pipeline for base-GRN construction.

Prerequisites
-------------

- This tutorial assumes that you have some Python programming experience. In particular, we assume you are familiar with Python data science libraries: jupyter, pandas, and matplotlib.

- Also, this tutorial assume that you are familiar with basic scRNA-seq data analysis. In particular, we assume you have some experience of scRNA-seq analysis using `Scanpy and Anndata <https://scanpy.readthedocs.io/en/stable/>`_ , which is a python toolkit for single-cell analysis. You can use scRNA-seq data processed with `Seurat <https://satijalab.org/seurat/>`_ . But the Seurat data need to be converted into Anndata format in advance to CellOracle analysis. See this  `tutorial  <https://xxx.com>`_ for detail.

- CellOracle provides pre-build base-GRN, and it is not necessary to construct custom base-GRN. But if you want to construct custom base-GRN from your scATAC-seq data, we recommend using `Cicero <https://cole-trapnell-lab.github.io/cicero-release/>`_ . In this case, please get used to Cicero, basic scATAC-seq data analysis, and TF motif analysis in advance to start constructing base-GRN.



Scripts and data
----------------

The jupyter notebook files and data used in this tutorial are available `here <https://github.com/morris-lab/CellOracle/tree/master/docs/notebooks>`_ .



Getting started
----------------

If you run CellOracle for the first time, please start with the :doc:`networkanalysis`.
We provide demo scRNA-seq dataset and base-GRN data as follows.
You can load these data using the CellOracle data loading function.

- scRNA-seq data: Hematopoiesis dataset published by `Paul et al (2015) <https://www.sciencedirect.com/science/article/pii/S0092867415014932?via%3Dihub>`_ .
- Base-GRN: Base-GRN generated from `Mouse sci-ATAC-seq atlas dataset <https://atlas.gs.washington.edu/mouse-atac/>`_ .

You can easily start CellOracle analysis with this dataset.
You can reproduce hematopoiesis network analysis and perturbation simulation results that are shown in `our bioarxiv preprint <https://www.biorxiv.org/content/10.1101/2020.02.17.947416v3>`_ .



GRN model construction and Network analysis
------------------------------------------

.. toctree::
   :maxdepth: 2
   :numbered: 1


   networkanalysis

in silico Gene perturbation with GRNs
-------------------------------------

.. toctree::
   :maxdepth: 2
   :numbered: 1

   simulation


Prepare input data
------------------

.. toctree::
   :maxdepth: 2
   :numbered: 1

   scrnaprocess
   base_grn

TF motif scan
-------------

.. toctree::
   :maxdepth: 2
   :numbered: 1

   motifscan
