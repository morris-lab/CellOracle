.. _base_grn:

Base GRN input data preparation
===============================


Overview
--------


There are several options for CellOracle base-GRN construction.
Here is the illustration of base-GRN construction workflow.


.. image:: ./base_GRN_workflow.png
   :scale: 23%



- In this documentation, we introduce details of option1 and option2.
- Option3 uses promoter database for the input of base-GRN construction. We provide pre-built promoter base-GRN for 10 species. You can load this base GRN using celloracle data loading function.
- In option4, any TF-target gene list can be used as a base-GRN. Here is an example notebook [link].


Option1. Data preprocessing of scATAC-seq data
----------------------------------------------

If you have scATAC-seq data, you can use scATAC-seq data to obtain the accessible promoter/enhancer DNA sequence.
To prepare input data of base-GRN construction, we need to get the accessible promoter/enhancer DNA sequence from scATAC-seq data.
This is just a data preparation step, and it uses Cicero and does NOT use celloracle.
Distal cis-regulatory elements can be picked up using  `Cicero <https://cole-trapnell-lab.github.io/cicero-release/docs/#installing-cicero>`_ .

Although we provide an example notebook here, you can analyze your data with Cicero in a different way if you are familiar with Cicero.
If you have a question about Cicero, please read `the documentation of Cicero <https://cole-trapnell-lab.github.io/cicero-release/>`_ for the detailed usage.

scATAC-seq analysis with Cicero
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The jupyter notebook files and data used in this tutorial are available `here <https://github.com/morris-lab/CellOracle/tree/master/docs/notebooks/01_ATAC-seq_data_processing>`_ .

R notebook

.. toctree::

   ../notebooks/01_ATAC-seq_data_processing/option1_scATAC-seq_data_analysis_with_cicero/01_atacdata_analysis_with_cicero_and_monocle3



TSS annotation
^^^^^^^^^^^^^^

We got active promoter / enhancer peaks. CellOracle annotate these peaks to get annotated promoter / enhancer peaks for target all genes.


The jupyter notebook files and data used in this tutorial are available `here <https://github.com/morris-lab/CellOracle/tree/master/docs/notebooks/01_ATAC-seq_data_processing>`_ .

Python notebook

.. toctree::

   ../notebooks/01_ATAC-seq_data_processing/option1_scATAC-seq_data_analysis_with_cicero/02_preprocess_peak_data



Option2. Data preprocessing of bulk ATAC-seq data
-------------------------------------------------
Bulk DNA-seq data can be used to get the accessible promoter/enhancer sequences.


The jupyter notebook files and data used in this tutorial are available `here <https://github.com/morris-lab/CellOracle/tree/master/docs/notebooks/01_ATAC-seq_data_processing>`_ .

Python notebook

.. toctree::

   ../notebooks/01_ATAC-seq_data_processing/option2_Bulk_ATAC-seq_data/01_preprocess_Bulk_ATAC_seq_peak_data
