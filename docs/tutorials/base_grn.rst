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

Here, we introduce an example method to extract active promoter / enhancer peaks from scATAC-seq data using Cicero.

.. note::
   Cicero is a R package for scATAC-seq data analysis. It can pick up distal cis-regulatory elements in scATAC-seq data.

.. warning::
   - Here, we intend to introduce an example of how to prepare input data. **This is not CellOracle analysis. We do NOT use celloracle in this step**.
   - This is just an example of data preparation step, you can analyze your data with Cicero in a different way if you are familiar with Cicero. If you have a question about Cicero, please read `the documentation of Cicero <https://cole-trapnell-lab.github.io/cicero-release/>`_ for the detailed usage.
   - If you have a favorite algorithm / software for scATAC-data analysis, you can use totally different software to pick up gene expression regulatory elements.


Step1. scATAC-seq analysis with Cicero
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The jupyter notebook file is available `here <https://github.com/morris-lab/CellOracle/blob/master/docs/notebooks/01_ATAC-seq_data_processing/option1_scATAC-seq_data_analysis_with_cicero/01_atacdata_analysis_with_cicero_and_monocle3.ipynb>`_ .
The R notebook file is available `here <https://github.com/morris-lab/CellOracle/blob/master/docs/notebooks/01_ATAC-seq_data_processing/option1_scATAC-seq_data_analysis_with_cicero/01_atacdata_analysis_with_cicero_and_monocle3.Rmd>`_ .

Or click below to see the contents.

.. toctree::
   :maxdepth: 1

   ../notebooks/01_ATAC-seq_data_processing/option1_scATAC-seq_data_analysis_with_cicero/01_atacdata_analysis_with_cicero_and_monocle3



Step2. TSS annotation
^^^^^^^^^^^^^^^^^^^^^

We can get active promoter / enhancer peaks in step1 above. Next, we will make gene annotations for these peaks.


The jupyter notebook file is available `here <https://github.com/morris-lab/CellOracle/blob/master/docs/notebooks/01_ATAC-seq_data_processing/option1_scATAC-seq_data_analysis_with_cicero/02_preprocess_peak_data.ipynb>`_ .


Or click below to see the contents.

.. toctree::
   :maxdepth: 1

   ../notebooks/01_ATAC-seq_data_processing/option1_scATAC-seq_data_analysis_with_cicero/02_preprocess_peak_data

Once you get the input data, please go to the Motif scan section.


Option2. Data preprocessing of bulk ATAC-seq data
-------------------------------------------------
Bulk DNA-seq data can be used to get the accessible promoter/enhancer sequences.


The jupyter notebook file is available `here <https://github.com/morris-lab/CellOracle/blob/master/docs/notebooks/01_ATAC-seq_data_processing/option2_Bulk_ATAC-seq_data/01_preprocess_Bulk_ATAC_seq_peak_data.ipynb>`_ .


Or click below to see the contents.

.. toctree::
   :maxdepth: 1

   ../notebooks/01_ATAC-seq_data_processing/option2_Bulk_ATAC-seq_data/01_preprocess_Bulk_ATAC_seq_peak_data
