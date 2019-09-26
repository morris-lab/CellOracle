.. _atac:

ATAC-seq data preprocessing
===========================

Process scATAC-seq data (or bulk ATAC-seq data) to get open accessible promoter/enhancer DNA sequence.
If you have a scATAC-seq data, proximal and distal cis-regulatory elements from TSS site can be picked up using `Cicero <https://cole-trapnell-lab.github.io/cicero-release/docs/#installing-cicero>`_ . If you have bulk ATAC-seq data, open accessible DNA elements around TSS site will be picked up.

A. Extract TF binding information from scATAC-seq data
----------------------------------------------------
If you want to use scATAC-seq data, you can start from Cicero analysis to get information of distal cis-regulatory elements.
This step use Cicero and do not use celloracle. Please refer to `the documentation of Cicero <https://cole-trapnell-lab.github.io/cicero-release/>`_ for the detailed usage.

R notebook

.. toctree::

   ../notebooks/01_ATAC-seq_data_processing/option1_scATAC-seq_data_analysis_with_cicero/01_atacdata_to_cicero

Next, the results of Cicero analysis will be processed with celloracle to make TSS annotations.

Python notebook

.. toctree::

   ../notebooks/01_ATAC-seq_data_processing/option1_scATAC-seq_data_analysis_with_cicero/02_preprocess_peak_data



B. Extract TF binding information from bulk ATAC-seq data or Chip-seq data
--------------------------------------------------------------------------
Instead of scATAC-seq data bulk DNA-seq data can be used.

Python notebook

.. toctree::

   ../notebooks/01_ATAC-seq_data_processing/option2_Bulk_ATAC-seq_data/01_preprocess_Bulk_ATAC_seq_peak_data
