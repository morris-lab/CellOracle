.. _atac:

ATAC-seq data preprocessing
===========================

In this step, we process scATAC-seq data (or bulk ATAC-seq data) to obtain the accessible promoter/enhancer DNA sequence.
We can get the active proximal promoter/enhancer genome sequences by picking up the ATAC-seq peaks that exist around the transcription starting site (TSS).
Distal cis-regulatory elements can be picked up using  `Cicero <https://cole-trapnell-lab.github.io/cicero-release/docs/#installing-cicero>`_ .
Cicero analyzes scATAC-seq data to calculate a co-accessible score between peaks.
We can identify cis-regulatory elements using Cicero's co-access score and TSS information.

If you have bulk ATAC-seq data instead of scATAC-data, we'll get only the proximal promoter/enhancer genome sequences.



A. Extract TF binding information from scATAC-seq data
------------------------------------------------------
If you have scATAC-seq data, you can get information on the distal cis-regulatory elements.
This step uses Cicero and does not use celloracle.
You need to get co-accessibility table in this analysis.
Although we provide an example notebook here, you can analyze your data with Cicero in a different way if you are familiar with Cicero.
If you have a question about Cicero, please read `the documentation of Cicero <https://cole-trapnell-lab.github.io/cicero-release/>`_ for the detailed usage.

scATAC-seq analysis with Cicero and Monocle3
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


The jupyter notebook files and data used in this tutorial are available `here <https://github.com/morris-lab/CellOracle/tree/master/docs/notebooks/01_ATAC-seq_data_processing>`_ .

R notebook

.. toctree::

   ../notebooks/01_ATAC-seq_data_processing/option1_scATAC-seq_data_analysis_with_cicero/01_atacdata_analysis_with_cicero_and_monocle3



TSS annotation
^^^^^^^^^^^^^^


The jupyter notebook files and data used in this tutorial are available `here <https://github.com/morris-lab/CellOracle/tree/master/docs/notebooks/01_ATAC-seq_data_processing>`_ .

Python notebook

.. toctree::

   ../notebooks/01_ATAC-seq_data_processing/option1_scATAC-seq_data_analysis_with_cicero/02_preprocess_peak_data



B. Extract TF binding information from bulk ATAC-seq data or Chip-seq data
--------------------------------------------------------------------------
Bulk DNA-seq data can be used to get the accessible promoter/enhancer sequences.


The jupyter notebook files and data used in this tutorial are available `here <https://github.com/morris-lab/CellOracle/tree/master/docs/notebooks/01_ATAC-seq_data_processing>`_ .

Python notebook

.. toctree::

   ../notebooks/01_ATAC-seq_data_processing/option2_Bulk_ATAC-seq_data/01_preprocess_Bulk_ATAC_seq_peak_data
