.. _motifscan:

Transcription factor binding motif scan
======================================



We identified accessible Promoter/enhancer DNA regions using ATAC-seq data.
Next, we will obtain a list of TFs for each target gene by scanning the regulatory genomic sequences for TF-binding motifs.
In the later GRN inference process, this list will be used to define potential regulatory connections.

The jupyter notebook files and data used in this tutorial are available `here <https://github.com/morris-lab/CellOracle/tree/master/docs/notebooks/02_motif_scan>`_ .


Scan DNA sequences searching for TF binding motifs
--------------------------------------------------

Python notebook

.. toctree::

   ../notebooks/02_motif_scan/02_atac_peaks_to_TFinfo_with_celloracle_20200801



[Optional step 1] How to use different motif data
-------------------------------------------------

Celloracle motif analysis pipeline provides several dafault motifs. If you don't enter motif data, celloracle automatically load default motifs for your species.
In that case, you don't need to make TF binding motifs.
But also, you can pick up TF binding motifs by yourself.
Here, we introduce how to find and load motifs for celloracle motif analysis.
Some codes for custom motif were based on the suggestion in `this post <https://https://github.com/morris-lab/CellOracle/issues/9>`_ . from KyleFerchen. Thank you Kyle!


.. toctree::

   ../notebooks/02_motif_scan/motif_data_preparation/01_How_to_load_motif_data




[Optional step 2] How to Make custom motifs for celloracle motif analysis
-------------------------------------------------------------------------
If you cannot find an appropriate motif dataset for your analysis and want to it by yourself. You can follow the instruction below.
We introduce an example way to make motifs using `CisDB TF binding database <https://http://cisbp.ccbr.utoronto.ca>`_ .



.. toctree::

   ../notebooks/02_motif_scan/motif_data_preparation/02_How_to_prepare_custom_motif_data
