.. _motifscan:

Transcription factor binding motif scan
======================================



We identified accessible Promoter/enhancer DNA regions using ATAC-seq data.
Next, we will obtain a list of TFs for each target gene by scanning the regulatory genomic sequences for TF-binding motifs.
In the later GRN inference process, this list will be used to define potential regulatory connections.

Python notebook

.. toctree::

   ../notebooks/02_motif_scan/02_atac_peaks_to_TFinfo_with_celloracle_190901
