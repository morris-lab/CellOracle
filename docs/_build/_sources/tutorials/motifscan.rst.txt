.. _motifscan:

Transcription factor binding motif scan
======================================



We identified open-accessible Promoter/enhancer DNAs using ATAC-seq data.
Next, we scan the regulatory genomic sequence searching for the TF-binding motifs to get the list of TFs for each target gene.
In the later GRN inference process, this TF list per target gene will be used as a potential regulatory connection.

Python notebook

.. toctree::

   ../notebooks/02_motif_scan/02_atac_peaks_to_TFinfo_with_celloracle_190901
