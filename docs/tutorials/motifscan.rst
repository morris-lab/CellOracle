.. _motifscan:

Transcription factor binding motif scan
=======================================


In the previous section, we got accessible Promoter/enhancer DNA regions using ATAC-seq data.
Next, we will obtain a base-GRN by scanning the regulatory genomic sequences for TF-binding motifs.
In the later GRN inference process, this list will be used to define potential regulatory connections.

The jupyter notebook files and data used in this tutorial are available `here <https://github.com/morris-lab/CellOracle/tree/master/docs/notebooks/02_motif_scan>`_ .


Scan DNA sequences searching for TF binding motifs
--------------------------------------------------

Python notebook

.. toctree::

   ../notebooks/02_motif_scan/02_atac_peaks_to_TFinfo_with_celloracle_20200801



How to use different motif data
===============================

Celloracle provides several dafault motifs. If you don't enter motif data, celloracle automatically load default motifs for your species.
In most case, you don't need to prepare TF binding motifs yourself.

But you can use another motif data.

gimmemotifs motif data
----------------------

Here is the notebook describing how to load a motif data from gimmemotifs database.
https://github.com/morris-lab/CellOracle/blob/master/docs/notebooks/02_motif_scan/motif_data_preparation/01_How_to_load_gimmemotifs_motif_data.ipynb

CellOracle motif dataset generated from CisBP version2 database
------------------------------------------------------------

Here is the notebook describing how to load a motif data from CisBP version2 database.
https://github.com/morris-lab/CellOracle/blob/master/docs/notebooks/02_motif_scan/motif_data_preparation/02_How_to_load_CisBPv2_motif_data.ipynb



How to create custom motif data
-------------------------------
We can creaet motif data by ourself.
Here is an example code.
https://github.com/morris-lab/CellOracle/blob/master/docs/notebooks/02_motif_scan/motif_data_preparation/03_How_to_make_custom_motif.ipynb
