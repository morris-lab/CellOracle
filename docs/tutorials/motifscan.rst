.. _motifscan:

Transcription factor binding motif scan
=======================================


In the previous section, we identified accessible promoter/enhancer DNA regions using ATAC-seq data.
Next, we will construct the base GRN by scanning the regulatory genomic sequences for TF-binding motifs.
The base GRN, potential TF-target gene connection list, will be used in the later GRN inference step.

The jupyter notebook files and data used in this tutorial are available `here <https://github.com/morris-lab/CellOracle/blob/master/docs/notebooks/02_motif_scan/02_atac_peaks_to_TFinfo_with_celloracle_20200801.ipynb>`_ .


Scan DNA sequences searching for TF binding motifs
--------------------------------------------------

Python notebook

.. toctree::

   ../notebooks/02_motif_scan/02_atac_peaks_to_TFinfo_with_celloracle_20200801



How to use custom motif data
===============================

CellOracle provides several default motif datasets.
If you do not specify the motif data, CellOracle automatically loads the default motifs for your species.
In most cases, you will not need to prepare your own TF binding motif dataset.

However, you do have the option to define or customize a motif dataset for your analysis.

gimmemotifs motif data
----------------------

Here is the notebook demonstrating how to load a motif data from gimmemotifs database.
https://github.com/morris-lab/CellOracle/blob/master/docs/notebooks/02_motif_scan/motif_data_preparation/01_How_to_load_gimmemotifs_motif_data.ipynb

CellOracle motif dataset generated from the CisBP version2 database
------------------------------------------------------------

Here is the notebook describing how to load a motif data from `CisBP version 2 database <http://cisbp.ccbr.utoronto.ca/index.php>`_ .
https://github.com/morris-lab/CellOracle/blob/master/docs/notebooks/02_motif_scan/motif_data_preparation/02_How_to_load_CisBPv2_motif_data.ipynb



How to create custom motif data
-------------------------------
We can also create new custom motif datasets defined by yourself.
Here is an example notebook.
https://github.com/morris-lab/CellOracle/blob/master/docs/notebooks/02_motif_scan/motif_data_preparation/03_How_to_make_custom_motif.ipynb
