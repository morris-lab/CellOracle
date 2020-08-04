.. _networkanalysis:

Network analysis
===========================

``celloracle`` imports the scRNA-seq dataset and TF binding information to find active regulatory connections for all genes, generating sample-specific GRNs.

The inferred GRN is analyzed with several network algorithms to get various network scores. The network score is useful to identify key regulatory genes.

Celloracle reconstructs a GRN for each cluster, enabling us to compare GRNs to each other. It is also possible to analyze how the GRN changes over differentiation.
The dynamics of the GRN structure can provide us insight into the context-dependent regulatory mechanisms.


The jupyter notebook files and data used in this tutorial are available `here <https://github.com/morris-lab/CellOracle/tree/master/docs/notebooks/04_Network_analysis>`_ .

Python notebook

.. toctree::

   ../notebooks/04_Network_analysis/Network_analysis_with_Paul_etal_2015_data
