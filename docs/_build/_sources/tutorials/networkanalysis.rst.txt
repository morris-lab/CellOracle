.. _networkanalysis:

Network analysis
===========================

``celloracle`` import scRNA-seq dataset and TF binding information to find active regulatory connections for all genes, generating sample-specific GRNs.

The inferred GRN is analyzed with several network algorithms to get various network scores. The network score is useful to identify key regulatory genes.

celloracle reconstructs GRN for each cluster enabling us to compare GRNs each other. It is also possible to analyze how the GRN change over differentiation.
The dynamics of the GRN structure provide us insight into the
context-dependent regulatory mechanism.

Python notebook

.. toctree::

   ../notebooks/04_Network_analysis/Network_analysis_with_with_Paul_etal_2015_data
