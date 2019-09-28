.. _simulation:

Simulation with GRNs
====================

``celloracle`` leverage GRNs to simulate signal propagation inside a cell.
We can estimate the effect of gene perturbation by the simulation with GRNs.

``celloracle`` leverage GRNs to simulate signal propagation inside a cell. We can estimate the effect of gene perturbation by the simulation with GRNs.

Besides, we will combine the simulation for the signal propagation and the simulation of cell state transition, which is performed by a python library for RNA-velocity analysis, ``velocyto`` .
A series of simulations visualizes a complex system in which TF controls enormous target genes to determines cell fate.

In short, celloracle provides insight into the regulatory mechanism of cell state from the viewpoint of TF and GRN.


Python notebook

.. toctree::

   ../notebooks/05_simulation/Gata1_KO_simulation_with_with_Paul_etal_2015_data
