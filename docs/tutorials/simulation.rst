.. _simulation:

Simulation with GRNs
====================

``celloracle`` leverage GRNs to simulate signal propagation inside a cell.
We can estimate the effect of gene perturbation by the simulation with GRNs.

Additonally, we will combine the signal propagation simulation with a cell state transition simulation. The latter simulation is performed by a python library for RNA-velocity analysis, called ``velocyto`` .
This analysis may provide an insight into a complex system how TF controls enormous target genes to determines cell fate.


Python notebook

.. toctree::

   ../notebooks/05_simulation/Gata1_KO_simulation_with_with_Paul_etal_2015_data
