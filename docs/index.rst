.. celloracle documentation master file, created by
   sphinx-quickstart on Thu Sep 19 17:42:52 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to CellOracle's documentation!
======================================

CellOracle is a python library for the in silico gene perturbation analysis using single-cell omics data and gene regulatory network models.

Source code is available at `CellOracle  GitHub repository <https://github.com/morris-lab/CellOracle>`_

For more information, please read our bioRxiv preprint: `CellOracle: Dissecting cell identity via network inference and in silico gene perturbation <https://www.biorxiv.org/content/10.1101/2020.02.17.947416v3>`_


.. warning::
   CellOracle is still under development. It is currently a beta version. Functions in this package may change in a future release.

News
====

Please look `Changelog page <https://morris-lab.github.io/CellOracle.documentation/changelog/index.html>`_ for all updates history of CellOracle package.

- 04/13/2022: We refactored the code to simplify the installation process and calculation speed. The new version, celloracle 0.10.0, does not require any R packages anymore. Also, codes were updated to comply with the recent numpy updates.

- 01/29/2022: Docker image was updated.

- 01/21/2022: We have changed the default color scheme for the inner product / perturbation score visualization in the in silico gene perturbation module. Please look at the tutorial notebook for detail.



Contents
========


.. toctree::
   :hidden:

   self


.. toctree::
   :maxdepth: 2

   installation/index

.. toctree::
   :maxdepth: 2

   tutorials/index


.. toctree::
   :maxdepth: 2

   modules/index

.. toctree::

   changelog/index

.. toctree::

   license/index

.. toctree::

   citation/index

.. toctree::

   contact/index


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
