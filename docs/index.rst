.. celloracle documentation master file, created by
   sphinx-quickstart on Thu Sep 19 17:42:52 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to CellOracle's documentation!
======================================
|GitHub Workflow Status| |PyPI| |PyPI - Python Version| |PyPI - Wheel|
|Downloads| |Docker Pulls|


CellOracle is a python library for the in silico gene perturbation analysis using single-cell omics data and gene regulatory network models.

Source code is available at `CellOracle  GitHub repository <https://github.com/morris-lab/CellOracle>`_ and `PyPI <https://pypi.org/project/celloracle/>`_.

For more information, please read our bioRxiv preprint: `CellOracle: Dissecting cell identity via network inference and in silico gene perturbation <https://www.biorxiv.org/content/10.1101/2020.02.17.947416v3>`_


.. warning::
   CellOracle is still under development. It is currently a beta version. Functions in this package may change in a future release.

News
====

Please look `Changelog page <https://morris-lab.github.io/CellOracle.documentation/changelog/index.html>`_ for all updates history of CellOracle package.


- 10/10/2022: We released celloracle version 0.10.12. We added a new function to check distribution range of simulated value. Tutorial is also updated.

- 09/27/2022: We released celloracle version 0.10.11. We added a new function to change cluster color.

- 09/09/2022: We released celloracle version 0.10.10. We updated some functions in the network score calculation and visualization.

- 07/21/2022: We released celloracle version 0.10.8. Simulation tutorial notebook is updated; We added a guidance to find best vm parameter for the simulation result visualization.



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


.. |GitHub Workflow Status| image:: https://img.shields.io/github/workflow/status/morris-lab/CellOracle/build_and_test
   :target: https://github.com/morris-lab/CellOracle/actions/workflows/build_check.yml
.. |PyPI| image:: https://img.shields.io/pypi/v/celloracle?color=blue
   :target: https://pypi.org/project/celloracle/
.. |PyPI - Python Version| image:: https://img.shields.io/pypi/pyversions/celloracle
   :target: https://pypi.org/project/celloracle/
.. |PyPI - Wheel| image:: https://img.shields.io/pypi/wheel/celloracle
   :target: https://pypi.org/project/celloracle/
.. |Downloads| image:: https://static.pepy.tech/personalized-badge/celloracle?period=total&units=international_system&left_color=grey&right_color=orange&left_text=Downloads
   :target: https://pepy.tech/project/celloracle
.. |Docker Pulls| image:: https://img.shields.io/docker/pulls/kenjikamimoto126/celloracle_ubuntu?color=orange
   :target: https://hub.docker.com/r/kenjikamimoto126/celloracle_ubuntu
