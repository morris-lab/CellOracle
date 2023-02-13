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

For more information, please read our paper: `Dissecting cell identity via network inference and in silico gene perturbation <https://www.nature.com/articles/s41586-022-05688-9>`_


.. raw:: html
   :file: test2.html


News
====

Please look `Changelog page <https://morris-lab.github.io/CellOracle.documentation/changelog/index.html>`_ for all updates history of CellOracle package.


- 02/08/2023: CellOracle paper has been published in `Nature article <https://www.nature.com/articles/s41586-022-05688-9>`_ ! It is also highlighted in `News and Views <https://www.nature.com/articles/d41586-023-00251-6>`_ . Thank you very much to all contributors and users for your help!!!

- 02/02/2023: We released celloracle version 0.10.15. We fixed a bug in Motif data loading function.

- 01/03/2023: We released celloracle version 0.10.14. We added a Xenopus Laevis data.

- 10/10/2022: We released celloracle version 0.10.12. We added a new function to check distribution range of simulated value. Tutorial is also updated.





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


.. |GitHub Workflow Status| image:: https://img.shields.io/github/actions/workflow/status/morris-lab/CellOracle/build_check.yml?branch=master
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
