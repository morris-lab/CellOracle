.. _check_installation:


Check python library installation status
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can check the installed library version as follows.

In python console,

.. code-block:: Python

   import celloracle as co
   co.check_python_requirements()


Check R library installation status
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Please make sure that all R libraries are installed using the following function.

.. code-block:: Python

  import celloracle as co
  co.test_R_libraries_installation()


The following message will be shown when all R libraries are appropriately installed.

| R path: /usr/lib/R/bin/R
| checking R library installation: igraph -> OK
| checking R library installation: linkcomm -> OK
| checking R library installation: rnetcarto -> OK


The first line above is your R path. If you want to use another R program installed at a different place, please set a new R path with the following command.

.. code-block:: Python

   co.network_analysis.set_R_path("ENTER YOUR R PATH HERE")
