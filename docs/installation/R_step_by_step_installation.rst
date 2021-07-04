.. _R_step_by_step_installation:

R dependent library installation troubleshooting
================================================


`igraph <https://igraph.org/r/>`_
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Please install ``igraph`` with the following r-script or `the author's instruction <https://igraph.org/r/>`_ .

In R console,

.. code-block:: r

   install.packages("igraph")

If you get an error during installation, please check compilers. `This GitHub issue page <https://github.com/igraph/rigraph/issues/275>`_ is helpful.


`linkcomm <https://cran.r-project.org/web/packages/linkcomm/index.html>`_
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Please install ``linkcomm`` with the following r-script or `the author's instruction <https://cran.r-project.org/web/packages/linkcomm/index.html>`_ .

In R console,

.. code-block:: r

   install.packages("linkcomm")


`rnetcarto <https://github.com/cran/rnetcarto/blob/master/src/rgraph/README.md>`_
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Please install ``rnetcarto`` with the following r-script or `the author's instruction <https://github.com/cran/rnetcarto/blob/master/src/rgraph/README.md>`_ .
``rnetcarto`` requires  `the GNU scientific libraries <https://www.gnu.org/software/gsl/>`_ .

If you use ubuntu, you can install the GNU scientific libraries as follows.


::

    sudo apt-get install libgsl-dev


In R console,

.. code-block:: r

   install.packages("rnetcarto")
