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

1. ``rnetcarto`` requires  `the GNU scientific libraries <https://www.gnu.org/software/gsl/>`_ .
 Please make sure your computational environment has GNU scientific libraries.
 In short, you can quickly install GNU scientific libraries as follows.

 For Linux users, please run the code below in the terminal.
 ::

     sudo apt-get install libgsl-dev
 For MacOS users, please run the code below in the terminal.
 ::

    brew install gsl

 For more infoamation about GNU scientific library installation, please look at `here <https://gist.github.com/TysonRayJones/af7bedcdb8dc59868c7966232b4da903>`_ .


2. ``rnetcarto`` is temporaliry unavailable through CRAN (at Jan 5 2021).
 But you can install ``rnetcarto`` with `source file  <https://cran.r-project.org/src/contrib/Archive/rnetcarto/>`_ .

 In short, please run the code below in R console.

 .. code-block:: r

    install.packages("https://cran.r-project.org/src/contrib/Archive/rnetcarto/rnetcarto_0.2.4.tar.gz",
                     repos = NULL, type = "source", configure.args = '--host=host')
