.. _python_step_by_step_installation:

Python dependent library installation troubleshooting
=====================================================


`velocyto <http://velocyto.org/velocyto.py/install/index.html>`_
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you failed CellOracle installation because of velocyto installation error, please try to install velocyto using the following commands or please follow `the author's instruction <http://velocyto.org/velocyto.py/install/index.html>`_ .

::

    conda install numpy scipy cython numba matplotlib scikit-learn h5py click pysam llvm louvain

Then

::

    pip install velocyto

Some MacOS users have encountered compiler errors during the installation of velocyto.
Several different errors have been reported.
You may find the solution for each of these in the links below.

- `Solution 1: Install Xcode <https://developer.apple.com/xcode/>`_. Please try this first.
- `Solution 2: Install macOS_SDK_headers <https://stackoverflow.com/a/53057706/10641716>`_. This solution is needed in addition to Solution-1 if your OS is macOS Mojave.
- `Solution 3 <https://github.com/morris-lab/CellOracle/issues/3>`_. This is the solution reported by a CellOracle user. Thank you very much!
- `Other solutions on Velocyto GitHub issue page <https://github.com/velocyto-team/velocyto.py/issues?q=>`_

`scanpy <https://scanpy.readthedocs.io/en/stable/installation.html>`_
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you failed CellOracle installation because of Scanpy installation error, please try to install Scanpy using the following commands or `the author's instruction <https://scanpy.readthedocs.io/en/stable/installation.html>`_ .

::

    conda install scanpy


Install gimmemotifs with conda
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you failed installation in the gimmemotifs installation step, please install gimmemotifs solely first.
Also, if your jupyter notebook kernel keep dying during TF motif scan, it means your gimmemotifs installation has a problem.
In that case, please re-install gimmemotifs using conda.

::

    conda install -c bioconda gimmemotifs




If you encounter an error related to "certifi"
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you encounter the following error, it means the error is caused by versiom mismatch of "certifi".
See `this page  <https://stackoverflow.com/questions/50129762/graphlab-create-2-1-installation-fails-to-uninstall-certifi-a-distutils-insta>`_. for more information.

::

    ERROR: Cannot uninstall 'certifi'. It is a distutils installed project and thus we cannot accurately determine which files belong to it which would lead to only a partial uninstall.

In this case, please add "--ignore-installed certifi " to the installation command.

::

    pip install git+https://github.com/morris-lab/CellOracle.git --ignore-installed certifi
