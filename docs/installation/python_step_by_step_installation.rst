.. _python_step_by_step_installation:

Python dependent library installation troubleshooting
=====================================================


Install `velocyto <http://velocyto.org/velocyto.py/install/index.html>`_
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you failed CellOracle installation because of velocyto installation error, please try to install velocyto with the following commands or `the author's instruction <http://velocyto.org/velocyto.py/install/index.html>`_ .

::

    conda install numpy scipy cython numba matplotlib scikit-learn h5py click pysam llvm louvain

Then

::

    pip install velocyto

It was reported that some compile errors might occur during the installation of velocyto on MacOS.
Various errors were reported, and you need to find the best solution depending on your error.
You may find the solution with these links below.

- `Solution 1: Install Xcode <https://developer.apple.com/xcode/>`_. Please try this first.
- `Solution 2: Install macOS_SDK_headers <https://stackoverflow.com/a/53057706/10641716>`_. This solution is needed in addition to Solution-1 if your OS is macOS Mojave.
- `Solution 3 <https://github.com/morris-lab/CellOracle/issues/3>`_. This is the solution reported by a CellOracle user. Thank you very much!
- `Other solutions on Velocyto GitHub issue page <https://github.com/velocyto-team/velocyto.py/issues?q=>`_

Install `scanpy <https://scanpy.readthedocs.io/en/stable/installation.html>`_
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you failed CellOracle installation because of scanpy installation error, please try to install scanpy with the following commands or `the author's instruction <https://scanpy.readthedocs.io/en/stable/installation.html>`_ .

::

    conda install scanpy



Install other python libraries
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Please install other python libraries below using conda prior to celloracle installation. It might solve some installation errors.

::

    conda install pybedtools pyarrow tqdm joblib jupyter gimmemotifs==0.14.4 genomepy==0.8.4


Install celloracle
^^^^^^^^^^^^^^^^^^
After installing the dependent libraries above, please install CellOracle again.

::

    pip install git+https://github.com/morris-lab/CellOracle.git


If you get error related to "certifi"
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you get the following error, it means the error is caused by versiom mismatch of "certifi".
See `this page  <https://stackoverflow.com/questions/50129762/graphlab-create-2-1-installation-fails-to-uninstall-certifi-a-distutils-insta>`_. for more information.

::

    ERROR: Cannot uninstall 'certifi'. It is a distutils installed project and thus we cannot accurately determine which files belong to it which would lead to only a partial uninstall.

In this case, please add "--ignore-installed certifi " to the installation command.

::

    pip install git+https://github.com/morris-lab/CellOracle.git --ignore-installed certifi
