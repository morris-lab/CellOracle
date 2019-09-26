API
===


Command Line API
----------------

CellOracle have one command line API. This command can be used to convert scRNA-seq data.
If you have a scRNA-seq data which was processed with seurat and saved as Rds file, you can use the following command to make anndata from seurat object.
The anndata produced by this command can be used for an input of celloracle. Please go to step3 and step4 of the tutorial to see an example.

::

    seuratToAnndata YOUR_SEURAT_OBJECT.Rds OUTPUT_PATH





Python API
-----------

.. toctree::
    celloracle
