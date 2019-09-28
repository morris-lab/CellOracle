API
===


Command Line API
----------------

CellOracle have a command line API.
This command can be used to convert scRNA-seq data.
If you have a scRNA-seq data which was processed with Seurat and saved as Rds file, you can use the following command to make anndata from Seurat object.
The anndata produced by this command can be used for input of celloracle.


::

    seuratToAnndata YOUR_SEURAT_OBJECT.Rds OUTPUT_PATH





Python API
-----------

.. toctree::
    celloracle
