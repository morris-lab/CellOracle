# -*- coding: utf-8 -*-


import pandas as pd
import numpy as np

import sys
import os
import shutil

import scanpy as sc

from ..utility import exec_process

rscript_folder = os.path.abspath(os.path.dirname(__file__))


# this is a function to integrate matrix and meta data and make AnnData object
def _constructAnnData(main_matrix, cell_ids, var_ids, meta_data,
                     categorical_in_meta, raw_data=None):
    '''
    Make anndata manually.

    Args:
        main_matrix (sparce matrix): this should be imported with sc.read_mtx() function.
        cell_ids (numpy.array): array of cell name. the shape should be match with main_matrix
        var_ids (numpy.array): array of variable name. e.g., gene name or peak name in genome.
        the shape should be match with main_matrix
        meta_data (pandas.dataframe): metadata. If the structure of meta_data does not match main_matrix,
        main_matrix will be reconstructed to fit the shape of meta_data.
        categorical_in_meta (array of str): some of meta_data should be categorical data rather than numeric. but such categorical data is might be imported as numeric data sometimes.
        this function convert numeric data into categorical data if you set the colname.

    Returns:
        anndata: anndata.

    '''

    main_matrix = main_matrix.transpose()

    # change dtyope in cluster info mation

    for i in meta_data.columns:
         if i in categorical_in_meta:
            meta_data[i] = meta_data[i].astype(np.object) # change dtype

    # integrate data.
    if main_matrix.shape[0] != meta_data.shape[0]:

        mat = sc.AnnData(main_matrix.X,
                         obs=pd.DataFrame(index=cell_ids),
                         var=pd.DataFrame(index=var_ids))

        cell_ids_in_mata = list(map(lambda x: x in meta_data.index, cell_ids))
        cells_ids_in_meta_index = np.arange(len(cell_ids))[cell_ids_in_mata]

        mat = sc.AnnData(main_matrix.X[cells_ids_in_meta_index,:],
                                          obs=meta_data, var=pd.DataFrame(index=var_ids))

    else:
        mat = sc.AnnData(main_matrix.X,
                         obs=meta_data,
                         var=pd.DataFrame(index=var_ids))

                         # add dimensional reduction information
    if ("tsne_1" in meta_data.columns):
        mat.obsm["X_tsne"] = pd.concat([mat.obs.tsne_1, mat.obs.tsne_2],axis=1).values

    if ("umap_1" in meta_data.columns):
        mat.obsm["X_umap"] = pd.concat([mat.obs.umap_1, mat.obs.umap_2],axis=1).values

    if not raw_data is None:
        mat_ = mat.copy()
        mat_.X = raw_data.transpose().X.copy()
        mat.raw = mat_

    return mat


def seurat_object_to_anndata(file_path_seurat_object, delete_tmp_file=True):
    """
    Convert seurat object into anndata.

    Args:
        file_path_seurat_object (str): File path of seurat object. Seurat object should be saved as Rds format.
        delete_tmp_file (bool): Whether to delete temporary file.
    Returns:
        anndata: anndata object.

    """


    # check file name
    print("input file name: " + file_path_seurat_object)
    if not file_path_seurat_object.lower().endswith(".rds"):
        raise ValueError("Seurat object should be saved as .Rds file")

    # run R script to extract information and make mtx files
    os.makedirs("tmp", exist_ok=True)
    command = f"Rscript {rscript_folder}/seurat_to_mtx.R {file_path_seurat_object}"
    #ret = os.system()
    #if ret == 0:
    #    pass
    #else:
    #    print("Error in R script")

    exec_process(command, message=True, wait_finished=True, return_process=False)

    print("making AnnData ...")

    folder = "./tmp"

    # load data
    mm = sc.read_mtx(os.path.join(folder, "data.mtx"))
    meta = pd.read_csv(os.path.join(folder, "meta_data.csv"),index_col=0)
    meta_dtype = pd.read_csv(os.path.join(folder, "meta_data_dtype.csv"),index_col=0)
    categorical_info = meta_dtype[meta_dtype["dtype"] == "factor"].index.values

    cell_ids = pd.read_csv(os.path.join(folder, "cells_index.csv")).x.values
    variable_ids = pd.read_csv(os.path.join(folder, "variables_index.csv")).x.values


    if "raw_data.mtx" in os.listdir(folder):
        raw_data = sc.read_mtx(os.path.join(folder, "raw_data.mtx"))
        mat = _constructAnnData(mm, cell_ids, variable_ids, meta, categorical_info, raw_data)
    else:
        mat = _constructAnnData(mm, cell_ids, variable_ids, meta, categorical_info)

    # add variable gene info
    if "var_genes.csv" in os.listdir(folder):
        variable_genes = pd.read_csv(os.path.join(folder, "var_genes.csv")).x.values
        mat.var["variable_gene"] = mat.var.index.isin(variable_genes)

    # add color data
    color_df = pd.read_csv("tmp/cluster_color_hex.csv", index_col=0)
    mat.uns["seurat_clusters_colors"] = color_df.colors_hex.values

    # delete temporary files
    if delete_tmp_file:
        shutil.rmtree(folder)


    return mat


def main():
    if len(sys.argv) == 3:
        output_file_name = sys.argv[2]
    else:
        output_file_name = sys.argv[1]

    if not output_file_name.endswith(".h5ad"):
         output_file_name += ".h5ad"


    #print(sys.argv[:])
    adata = seurat_object_to_anndata(file_path_seurat_object=sys.argv[1],
                                     delete_tmp_file=True)

    # save
    print("saving AnnData as h5ad file ...")
    adata.write(filename=output_file_name)

    print("finished")
if __name__ == "__main__":
    main()
