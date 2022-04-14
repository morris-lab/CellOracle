# -*- coding: utf-8 -*-


import pandas as pd
import numpy as np

import sys
import os
import shutil
import glob

import scanpy as sc

from ..utility import exec_process

rscript_folder = os.path.abspath(os.path.dirname(__file__))


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
    #os.makedirs("tmp", exist_ok=True)
    command = f"Rscript {rscript_folder}/seurat_to_mtx.R {file_path_seurat_object}"


    exec_process(command, message=True, wait_finished=True, return_process=False)

    print("making AnnData ...")


    folder = "./tmp"
    assay_files = glob.glob(folder + "/assay_*")
    assay_names = np.unique([os.path.basename(i).split("_")[1] for i in assay_files])

    if len(assay_names) > 1:
        print(f"{len(assay_names)} assays found in the seurat object.")
        print(f"Data is exported as multiple files.")
        print("If the seurat object was made by integrating multiple files, please be careful about 'Simpson's effect' in the inferred GRN.")
        print("Go to CellOracle web documentation for detailed information about this issue.")

    results = {}
    for assay_name in assay_names:
        results[assay_name] = make_anndata_from_files(assay_name=assay_name,
                                                      folder=folder)
    # Delete temporary files
    if delete_tmp_file:
        shutil.rmtree(folder)

    return results


def make_anndata_from_files(assay_name, folder="./tmp"):

    # Load gene expression matrix
    mm = sc.read_mtx(folder + f"/assay_{assay_name}_data.mtx")
    cell_ids = pd.read_csv(folder + f"/assay_{assay_name}_cells.csv").x.values
    gene_ids = pd.read_csv(folder + f"/assay_{assay_name}_genes.csv").x.values

    # Load meta data
    meta_data = pd.read_csv(os.path.join(folder, "meta_data.csv"),index_col=0)
    meta_dtype = pd.read_csv(os.path.join(folder, "meta_data_dtype.csv"),index_col=0)
    categorical_in_meta = meta_dtype[meta_dtype["dtype"] == "factor"].index.values
    for i in meta_data.columns:
        if i in categorical_in_meta:
            meta_data[i] = meta_data[i].astype(np.object) # change dtype

    # Make anndata
    adata = sc.AnnData(mm.X,
                       obs=meta_data.reindex(cell_ids),
                       var=pd.DataFrame(index=gene_ids))

    # Add raw count data if available
    if f"assay_{assay_name}_rawdata.mtx" in os.listdir(folder):
        raw_data = sc.read_mtx(folder + f"/assay_{assay_name}_rawdata.mtx")
        adata.layers["raw_count"] = raw_data.X.copy()


    # Add dimensional reduction data into anndata
    dr_files = glob.glob(folder + "/reduction_*")
    for i in dr_files:
        # Load reduction data
        dr = pd.read_csv(i, index_col=0)
        reduction_name = os.path.basename(i)[10:-4]
        if reduction_name == "FItSNE":
            reduction_name = "tsne"
        adata.obsm[f"X_{reduction_name}"] = dr.values

    # Add variable gene info
    if "var_genes.csv" in os.listdir(folder):
        variable_genes = pd.read_csv(os.path.join(folder, "var_genes.csv")).x.values
        adata.var["variable_gene"] = adata.var.index.isin(variable_genes)

    # Add color data
    color_df = pd.read_csv(os.path.join(folder, "cluster_color_hex.csv"), index_col=0)
    for i in adata.obs.columns:
        if all(adata.obs["active_ident"].values.astype("object") == \
               adata.obs[i].values.astype("object")):
            adata.uns[f"{i}_colors"] = color_df.colors_hex.values

    return adata


def main():
    if len(sys.argv) == 3:
        output_file_name = sys.argv[2]
    else:
        output_file_name = sys.argv[1]

    if not output_file_name.endswith(".h5ad"):
         output_file_name += ".h5ad"

    #print(sys.argv[:])
    adata_dictionary = seurat_object_to_anndata(file_path_seurat_object=sys.argv[1],
                                                delete_tmp_file=True)

    # save
    print("saving AnnData as h5ad file ...")
    if len(adata_dictionary.keys()) >= 2:
        for assay_name, adata in adata_dictionary.items():
            adata.write(filename=output_file_name.replace(".h5ad", f"_{assay_name}.h5ad"))
    else:
        for assay_name, adata in adata_dictionary.items():
            adata.write(filename=output_file_name)
    print("finished")
if __name__ == "__main__":
    main()
