# -*- coding: utf-8 -*-

import warnings
from copy import deepcopy
import math
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from tqdm.auto import tqdm

from ..utility.hdf5_processing import dump_hdf5, load_hdf5

from sklearn.neighbors import NearestNeighbors

from .sankey import sankey
from .markov_simulation import _walk
from .oracle_utility import (
    _adata_to_df,
    _adata_to_color_dict,
    _get_clustercolor_from_anndata,
    _linklist2dict,
    _decompose_TFdict,
    _is_perturb_condition_valid,
    _calculate_relative_ratio_of_simulated_value,
    _check_color_information_and_create_if_not_found,
)
from .oracle_GRN import (
    _do_simulation,
    _getCoefMatrix,
    _coef_to_active_gene_list,
    _shuffle_celloracle_GRN_coef_table,
)
from .modified_VelocytoLoom_class import modified_VelocytoLoom
from ..network_analysis.network_construction import get_links
from ..visualizations.oracle_object_visualization import Oracle_visualization
from ..version import __version__

CONFIG = {"N_PROP_MIN": 1, "N_PROP_MAX": 5, "OOD_WARNING_EXCEEDING_PERCENTAGE": 50}


def update_adata(adata):
    # Update Anndata
    # Anndata generated with Scanpy 1.4 or less should be updated with this function
    # This function will be depricated in the future.

    try:
        lo = adata.uns["draw_graph"]["params"]["layout"]
        if isinstance(lo, np.ndarray):
            lo = lo[0]
        adata.uns["draw_graph"]["params"]["layout"] = lo
    except KeyError:
        pass


def load_oracle(file_path):
    """
    Load oracle object saved as hdf5 file.

    Args:
        file_path (str): File path to the hdf5 file.


    """

    if os.path.exists(file_path):
        pass
    else:
        raise ValueError("File not found. Please check if the file_path is correct.")

    obj = load_hdf5(
        filename=file_path,
        obj_class=Oracle,
        ignore_attrs_if_err=["knn", "knn_smoothing_w", "pca"],
    )

    # Update Anndata
    update_adata(obj.adata)

    return obj


class Oracle(modified_VelocytoLoom, Oracle_visualization):
    """
    Oracle is the main class in CellOracle. Oracle object imports scRNA-seq data (anndata) and TF
    information to infer cluster-specific GRNs. It can predict the future gene expression patterns
    and cell state transitions in response to  the perturbation of TFs. Please see the CellOracle
    paper for details.

    The code of the Oracle class was made of the three components below.

    (1) Anndata: Gene expression matrix and metadata from single-cell RNA-seq are stored
        in the anndata object. Processed values, such as normalized counts and simulated
        values, are stored as layers of anndata. Metadata (i.e., Cluster info) are saved
        in anndata.obs. Refer to scanpy/anndata documentation for detail.

    (2) Net: Net is a custom class in celloracle. Net object processes several data to
        infer GRN. See the Net class documentation for details.

    (3) VelycytoLoom: Calculation of transition probability and visualization of directed
        trajectory graph will be performed in the same way as velocytoloom. VelocytoLoom
        is class from Velocyto, a python library for RNA-velocity analysis. In celloracle,
        we use some functions in velocytoloom for the visualization.


    Attributes:
        adata (anndata): Imported anndata object
        cluster_column_name (str): The column name in adata.obs containing cluster info
        embedding_name (str): The key name in adata.obsm containing dimensional reduction cordinates

    """

    def __init__(self):
        self.celloracle_version_used = __version__
        self.adata = None
        self.cluster_column_name = None
        self.embedding_name = None
        self.ixs_mcmc = None
        self.cluster_specific_TFdict = None
        self.cv_mean_selected_genes = None
        self.TFdict = {}

    ############################
    ### 0. utility functions ###
    ############################
    def copy(self):
        """
        Deepcopy itself.
        """
        return deepcopy(self)

    def to_hdf5(self, file_path):
        """
        Save object as hdf5.

        Args:
            file_path (str): file path to save file. Filename needs to end with '.celloracle.oracle'
        """
        if file_path.endswith(".celloracle.oracle"):
            pass
        else:
            raise ValueError("Filename needs to end with '.celloracle.oracle'")

        compression_opts = 5
        dump_hdf5(
            obj=self,
            filename=file_path,
            data_compression=compression_opts,
            chunks=(2048, 2048),
            noarray_compression=compression_opts,
            pickle_protocol=4,
        )

    def _generate_meta_data(self):
        info = {}
        if hasattr(self, "celloracle_version_used"):
            info[
                "celloracle version used for instantiation"
            ] = self.celloracle_version_used
        else:
            info["celloracle version used for instantiation"] = "NA"

        if self.adata is not None:
            info["n_cells"] = self.adata.shape[0]
            info["n_genes"] = self.adata.shape[1]
            info["status - Gene expression matrix"] = "Ready"
        else:
            info["n_cells"] = "NA"
            info["n_genes"] = "NA"
            info["status - Gene expression matrix"] = "Not imported"

        info["cluster_name"] = self.cluster_column_name
        info["dimensional_reduction_name"] = self.embedding_name

        if hasattr(self, "all_target_genes_in_TFdict"):
            pass
        else:
            if self.adata is not None:
                self._process_TFdict_metadata(verbose=False)
        if self.adata is not None:
            info[
                "n_target_genes_in_TFdict"
            ] = f"{len(self.all_target_genes_in_TFdict)} genes"
            info[
                "n_regulatory_in_TFdict"
            ] = f"{len(self.all_regulatory_genes_in_TFdict)} genes"
            info[
                "n_regulatory_in_both_TFdict_and_scRNA-seq"
            ] = f"{self.adata.var['isin_TFdict_regulators'].sum()} genes"
            info[
                "n_target_genes_both_TFdict_and_scRNA-seq"
            ] = f"{self.adata.var['isin_TFdict_targets'].sum()} genes"

        if len(self.TFdict.keys()) > 0:
            info["status - BaseGRN"] = "Ready"
        else:
            info["status - BaseGRN"] = "Not imported"

        if hasattr(self, "pcs"):
            info["status - PCA calculation"] = "Done"
        else:
            info["status - PCA calculation"] = "Not finished"

        if hasattr(self, "knn"):
            info["status - Knn imputation"] = "Done"
        else:
            info["status - Knn imputation"] = "Not finished"

        if hasattr(self, "k_knn_imputation"):
            info["k_for_knn_imputation"] = self.k_knn_imputation
        else:
            info["k_for_knn_imputation"] = "NA"

        if hasattr(self, "coef_matrix_per_cluster") | hasattr(self, "coef_matrix"):
            info["status - GRN calculation for simulation"] = "Done"
        else:
            info["status - GRN calculation for simulation"] = "Not finished"
        return info

    def __repr__(self):
        info = self._generate_meta_data()

        message = "Oracle object\n\n"
        message += "Meta data\n"
        message_status = "Status\n"
        for key, value in info.items():
            if key.startswith("status"):
                message_status += (
                    "    " + key.replace("status - ", "") + ": " + str(value) + "\n"
                )
            else:
                message += "    " + key + ": " + str(value) + "\n"

        message += message_status

        return message

    ###################################
    ### 1. Methods for loading data ###
    ###################################
    def _process_TFdict_metadata(self, verbose=True):
        # Make list of all target genes and all reguolatory genes in the TFdict
        (
            self.all_target_genes_in_TFdict,
            self.all_regulatory_genes_in_TFdict,
        ) = _decompose_TFdict(TFdict=self.TFdict)

        # Intersect gene between the list above and gene expression matrix.
        self.adata.var["symbol"] = self.adata.var.index.values
        self.adata.var["isin_TFdict_targets"] = self.adata.var.symbol.isin(
            self.all_target_genes_in_TFdict
        )
        self.adata.var["isin_TFdict_regulators"] = self.adata.var.symbol.isin(
            self.all_regulatory_genes_in_TFdict
        )

        # n_target = self.adata.var["isin_TFdict_targets"].sum()
        # if n_target == 0:
        # print("Found no overlap between TF info (base GRN) and your scRNA-seq data. Please check your data format and species.")
        if verbose:
            n_reg = self.adata.var["isin_TFdict_regulators"].sum()
            if n_reg == 0:
                print(
                    "Found No overlap between TF info (base GRN) and your scRNA-seq data. Please check your data format and species."
                )

            elif n_reg < 50:
                print(
                    f"Total number of TF was {n_reg}. Although we can go to the GRN calculation with this data, but the TF number is small."
                )

    def import_TF_data(
        self, TF_info_matrix=None, TFdict=None, TF_info_matrix_path=None
    ):
        """
        Load data about potential-regulatory TFs.
        You can import either TF_info_matrix or TFdict.
        For more information on how to make these files, please see the motif analysis module within the celloracle tutorial.

        Args:
            TF_info_matrix (pandas.DataFrame, str): TF_info_matrix or path to a TF_info_matrix.

            TFdict (dictionary): Python dictionary of TF info.
        """
        if TF_info_matrix_path is not None:
            import warnings

            warnings.warn(
                "DeprecationWarning: In the future, please pass path to TF matrix to `TF_info_matrix`"
            )
            TF_info_matrix = TF_info_matrix_path

        if self.adata is None:
            raise ValueError("Please import scRNA-seq data first.")

        if len(self.TFdict) != 0:
            print("TF dict already exists. The old TF dict data was deleted. \n")

        if isinstance(TF_info_matrix, pd.DataFrame):
            tmp = TF_info_matrix.copy()
            tmp = tmp.drop(["peak_id"], axis=1)
            tmp = tmp.groupby(by="gene_short_name").sum()
            self.TFdict = dict(tmp.apply(lambda x: x[x > 0].index.values, axis=1))
        elif type(TF_info_matrix_path) == str:
            tmp = pd.read_parquet(TF_info_matrix_path)
            tmp = tmp.drop(["peak_id"], axis=1)
            tmp = tmp.groupby(by="gene_short_name").sum()
            self.TFdict = dict(tmp.apply(lambda x: x[x > 0].index.values, axis=1))
        else:
            raise ValueError(
                "`TF_info_matrix` must be a pandas.DataFrame or a file path"
            )

        if not TFdict is None:
            self.TFdict = TFdict.copy()

        # Update summary of TFdata
        self._process_TFdict_metadata()

    def updateTFinfo_dictionary(self, TFdict={}):
        """
        Update a TF dictionary.
        If a key in the new TF dictionary already exists in the old TF dictionary, old values will be replaced with a new one.

        Args:
            TFdict (dictionary): Python dictionary of TF info.
        """

        self.TFdict.update(TFdict)

        # Update summary of TFdata
        self._process_TFdict_metadata()

    def addTFinfo_dictionary(self, TFdict):
        """
        Add new TF info to pre-existing TFdict.
        Values in the old TF dictionary will remain.

        Args:
            TFdict (dictionary): Python dictionary of TF info.
        """
        for tf in TFdict:
            if tf in self.TFdict.keys():
                targets = self.TFdict[tf]
                targets = list(TFdict[tf]) + list(targets)
                targets = np.unique(targets)
                self.TFdict.update({tf: targets})
            else:
                self.TFdict.update({tf: TFdict[tf]})

        # Update summary of TFdata
        self._process_TFdict_metadata()

    def get_cluster_specific_TFdict_from_Links(
        self, links_object, ignore_warning=False
    ):
        """
        Extract TF and its target gene information from Links object.
        This function can be used to reconstruct GRNs based on pre-existing GRNs saved in Links object.

        Args:
            links_object (Links): Please see the explanation of Links class.

        """
        # Check cluster unit in oracle object is same as cluster unit in links_object
        clusters_in_oracle_object = sorted(
            self.adata.obs[self.cluster_column_name].unique()
        )
        clusters_in_link_object = sorted(links_object.cluster)
        if (self.cluster_column_name == links_object.name) & (
            clusters_in_link_object == clusters_in_oracle_object
        ):
            pass
        else:
            if ignore_warning:
                pass
            else:
                raise ValueError(
                    "Clustering unit does not match. Please prepare links object that was made with same cluster data."
                )

        self.cluster_specific_TFdict = {}

        for i in links_object.filtered_links:
            self.cluster_specific_TFdict[i] = _linklist2dict(
                links_object.filtered_links[i]
            )

    def import_anndata_as_raw_count(
        self,
        adata,
        cluster_column_name=None,
        embedding_name=None,
        transform="natural_log",
    ):
        """
        Load scRNA-seq data. scRNA-seq data should be prepared as an anndata object.
        Preprocessing (cell and gene filtering, dimensional reduction, clustering, etc.) should be done before loading data.
        The method imports RAW GENE COUNTS because unscaled and uncentered gene expression data are required for the GRN inference and simulation.
        See tutorial notebook for the details about how to process scRNA-seq data.

        Args:
            adata (anndata): anndata object that stores scRNA-seq data.

            cluster_column_name (str): the name of column containing cluster information in anndata.obs.
                Clustering data should be in anndata.obs.

            embedding_name (str): the key name for dimensional reduction information in anndata.obsm.
                Dimensional reduction (or 2D trajectory graph) should be in anndata.obsm.

            transform (str): The method for log-transformation. Chose one from "natural_log" or "log2".

        """
        if adata.X.min() < 0:
            raise ValueError(
                "gene expression matrix (adata.X) does not seems to be raw_count because it contains negavive values."
            )

        if (adata.shape[1] < 1000) | (adata.shape[1] > 4000):
            print(
                f"{adata.shape[1]} genes were found in the adata. Note that Celloracle is intended to use around 1000-3000 genes, so the behavior with this number of genes may differ from what is expected."
            )

        # store data
        self.adata = adata.copy()

        # update anndata format
        update_adata(self.adata)

        self.cluster_column_name = cluster_column_name
        self.embedding_name = embedding_name
        self.embedding = self.adata.obsm[embedding_name].copy()

        # if hasattr(self.adata, "raw"):
        #    self.adata.X = self.adata.raw.X.copy()

        # store raw count data
        self.adata.layers["raw_count"] = self.adata.X.copy()

        # log transformation
        if transform == "log2":
            self.adata.X = np.log2(self.adata.X + 1)
        elif transform == "natural_log":
            sc.pp.log1p(self.adata)

        self.adata.layers["normalized_count"] = self.adata.X.copy()

        # update color information
        _check_color_information_and_create_if_not_found(
            adata=self.adata,
            cluster_column_name=cluster_column_name,
            embedding_name=embedding_name,
        )
        col_dict = _get_clustercolor_from_anndata(
            adata=self.adata, cluster_name=self.cluster_column_name, return_as="dict"
        )
        self.colorandum = np.array(
            [col_dict[i] for i in self.adata.obs[self.cluster_column_name]]
        )

        # variable gene detection for the QC of simulation
        n = min(adata.shape[1], 1000) - 1

        self.score_cv_vs_mean(n, plot=False, max_expr_avg=35)
        self.high_var_genes = self.cv_mean_selected_genes.copy()
        self.cv_mean_selected_genes = None

        self.adata.var["symbol"] = self.adata.var.index.values
        self.adata.var["isin_top1000_var_mean_genes"] = self.adata.var.symbol.isin(
            self.high_var_genes
        )

    def import_anndata_as_normalized_count(
        self, adata, cluster_column_name=None, embedding_name=None, test_mode=False
    ):
        """
        Load scRNA-seq data. scRNA-seq data should be prepared as an anndata object.
        Preprocessing (cell and gene filtering, dimensional reduction, clustering, etc.) should be done before loading data.
        The method will import NORMALIZED and LOG TRANSFORMED data but NOT SCALED and NOT CENTERED data.
        See the tutorial for more details on how to process scRNA-seq data.

        Args:
            adata (anndata): anndata object containing scRNA-seq data.

            cluster_column_name (str): the name of column containing cluster information in anndata.obs.
                Clustering data should be in anndata.obs.

            embedding_name (str): the key name for dimensional reduction information in anndata.obsm.
                Dimensional reduction (or 2D trajectory graph) should be in anndata.obsm.

            transform (str): The method for log-transformation. Chose one from "natural_log" or "log2".
        """
        if adata.X.min() < 0:
            raise ValueError(
                "Gene expression matrix (adata.X) contains negavive values. Please use UNSCALED and UNCENTERED data."
            )

        if (adata.shape[1] < 1000) | (adata.shape[1] > 4000):
            print(
                f"{adata.shape[1]} genes were found in the adata. Note that Celloracle is intended to use around 1000-3000 genes, so the behavior with this number of genes may differ from what is expected."
            )

        # Store data
        self.adata = adata.copy()

        # Update anndata format
        update_adata(self.adata)

        self.cluster_column_name = cluster_column_name
        self.embedding_name = embedding_name
        self.embedding = self.adata.obsm[embedding_name].copy()

        # store raw count data
        # self.adata.layers["raw_count"] = adata.X.copy()

        # normalization and log transformation
        self.adata.layers["normalized_count"] = self.adata.X.copy()

        # update color information
        if not test_mode:
            _check_color_information_and_create_if_not_found(
                adata=self.adata,
                cluster_column_name=cluster_column_name,
                embedding_name=embedding_name,
            )
            col_dict = _get_clustercolor_from_anndata(
                adata=self.adata,
                cluster_name=self.cluster_column_name,
                return_as="dict",
            )
            self.colorandum = np.array(
                [col_dict[i] for i in self.adata.obs[self.cluster_column_name]]
            )

            # variable gene detection for the QC of simulation
            n = min(adata.shape[1], 1000) - 1

            self.score_cv_vs_mean(n, plot=False, max_expr_avg=35)
            self.high_var_genes = self.cv_mean_selected_genes.copy()
            self.cv_mean_selected_genes = None

            self.adata.var["symbol"] = self.adata.var.index.values
            self.adata.var["isin_top1000_var_mean_genes"] = self.adata.var.symbol.isin(
                self.high_var_genes
            )

    def change_cluster_unit(self, new_cluster_column_name):
        """
        Change clustering unit.
        If you change cluster, previous GRN data and simulation data will be delated.
        Please re-calculate GRNs.
        """

        # 1. Check new cluster information exists in anndata.
        if new_cluster_column_name in self.adata.obs.columns:
            _check_color_information_and_create_if_not_found(
                adata=self.adata,
                cluster_column_name=new_cluster_column_name,
                embedding_name=self.embedding_name,
            )
        else:
            raise ValueError(f"{new_cluster_column_name} was not found in anndata")

        # 2. Reset previous GRN data and simoulation data
        attributes_delete = [
            "ixs_mcmc",
            "colorandum",
            "alpha_for_trajectory_GRN",
            "GRN_unit",
            "coef_matrix_per_cluster",
            "perturb_condition",
            "corr_calc",
            "embedding_knn",
            "sampling_ixs",
            "corrcoef",
            "corrcoef_random",
            "transition_prob",
            "transition_prob_random",
            "delta_embedding",
            "delta_embedding_random",
            "total_p_mass",
            "flow_embedding",
            "flow_grid",
            "flow",
            "flow_norm",
            "flow_norm_magnitude",
            "flow_rndm",
            "flow_norm_rndm",
            "flow_norm_magnitude_rndm",
        ]

        attributes = list(self.__dict__.keys())
        for i in attributes:
            if i in attributes_delete:
                delattr(self, i)

        # 4. Update cluster info
        self.cluster_column_name = new_cluster_column_name

        # 3. Update color information
        col_dict = _get_clustercolor_from_anndata(
            adata=self.adata, cluster_name=new_cluster_column_name, return_as="dict"
        )
        self.colorandum = np.array(
            [col_dict[i] for i in self.adata.obs[new_cluster_column_name]]
        )

    def update_cluster_colors(self, palette):
        """
        Update color information stored in the oracle object.
        The color information is overwritten.
        """

        sc.pl.embedding(
            self.adata,
            basis=self.embedding_name,
            color=self.cluster_column_name,
            palette=palette,
        )

        col_dict = _get_clustercolor_from_anndata(
            adata=self.adata, cluster_name=self.cluster_column_name, return_as="dict"
        )
        self.colorandum = np.array(
            [col_dict[i] for i in self.adata.obs[self.cluster_column_name]]
        )

    ####################################
    ### 2. Methods for GRN inference ###
    ####################################
    def fit_GRN_for_simulation(
        self,
        GRN_unit="cluster",
        alpha=1,
        use_cluster_specific_TFdict=False,
        verbose_level=1,
    ):
        """
        Do GRN inference.
        Please see the paper of CellOracle paper for details.

        GRN can be constructed for the entire population or each clusters.
        If you want to infer cluster-specific GRN, please set [GRN_unit="cluster"].
        You can select cluster information when you import data.

        If you set [GRN_unit="whole"], GRN will be made using all cells.

        Args:
            GRN_unit (str): Select "cluster" or "whole"

            alpha (float or int): The strength of regularization.
                If you set a lower value, the sensitivity increases, and you can detect weaker network connections. However, there may be more noise.
                If you select a higher value, it will reduce the chance of overfitting.

            verbose_level (int): if [verbose_level>1], most detailed progress information will be shown.
                if [1 >= verbose_level > 0], one progress bar will be shown.
                if [verbose_level == 0], no progress bar will be shown.

        """

        if verbose_level > 1:
            verbose_cluster = True
            verbose_gene = True
        elif 0 < verbose_level <= 1:
            verbose_cluster = True
            verbose_gene = False
        else:
            verbose_cluster = False
            verbose_gene = False

        # prepare data for GRN calculation
        gem_imputed = _adata_to_df(self.adata, "imputed_count")
        self.adata.layers["simulation_input"] = self.adata.layers[
            "imputed_count"
        ].copy()
        self.alpha_for_trajectory_GRN = alpha
        self.GRN_unit = GRN_unit

        if use_cluster_specific_TFdict & (self.cluster_specific_TFdict is not None):
            self.coef_matrix_per_cluster = {}
            cluster_info = self.adata.obs[self.cluster_column_name]
            with tqdm(
                np.unique(cluster_info), disable=(verbose_cluster == False)
            ) as pbar:
                for cluster in pbar:
                    pbar.set_postfix(cluster=f"{cluster}")
                    cells_in_the_cluster_bool = cluster_info == cluster
                    gem_ = gem_imputed[cells_in_the_cluster_bool]
                    self.coef_matrix_per_cluster[cluster] = _getCoefMatrix(
                        gem=gem_,
                        TFdict=self.cluster_specific_TFdict[cluster],
                        alpha=alpha,
                        verbose=verbose_gene,
                    )
        else:
            if GRN_unit == "whole":
                self.coef_matrix = _getCoefMatrix(
                    gem=gem_imputed,
                    TFdict=self.TFdict,
                    alpha=alpha,
                    verbose=verbose_gene,
                )
            if GRN_unit == "cluster":
                self.coef_matrix_per_cluster = {}
                cluster_info = self.adata.obs[self.cluster_column_name]
                with tqdm(
                    np.unique(cluster_info), disable=(verbose_cluster == False)
                ) as pbar:
                    for cluster in pbar:
                        pbar.set_postfix(cluster=f"{cluster}")
                        cells_in_the_cluster_bool = cluster_info == cluster
                        gem_ = gem_imputed[cells_in_the_cluster_bool]
                        self.coef_matrix_per_cluster[cluster] = _getCoefMatrix(
                            gem=gem_,
                            TFdict=self.TFdict,
                            alpha=alpha,
                            verbose=verbose_gene,
                        )

        self.extract_active_gene_lists(verbose=False)

    def extract_active_gene_lists(self, return_as=None, verbose=False):
        """
        Args:
            return_as (str): If not None, it returns dictionary or list. Chose either "indivisual_dict" or "unified_list".
            verbose (bool): Whether to show progress bar.

        Returns:
            dictionary or list: The format depends on the argument, "return_as".

        """
        if return_as not in ["indivisual_dict", "unified_list", None]:
            raise ValueError(
                "return_as should be either 'indivisual_dict' or 'unified_list'."
            )

        if not hasattr(self, "GRN_unit"):
            try:
                loop = self.coef_matrix_per_cluster.items()
                self.GRN_unit = "cluster"
                print("Currently selected GRN_unit: ", self.GRN_unit)

            except:
                try:
                    loop = {"whole_cell": self.coef_matrix}.items()
                    self.GRN_unit = "whole"
                    print("Currently selected GRN_unit: ", self.GRN_unit)
                except:
                    raise ValueError(
                        "GRN is not ready. Please run 'fit_GRN_for_simulation' first."
                    )

        elif self.GRN_unit == "cluster":
            loop = self.coef_matrix_per_cluster.items()
        elif self.GRN_unit == "whole":
            loop = {"whole_cell": self.coef_matrix}.items()

        if verbose:
            loop = tqdm(loop)

        unified_list = []
        indivisual_dict = {}
        for cluster, coef_matrix in loop:
            active_genes = _coef_to_active_gene_list(coef_matrix)
            unified_list += active_genes
            indivisual_dict[cluster] = active_genes

        unified_list = list(np.unique(unified_list))

        # Store data
        self.active_regulatory_genes = unified_list.copy()
        self.adata.var["symbol"] = self.adata.var.index.values
        if "isin_top1000_var_mean_genes" not in self.adata.var.columns:
            self.adata.var["isin_top1000_var_mean_genes"] = self.adata.var.symbol.isin(
                self.high_var_genes
            )
        self.adata.var["isin_actve_regulators"] = self.adata.var.symbol.isin(
            unified_list
        )

        if return_as == "indivisual_dict":
            return indivisual_dict

        elif return_as == "unified_list":
            return unified_list

    #######################################################
    ### 3. Methods for simulation of signal propagation ###
    #######################################################

    def simulate_shift(
        self,
        perturb_condition=None,
        GRN_unit=None,
        n_propagation=3,
        ignore_warning=False,
        use_randomized_GRN=False,
        clip_delta_X=False,
    ):
        """
        Simulate signal propagation with GRNs. Please see the CellOracle paper for details.
        This function simulates a gene expression pattern in the near future.
        Simulated values will be stored in anndata.layers: ["simulated_count"]


        The simulation use three types of data.
        (1) GRN inference results (coef_matrix).
        (2) Perturb_condition: You can set arbitrary perturbation condition.
        (3) Gene expression matrix: The simulation starts from imputed gene expression data.

        Args:
            perturb_condition (dictionary): condition for perturbation.
               if you want to simulate knockout for GeneX, please set [perturb_condition={"GeneX": 0.0}]
               Although you can set any non-negative values for the gene condition, avoid setting biologically infeasible values for the perturb condition.
               It is strongly recommended to check gene expression values in your data before selecting the perturb condition.

            GRN_unit (str): GRN type. Please select either "whole" or "cluster". See the documentation of "fit_GRN_for_simulation" for the detailed explanation.

            n_propagation (int): Calculation will be performed iteratively to simulate signal propagation in GRN.
                You can set the number of steps for this calculation.
                With a higher number, the results may recapitulate signal propagation for many genes.
                However, a higher number of propagation may cause more error/noise.

            clip_delta_X (bool): If simulated gene expression shift can lead to gene expression value that is outside of WT distribution, such gene expression is clipped to WT range.
        """
        self.__simulate_shift(
            perturb_condition=perturb_condition,
            GRN_unit=GRN_unit,
            n_propagation=n_propagation,
            ignore_warning=ignore_warning,
            use_randomized_GRN=use_randomized_GRN,
            clip_delta_X=clip_delta_X,
        )

    def __simulate_shift(
        self,
        perturb_condition=None,
        GRN_unit=None,
        n_propagation=3,
        ignore_warning=False,
        use_randomized_GRN=False,
        n_min=None,
        n_max=None,
        clip_delta_X=False,
    ):
        """
        Simulate signal propagation with GRNs. Please see the CellOracle paper for details.
        This function simulates a gene expression pattern in the near future.
        Simulated values will be stored in anndata.layers: ["simulated_count"]


        The simulation use three types of data.
        (1) GRN inference results (coef_matrix).
        (2) Perturb_condition: You can set arbitrary perturbation condition.
        (3) Gene expression matrix: The simulation starts from imputed gene expression data.

        Args:
            perturb_condition (dictionary): condition for perturbation.
               if you want to simulate knockout for GeneX, please set [perturb_condition={"GeneX": 0.0}]
               Although you can set any non-negative values for the gene condition, avoid setting biologically infeasible values for the perturb condition.
               It is strongly recommended to check gene expression values in your data before selecting the perturb condition.

            GRN_unit (str): GRN type. Please select either "whole" or "cluster". See the documentation of "fit_GRN_for_simulation" for the detailed explanation.

            n_propagation (int): Calculation will be performed iteratively to simulate signal propagation in GRN.
                You can set the number of steps for this calculation.
                With a higher number, the results may recapitulate signal propagation for many genes.
                However, a higher number of propagation may cause more error/noise.
        """

        # 0. Reset previous simulation results if it exist
        # self.ixs_mcmc = None
        # self.mcmc_transition_id = None
        # self.corrcoef = None
        # self.transition_prob = None
        # self.tr = None
        if n_min is None:
            n_min = CONFIG["N_PROP_MIN"]
        if n_max is None:
            n_max = CONFIG["N_PROP_MAX"]
        self._clear_simulation_results()

        if GRN_unit is not None:
            self.GRN_unit = GRN_unit
        elif hasattr(self, "GRN_unit"):
            GRN_unit = self.GRN_unit
            # print("Currently selected GRN_unit: ", self.GRN_unit)
        elif hasattr(self, "coef_matrix_per_cluster"):
            GRN_unit = "cluster"
            self.GRN_unit = GRN_unit
        elif hasattr(self, "coef_matrix"):
            GRN_unit = "whole"
            self.GRN_unit = GRN_unit
        else:
            raise ValueError(
                "GRN is not ready. Please run 'fit_GRN_for_simulation' first."
            )

        if use_randomized_GRN:
            print("Attention: Using randomized GRN for the perturbation simulation.")

        # 1. prepare perturb information

        self.perturb_condition = perturb_condition.copy()

        # Prepare metadata before simulation
        if not hasattr(self, "active_regulatory_genes"):
            self.extract_active_gene_lists(verbose=False)

        if not hasattr(self, "all_regulatory_genes_in_TFdict"):
            self._process_TFdict_metadata()

        for i, value in perturb_condition.items():
            # 1st Sanity check
            if not i in self.adata.var.index:
                raise ValueError(
                    f"Gene {i} is not included in the Gene expression matrix."
                )

            # 2nd Sanity check
            if i not in self.all_regulatory_genes_in_TFdict:
                raise ValueError(
                    f"Gene {i} is not included in the base GRN; It is not TF or TF motif information is not available. Cannot perform simulation."
                )

            # 3rd Sanity check
            if i not in self.active_regulatory_genes:
                raise ValueError(
                    f"Gene {i} does not have enough regulatory connection in the GRNs. Cannot perform simulation."
                )

            # 4th Sanity check
            if i not in self.high_var_genes:
                if ignore_warning:
                    pass
                    # print(f"Variability score of Gene {i} is too low. Simulation accuracy may be poor with this gene.")
                else:
                    pass
                    # print(f"Variability score of Gene {i} is too low. Simulation accuracy may be poor with this gene.")
                    # raise ValueError(f"Variability score of Gene {i} is too low. Cannot perform simulation.")

            # 5th Sanity check
            if value < 0:
                raise ValueError(f"Negative gene expression value is not allowed.")

            # 6th Sanity check
            safe = _is_perturb_condition_valid(
                adata=self.adata, goi=i, value=value, safe_range_fold=2
            )
            if not safe:
                if ignore_warning:
                    pass
                else:
                    raise ValueError(
                        f"Input perturbation condition is far from actural gene expression value. Please follow the recommended usage. "
                    )
            # 7th QC
            if n_min <= n_propagation <= n_max:
                pass
            else:
                raise ValueError(
                    f"n_propagation value error. It should be an integer from {n_min} to {n_max}."
                )

        # reset simulation initiation point
        self.adata.layers["simulation_input"] = self.adata.layers[
            "imputed_count"
        ].copy()
        simulation_input = _adata_to_df(self.adata, "simulation_input")
        for i in perturb_condition.keys():
            simulation_input[i] = perturb_condition[i]

        # 2. load gene expression matrix (initiation information for the simulation)
        gem_imputed = _adata_to_df(self.adata, "imputed_count")

        # 3. do simulation for signal propagation within GRNs
        if GRN_unit == "whole":
            if use_randomized_GRN == False:
                coef_matrix = self.coef_matrix.copy()
            else:
                if hasattr(self, "coef_matrix_randomized") == False:
                    print("The random coef matrix was calculated.")
                    self.calculate_randomized_coef_table()
                coef_matrix = self.coef_matrix_randomized.copy()
            gem_simulated = _do_simulation(
                coef_matrix=coef_matrix,
                simulation_input=simulation_input,
                gem=gem_imputed,
                n_propagation=n_propagation,
            )

        elif GRN_unit == "cluster":
            simulated = []
            cluster_info = self.adata.obs[self.cluster_column_name]
            for cluster in np.unique(cluster_info):
                if use_randomized_GRN == False:
                    coef_matrix = self.coef_matrix_per_cluster[cluster].copy()
                else:
                    if hasattr(self, "coef_matrix_per_cluster_randomized") == False:
                        print("The random coef matrix was calculated.")
                        self.calculate_randomized_coef_table()
                    coef_matrix = self.coef_matrix_per_cluster_randomized[
                        cluster
                    ].copy()
                cells_in_the_cluster_bool = cluster_info == cluster
                simulation_input_ = simulation_input[cells_in_the_cluster_bool]
                gem_ = gem_imputed[cells_in_the_cluster_bool]

                simulated_in_the_cluster = _do_simulation(
                    coef_matrix=coef_matrix,
                    simulation_input=simulation_input_,
                    gem=gem_,
                    n_propagation=n_propagation,
                )

                simulated.append(simulated_in_the_cluster)
            gem_simulated = pd.concat(simulated, axis=0)
            gem_simulated = gem_simulated.reindex(gem_imputed.index)

        else:
            raise ValueError("GRN_unit shold be either of 'whole' or 'cluster'")

        # 4. store simulation results
        #  simulated future gene expression matrix
        self.adata.layers["simulated_count"] = gem_simulated.values

        #  difference between simulated values and original values
        self.adata.layers["delta_X"] = (
            self.adata.layers["simulated_count"] - self.adata.layers["imputed_count"]
        )

        # Clip simulated gene expression to avoid out of distribution prediction.
        if clip_delta_X:
            self.clip_delta_X()

        # Sanity check; check distribution of simulated values. If the value is far from original gene expression range, it will give warning.
        if ignore_warning:
            pass
        else:
            ood_stat = self.evaluate_simulated_gene_distribution_range()
            ood_stat = ood_stat[
                ood_stat.Max_exceeding_ratio
                > CONFIG["OOD_WARNING_EXCEEDING_PERCENTAGE"] / 100
            ]
            if len(ood_stat) > 0:
                message = f"There may be out of distribution prediction in {len(ood_stat)} genes. It is recommended to set `clip_delta_X=True` to avoid the out of distribution prediction."
                message += "\n To see the detail, please run `oracle.evaluate_simulated_gene_distribution_range()`"
                warnings.warn(message, UserWarning, stacklevel=2)

    def _clear_simulation_results(self):
        att_list = [
            "flow_embedding",
            "flow_grid",
            "flow",
            "flow_norm_magnitude",
            "flow_rndm",
            "flow_norm_rndm",
            "flow_norm_magnitude_rndm",
            "corrcoef",
            "corrcoef_random",
            "transition_prob",
            "transition_prob_random",
            "delta_embedding",
            "delta_embedding_random",
            "ixs_mcmc",
            "mcmc_transition_id",
            "tr",
        ]

        for i in att_list:
            if hasattr(self, i):
                setattr(self, i, None)

    def evaluate_simulated_gene_distribution_range(self):
        """
        CellOracle does not intend to simulate out-of-distribution simulation.
        This function evaluates how the simulated gene expression values differ from the undisturbed gene expression distribution range.
        """

        exceedance = self._calculate_potential_OOD_exceedance_ratio()
        statistics = pd.concat(
            [exceedance.max(), (exceedance != 0).mean(axis=0)], axis=1
        )
        statistics.columns = ["Max_exceeding_ratio", "OOD_cell_ratio"]

        statistics = statistics.sort_values(by="Max_exceeding_ratio", ascending=False)

        return statistics

    def _calculate_potential_OOD_exceedance_ratio(self):
        """
        CellOracle does not intend to simulate out-of-distribution simulation.
        This function evaluates how the simulated gene expression values differ from the undisturbed gene expression distribution range.

        Args:

            pandas.DataFrame: The value is exceeding ratio.
        """
        # Sanity check
        if "simulated_count" in self.adata.layers.keys():
            pass
        else:
            raise ValueError("Simulation results not found. Run simulation first.")

        simulated_count = self.adata.to_df(layer="simulated_count")
        imputed_count = self.adata.to_df(layer="imputed_count")

        relative_ratio = _calculate_relative_ratio_of_simulated_value(
            simulated_count=simulated_count, reference_count=imputed_count
        )

        lower_exceedance = np.clip(relative_ratio, -np.inf, 0).abs()
        higer_exceedance = np.clip(relative_ratio - 1, 0, np.inf)
        exceedance = pd.DataFrame(
            np.maximum(lower_exceedance.values, higer_exceedance.values),
            index=relative_ratio.index,
            columns=relative_ratio.columns,
        )
        return exceedance

    def evaluate_and_plot_simulation_value_distribution(
        self, n_genes=4, n_bins=10, alpha=0.5, figsize=[5, 3], save=None
    ):
        """
        This function will visualize distribution of original gene expression value and simulated values.
        This cunction is built to confirm there is no significant out-of-distribution in the simulation results.

        Args:
            n_genes (int): Number of genes to show. This functin will show the results of top n_genes with large difference between original and simulation values.
            n_bins (int): Number of bins.
            alpha (float): Transparency.
            figsize ([float, float]): Figure size.
            save (str): Folder path to save your plots. If it is not specified, no figure is saved.
        Return:
            None
        """

        simulated_count = self.adata.to_df(layer="simulated_count")
        imputed_count = self.adata.to_df(layer="imputed_count")

        ood_stats = self.evaluate_simulated_gene_distribution_range()

        if save is not None:
            os.makedirs(save, exist_ok=True)

        for goi, val in ood_stats[:4].iterrows():
            fig, ax = plt.subplots(figsize=figsize)
            in_range_cell_ratio = 1 - val["OOD_cell_ratio"]
            ax.hist(
                imputed_count[goi], label="Original value", alpha=alpha, bins=n_bins
            )
            ax.hist(
                simulated_count[goi], label="Simulation value", alpha=alpha, bins=n_bins
            )
            message = f"Gene: ${goi}$, "
            message += f"Cells in original gene range: {in_range_cell_ratio*100:.5g}%, "
            message += f"\nMax exceedance: {val['Max_exceeding_ratio']*100:.3g}%"
            plt.title(message)
            plt.legend()
            plt.xlabel("Gene expression")
            plt.ylabel("Count")
            plt.subplots_adjust(left=0.2, right=0.8, top=0.8, bottom=0.15)
            if save is not None:
                fig.savefig(
                    os.path.join(save, f"gene_expression_histogram_Spi1_KO_{goi}.png"),
                    transparent=True,
                )
            plt.show()

    def clip_delta_X(self):
        """
        To avoid issue caused by out-of-distribution prediction, this function clip simulated gene expression value to the unperturbed gene expression range.
        """
        # Sanity check
        if "simulated_count" in self.adata.layers.keys():
            pass
        else:
            raise ValueError("Simulation results not found. Run simulation first.")

        simulated_count = self.adata.to_df(layer="simulated_count").copy()
        imputed_count = self.adata.to_df(layer="imputed_count").copy()

        min_ = imputed_count.min(axis=0)
        max_ = imputed_count.max(axis=0)

        for goi in simulated_count.columns:
            simulated_count[goi] = np.clip(simulated_count[goi], min_[goi], max_[goi])

        self.adata.layers["simulated_count"] = simulated_count.values
        self.adata.layers["delta_X"] = (
            self.adata.layers["simulated_count"] - self.adata.layers["imputed_count"]
        )

    def estimate_impact_of_perturbations_under_various_ns(
        self, perturb_condition, order=1, n_prop_max=5, GRN_unit=None, figsize=[7, 3]
    ):
        """
        This function is designed to help user to estimate appropriate n for signal propagation.
        The function will do the following calculation for each n and plot results as line plot.
        1. Calculate signal propagation.
        2. Calculate the vector length of delta_X, which represents the simulated shift vector for each cell in the multi dimensional gene expression space.
        3. Calculate mean of delta_X for each cluster.
        Repeat step 1~3 for each n and plot results as a line plot.

        Args:
            perturb_condition (dictionary): Please refer to the function 'simulate_shift' for detail of this.
            order (int): If order=1, this function calculate l1 norm. If order=2, it calculate l2 norm.
            n_prop_max (int): Max of n to try.
        Return
            matplotlib figure
        """
        lengths = []
        for i in tqdm(range(0, n_prop_max + 1)):
            self.__simulate_shift(
                perturb_condition=perturb_condition,
                GRN_unit=None,
                n_propagation=i,
                ignore_warning=False,
                use_randomized_GRN=False,
                n_min=0,
                n_max=n_prop_max + 1,
            )

            delta = self.adata.to_df(layer="delta_X")
            length = np.linalg.norm(delta, ord=order, axis=1)
            lengths.append(length)

        lengths = pd.DataFrame(lengths).transpose()
        lengths.columns = [f"{i}" for i in range(0, n_prop_max + 1)]
        lengths["group"] = self.adata.obs[self.cluster_column_name].values

        # Plot results
        fig, ax = plt.subplots(figsize=figsize)
        lengths.groupby("group").mean().transpose().plot(ax=ax)
        plt.xlabel("n_propagation")
        plt.ylabel(f"Mean delta_X length")
        plt.legend(
            bbox_to_anchor=(1.05, 0), loc="lower left", borderaxespad=0, fontsize=13
        )
        plt.subplots_adjust(right=0.5, bottom=0.15)

        return fig

    def calculate_p_mass(self, smooth=0.8, n_grid=40, n_neighbors=200, n_jobs=-1):
        self.calculate_grid_arrows(
            smooth=0.8, steps=(n_grid, n_grid), n_neighbors=n_neighbors, n_jobs=-1
        )

    def suggest_mass_thresholds(self, n_suggestion=12, s=1, n_col=4):
        min_ = self.total_p_mass.min()
        max_ = self.total_p_mass.max()
        suggestions = np.linspace(min_, max_ / 2, n_suggestion)

        n_rows = math.ceil(n_suggestion / n_col)

        fig, ax = plt.subplots(n_rows, n_col, figsize=[5 * n_col, 5 * n_rows])
        if n_rows == 1:
            ax = ax.reshape(1, -1)

        row = 0
        col = 0
        for i in range(n_suggestion):
            ax_ = ax[row, col]

            col += 1
            if col == n_col:
                col = 0
                row += 1

            idx = self.total_p_mass > suggestions[i]

            # ax_.scatter(gridpoints_coordinates[mass_filter, 0], gridpoints_coordinates[mass_filter, 1], s=0)
            ax_.scatter(self.embedding[:, 0], self.embedding[:, 1], c="lightgray", s=s)
            ax_.scatter(self.flow_grid[idx, 0], self.flow_grid[idx, 1], c="black", s=s)
            ax_.set_title(f"min_mass: {suggestions[i]: .2g}")
            ax_.axis("off")

    def calculate_mass_filter(self, min_mass=0.01, plot=False):
        self.min_mass = min_mass
        self.mass_filter = self.total_p_mass < min_mass

        if plot:
            fig, ax = plt.subplots(figsize=[5, 5])

            # ax_.scatter(gridpoints_coordinates[mass_filter, 0], gridpoints_coordinates[mass_filter, 1], s=0)
            ax.scatter(self.embedding[:, 0], self.embedding[:, 1], c="lightgray", s=10)
            ax.scatter(
                self.flow_grid[~self.mass_filter, 0],
                self.flow_grid[~self.mass_filter, 1],
                c="black",
                s=0.5,
            )
            ax.set_title("Grid points selected")
            ax.axis("off")

    ## Get randomized GRN coef to do randomized perturbation simulation
    def calculate_randomized_coef_table(self, random_seed=123):
        "Calculate randomized GRN coef table."

        if hasattr(self, "coef_matrix_per_cluster"):
            coef_matrix_per_cluster_randomized = {}
            for key, val in self.coef_matrix_per_cluster.items():
                coef_matrix_per_cluster_randomized[
                    key
                ] = _shuffle_celloracle_GRN_coef_table(
                    coef_dataframe=val, random_seed=random_seed
                )
            self.coef_matrix_per_cluster_randomized = coef_matrix_per_cluster_randomized

        if hasattr(self, "coef_matrix"):
            self.coef_matrix_randomized = _shuffle_celloracle_GRN_coef_table(
                coef_dataframe=self.coef_matrix, random_seed=random_seed
            )

        if (hasattr(self, "coef_matrix_per_cluster") == False) and (
            hasattr(self, "coef_matrix") == False
        ):
            print(
                "GRN calculation for simulation is not finished. Run fit_GRN_for_simulation() first."
            )

    ########################################
    ### 4. Methods for Markov simulation ###
    ########################################
    def prepare_markov_simulation(self, verbose=False):
        """
        Pick up cells for Markov simulation.

        Args:
            verbose (bool): If True, it plots selected cells.

        """
        # Sample uniformly the points to avoid density driven effects - Should reimplement as a method
        steps = 100, 100
        grs = []
        for dim_i in range(self.embedding.shape[1]):
            m, M = np.min(self.embedding[:, dim_i]), np.max(self.embedding[:, dim_i])
            m = m - 0.025 * np.abs(M - m)
            M = M + 0.025 * np.abs(M - m)
            gr = np.linspace(m, M, steps[dim_i])
            grs.append(gr)

        meshes_tuple = np.meshgrid(*grs)
        gridpoints_coordinates = np.vstack([i.flat for i in meshes_tuple]).T

        nn = NearestNeighbors()
        nn.fit(self.embedding)
        dist, ixs = nn.kneighbors(gridpoints_coordinates, 1)

        diag_step_dist = np.sqrt(
            (meshes_tuple[0][0, 0] - meshes_tuple[0][0, 1]) ** 2
            + (meshes_tuple[1][0, 0] - meshes_tuple[1][1, 0]) ** 2
        )
        min_dist = diag_step_dist / 2
        ixs = ixs[dist < min_dist]
        gridpoints_coordinates = gridpoints_coordinates[dist.flat[:] < min_dist, :]
        dist = dist[dist < min_dist]

        ixs = np.unique(ixs)
        self.ixs_mcmc = ixs

        if verbose:
            plt.scatter(
                self.embedding[ixs, 0],
                self.embedding[ixs, 1],
                c=self.colorandum[ixs],
                alpha=1,
                s=30,
                lw=0.4,
                edgecolor="0.4",
            )

        self.prepare_markov(
            sigma_D=diag_step_dist,
            sigma_W=diag_step_dist / 2.0,
            direction="forward",
            cells_ixs=ixs,
        )

    def run_markov_chain_simulation(
        self, n_steps=500, n_duplication=5, seed=123, calculate_randomized=True
    ):
        """
        Do Markov simlations to predict cell transition after perturbation.
        The transition probability between cells has been calculated
        based on simulated gene expression values in the signal propagation process.
        The cell state transition will be simulated based on the probability.
        You can simulate the process multiple times to get a robust outcome.

        Args:
            n_steps (int): steps for Markov simulation. This value is equivalent to the amount of time after perturbation.

            n_duplication (int): the number for multiple calculations.

        """
        np.random.seed(seed)
        # _numba_random_seed(seed)

        self.prepare_markov_simulation()

        transition_prob = self.tr.toarray()

        #
        transition_prob = _deal_with_na(transition_prob)  # added 20200607

        n_cells = transition_prob.shape[0]

        start_cell_id_array = np.repeat(np.arange(n_cells), n_duplication)

        transition = _walk(start_cell_id_array, transition_prob, n_steps)
        transition = self.ixs_mcmc[transition]

        li = None

        ind = np.repeat(self.ixs_mcmc, n_duplication)
        self.mcmc_transition_id = pd.DataFrame(transition, ind)

        if calculate_randomized:
            transition_prob_random = self.tr_random.toarray()
            #
            transition_prob_random = _deal_with_na(
                transition_prob_random
            )  # added 20200607

            n_cells = transition_prob_random.shape[0]

            start_cell_id_array = np.repeat(np.arange(n_cells), n_duplication)

            transition_random = _walk(
                start_cell_id_array, transition_prob_random, n_steps
            )
            transition_random = self.ixs_mcmc[transition_random]

            li = None

            ind = np.repeat(self.ixs_mcmc, n_duplication)
            self.mcmc_transition_random_id = pd.DataFrame(transition_random, ind)

    def summarize_mc_results_by_cluster(self, cluster_use, random=False):
        """
        This function summarizes the simulated cell state-transition by groping the results into each cluster.
        It returns sumarized results as a pandas.DataFrame.

        Args:
            cluster_use (str): cluster information name in anndata.obs.
               You can use any arbitrary cluster information in anndata.obs.
        """
        if random:
            transition = self.mcmc_transition_random_id.values
        else:
            transition = self.mcmc_transition_id.values

        mcmc_transition_cluster = np.array(self.adata.obs[cluster_use])[transition]
        mcmc_transition_cluster = pd.DataFrame(
            mcmc_transition_cluster, index=self.mcmc_transition_id.index
        )
        return mcmc_transition_cluster

    def plot_mc_results_as_sankey(
        self, cluster_use, start=0, end=-1, order=None, font_size=10
    ):
        """
        Plot the simulated cell state-transition as a Sankey-diagram after groping by the cluster.

        Args:
            cluster_use (str): cluster information name in anndata.obs.
               You can use any cluster information in anndata.obs.

            start (int): The starting point of Sankey-diagram. Please select a  step in the Markov simulation.

            end (int): The end point of Sankey-diagram. Please select a  step in the Markov simulation.
                if you set [end=-1], the final step of Markov simulation will be used.

            order (list of str): The order of cluster name in the Sankey-diagram.

            font_size (int): Font size for cluster name label in the Sankey diagram.

        """
        mcmc_transition_cluster = self.summarize_mc_results_by_cluster(cluster_use)
        mcmc_color_dict = _adata_to_color_dict(self.adata, cluster_use)

        df = mcmc_transition_cluster.iloc[:, [start, end]]
        df.columns = ["start", "end"]

        if not order is None:
            order_ = order.copy()
            order_.reverse()
            order_left = [i for i in order_ if i in df.start.unique()]
            order_right = [i for i in order_ if i in df.end.unique()]
        else:
            order_left = list(df.start.unique())
            order_right = list(df.end.unique())

        sankey(
            left=df["start"],
            right=df["end"],
            aspect=2,
            fontsize=font_size,
            colorDict=mcmc_color_dict,
            leftLabels=order_left,
            rightLabels=order_right,
        )

    def plot_mc_results_as_kde(self, n_time, args={}):
        """
        Pick up one timepoint in the cell state-transition simulation and plot as a kde plot.

        Args:
            n_time (int): the number in Markov simulation

            args (dictionary): An argument for seaborn.kdeplot.
                See seaborn documentation for details (https://seaborn.pydata.org/generated/seaborn.kdeplot.html#seaborn.kdeplot).

        """
        cell_ix = self.mcmc_transition_id.iloc[:, n_time].values

        x = self.embedding[cell_ix, 0]
        y = self.embedding[cell_ix, 1]

        sns.kdeplot(x, y, **args)

    def plot_mc_results_as_trajectory(self, cell_name, time_range, args={}):
        """
        Pick up several timepoints in the cell state-transition simulation and plot as a line plot.
        This function can be used to visualize how cell-state changes after perturbation focusing on a specific cell.

        Args:
            cell_name (str): cell name. chose from adata.obs.index

            time_range (list of int): the list of index in Markov simulation

            args (dictionary): dictionary for the arguments for matplotlib.pyplit.plot.
                See matplotlib documentation for details (https://matplotlib.org/api/_as_gen/matplotlib.pyplot.plot.html#matplotlib.pyplot.plot).

        """
        cell_ix = np.where(self.adata.obs.index == cell_name)[0][0]
        cell_ix_in_mcmctid = np.where(self.mcmc_transition_id.index == cell_ix)[0]

        # plot all cells in gray color
        plt.scatter(self.embedding[:, 0], self.embedding[:, 1], s=1, c="lightgray")

        for i in cell_ix_in_mcmctid:
            self._plot_one_trajectory(i, time_range, args)

        # plot cell of interest (initiation point of simulation) in red color
        plt.scatter(
            self.embedding[cell_ix, 0], self.embedding[cell_ix, 1], s=50, c="red"
        )

    def _plot_one_trajectory(self, cell_ix_in_mcmctid, time_range, args={}):
        tt = self.mcmc_transition_id.iloc[cell_ix_in_mcmctid, :].values[time_range]
        plt.plot(self.embedding[:, 0][tt], self.embedding[:, 1][tt], **args)

    def count_cells_in_mc_resutls(self, cluster_use, end=-1, order=None):
        """
        Count the simulated cell by the cluster.

        Args:
            cluster_use (str): cluster information name in anndata.obs.
               You can use any cluster information in anndata.obs.

            end (int): The end point of Sankey-diagram. Please select a  step in the Markov simulation.
                if you set [end=-1], the final step of Markov simulation will be used.
        Returns:
            pandas.DataFrame : Number of cells before / after simulation

        """
        mcmc_transition_cluster = self.summarize_mc_results_by_cluster(
            cluster_use, random=False
        )

        if hasattr(self, "mcmc_transition_random_id"):
            mcmc_transition_cluster_random = self.summarize_mc_results_by_cluster(
                cluster_use, random=True
            )

            df = pd.DataFrame(
                {
                    "original": mcmc_transition_cluster.iloc[:, 0],
                    "simulated": mcmc_transition_cluster.iloc[:, end],
                    "randomized": mcmc_transition_cluster_random.iloc[:, end],
                }
            )
        else:
            df = pd.DataFrame(
                {
                    "original": mcmc_transition_cluster.iloc[:, 0],
                    "simulated": mcmc_transition_cluster.iloc[:, end],
                }
            )

        # Post processing
        n_duplicated = df.index.value_counts().values[0]
        df["simulation_batch"] = [i % n_duplicated for i in np.arange(len(df))]

        df = df.melt(id_vars="simulation_batch")
        df["count"] = 1
        df = df.groupby(["value", "variable", "simulation_batch"]).count()
        df = df.reset_index(drop=False)

        df = df.rename(columns={"value": "cluster", "variable": "data"})
        df["simulation_batch"] = df["simulation_batch"].astype(np.object)

        return df

    def get_mcmc_cell_transition_table(self, cluster_column_name=None, end=-1):
        """
        Return cell count in the initial state and final state after mcmc.
        Cell counts are grouped by the cluster of interest.
        Result will be returned as 2D matrix.
        """

        if cluster_column_name is None:
            cluster_column_name = self.cluster_column_name

        start = 0

        mcmc_transition = self.summarize_mc_results_by_cluster(
            cluster_column_name, random=False
        )
        mcmc_transition = mcmc_transition.iloc[:, [start, end]]
        mcmc_transition.columns = ["start", "end"]
        mcmc_transition["count"] = 1
        mcmc_transition = pd.pivot_table(
            mcmc_transition,
            values="count",
            index=["start"],
            columns=["end"],
            aggfunc=np.sum,
            fill_value=0,
        )

        mcmc_transition_random = self.summarize_mc_results_by_cluster(
            cluster_column_name, random=True
        )
        mcmc_transition_random = mcmc_transition_random.iloc[:, [start, end]]
        mcmc_transition_random.columns = ["start", "end"]
        mcmc_transition_random["count"] = 1
        mcmc_transition_random = pd.pivot_table(
            mcmc_transition_random,
            values="count",
            index=["start"],
            columns=["end"],
            aggfunc=np.sum,
            fill_value=0,
        )

        # store data
        mcmc_transition_random.index.name = None
        mcmc_transition_random.columns.name = None
        mcmc_transition.index.name = None
        mcmc_transition.columns.name = None

        self.mcmc_transition = mcmc_transition
        self.mcmc_transition_random = mcmc_transition_random

    ###################################################
    ### 5. GRN inference for Network score analysis ###
    ###################################################
    def get_links(
        self,
        cluster_name_for_GRN_unit=None,
        alpha=10,
        bagging_number=20,
        verbose_level=1,
        test_mode=False,
        model_method="bagging_ridge",
        ignore_warning=False,
        n_jobs=-1,
    ):
        """
        Makes GRN for each cluster and returns results as a Links object.
        Several preprocessing should be done before using this function.

        Args:
            cluster_name_for_GRN_unit (str): Cluster name for GRN calculation. The cluster information should be stored in Oracle.adata.obs.

            alpha (float or int): The strength of regularization.
                If you set a lower value, the sensitivity increases, and you can detect weaker network connections. However, there may be more noise.
                If you select a higher value, it will reduce the chance of overfitting.

            bagging_number (int): The number used in bagging calculation.


            verbose_level (int): if [verbose_level>1], most detailed progress information will be shown.
                if [verbose_level > 0], one progress bar will be shown.
                if [verbose_level == 0], no progress bar will be shown.

            test_mode (bool): If test_mode is True, GRN calculation will be done for only one cluster rather than all clusters.

            model_method (str): Chose modeling algorithm. "bagging_ridge" or "bayesian_ridge"

            n_jobs (int): Number of cpu cores for parallel calculation.  -1 means using all available cores. Default is -1.


        """

        ## Check data
        info = self._generate_meta_data()

        if ignore_warning:
            pass
        else:
            if info["status - Gene expression matrix"] != "Ready":
                raise ValueError("scRNA-seq data is not imported.")

            if info["status - PCA calculation"] != "Done":
                raise ValueError(
                    "Preprocessing is not done. Do PCA and Knn imputation."
                )

            if info["status - Knn imputation"] != "Done":
                raise ValueError("Preprocessing is not done. Do Knn imputation.")

            if info["status - BaseGRN"] != "Ready":
                raise ValueError(
                    "Found No TF information. Please import TF data (base-GRN) first."
                )

            if info["n_regulatory_in_both_TFdict_and_scRNA-seq"] == "0 genes":
                raise ValueError(
                    "Found No overlap between TF info (base GRN) and your scRNA-seq data. Please check your data format and species."
                )

            if info["n_target_genes_both_TFdict_and_scRNA-seq"] == "0 genes":
                raise ValueError(
                    "Found No overlap between TF info (base GRN) and your scRNA-seq data. Please check your data format and species."
                )

        links = get_links(
            oracle_object=self,
            cluster_name_for_GRN_unit=cluster_name_for_GRN_unit,
            alpha=alpha,
            bagging_number=bagging_number,
            verbose_level=verbose_level,
            test_mode=test_mode,
            model_method=model_method,
            n_jobs=n_jobs,
        )
        return links


def _deal_with_na(transition_prob):
    tr = transition_prob.copy()

    # remove nan
    tr = np.nan_to_num(tr, copy=True, nan=0)

    # if transition prob is 0 in all row, assign transitionprob = 1 to self row.
    no_transition_ids = tr.sum(axis=1) == 0
    tr[no_transition_ids, no_transition_ids] = 1

    return tr
