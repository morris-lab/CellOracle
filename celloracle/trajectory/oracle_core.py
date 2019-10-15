# -*- coding: utf-8 -*-

import logging
import warnings
from copy import deepcopy

import matplotlib.pyplot as plt
import numpy as np
from typing import Dict, Any, List, Union, Tuple
import pandas as pd
import scanpy as sc
import seaborn as sns
from ..utility.hdf5_processing import dump_hdf5, load_hdf5

from sklearn.neighbors import NearestNeighbors

from .sankey import sankey
from .markov_simulation import _walk
from .oracle_utility import (_adata_to_matrix, _adata_to_df,
                             _adata_to_color_dict, _get_clustercolor_from_anndata,
                             _numba_random_seed, _linklist2dict)
from .oracle_GRN import _do_simulation, _getCoefMatrix
from .modified_VelocytoLoom_class import modified_VelocytoLoom
from ..network_analysis.network_construction import get_links



def load_oracle(file_path):
    """
    Load oracle object saved as hdf5 file.

    Args:
        file_path (str): File path to the hdf5 file.

    """

    return load_hdf5(filename=file_path, obj_class=Oracle)


class Oracle(modified_VelocytoLoom):
    """
    Oracle is the main class in CellOracle.
    Oracle object imports scRNA-seq data (anndata) and TF information to infer cluster-specific GRNs.
    It can predict future gene expression patterns and cell state transition after perturbations of TFs.
    Please see the paper of CellOracle for details.

    The code of the Oracle class was made of three components below.

    (1) Anndata:  Gene expression matrix and metadata in single-cell RNA-seq are stored in anndata object.
    Processed values, such as normalized counts and simulated values, are stored as layers of anndata.
    Metadata (i.e., Cluster info) are saved in anndata.obs. Refer to scanpy/anndata documentation for detail.


    (2) Net: Net is a custom class in celloracle.
    Net object process several data to infer GRN.
    See the documentation of Net class for detail.

    (3) VelycytoLoom: Calculation of transition probability and visualization of directed trajectory graph will be performed in the same way as velocytoloom.
    VelocytoLoom is class for the Velocyto, which is a python library for RNA-velocity analysis.
    In the celloracle, we use almost the same functions as velocytoloom for the visualization, but celloracle use simulated gene expression values instead of RNA-velocity data.

    Some CellOracle's methods were inspired by velocyto analysis and some codes were made by modifying VelocytoLoom class.

    Attributes:
        adata (anndata): Imported anndata object
        cluster_column_name (str): The column name in adata.obs about cluster info
        embedding_name (str): The key name in adata.obsm about dimensional reduction cordinates

    """

    def __init__(self):

        self.adata = None

        self.cluster_column_name = None
        self.embedding_name = None
        self.ixs_mcmc = None
        self.cluster_specific_TFdict = None
        self.cv_mean_selected_genes = None

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
        dump_hdf5(obj=self, filename=file_path,
                  data_compression=compression_opts,  chunks=(2048, 2048),
                  noarray_compression=compression_opts, pickle_protocol=2)


    ###################################
    ### 1. Methods for loading data ###
    ###################################

    def import_TF_data(self, TF_info_matrix=None, TF_info_matrix_path=None, TFdict=None):
        """
        Load data about potential-regulatory TFs.
        You can import either TF_info_matrix or TFdict.
        See the tutorial of celloracle or motif_analysis module for an example to make such files.

        Args:
            TF_info_matrix (pandas.DataFrame): TF_info_matrix.

            TF_info_matrix_path (str): File path for TF_info_matrix (pandas.DataFrame).

            TFdict (dictionary): Python dictionary of TF info.
        """
        if not TF_info_matrix is None:
            tmp = TF_info_matrix.copy()
            tmp = tmp.drop(["peak_id"], axis=1)
            tmp = tmp.groupby(by="gene_short_name").sum()
            self.TFdict = dict(tmp.apply(lambda x: x[x>0].index.values, axis=1))

        if not TF_info_matrix_path is None:
            tmp = pd.read_parquet(TF_info_matrix_path)
            tmp = tmp.drop(["peak_id"], axis=1)
            tmp = tmp.groupby(by="gene_short_name").sum()
            self.TFdict = dict(tmp.apply(lambda x: x[x>0].index.values, axis=1))

        if not TFdict is None:
            self.TFdict=TFdict.copy()

    def updateTFinfo_dictionary(self, TFdict):
        """
        Update a TF dictionary.
        If a key in the new TF dictionary already existed in the old TF dictionary, old values will be replaced with a new one.

        Args:
            TFdict (dictionary): Python dictionary of TF info.
        """

        self.TFdict.update(TFdict)

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

    def get_cluster_specific_TFdict_from_Links(self, links_object):

        """
        Extract TF and its target gene information from Links object.
        This function can be used to reconstruct GRNs based on pre-existing GRNs saved in Links object.

        Args:
            links_object (Links): Please see the explanation of Links class.

        """
        self.cluster_specific_TFdict = {}

        for i in links_object.filtered_links:
            self.cluster_specific_TFdict[i] = _linklist2dict(links_object.filtered_links[i])

    def import_anndata_as_raw_count(self, adata, cluster_column_name=None, embedding_name=None,
                                    transform="natural_log"):
        """
        Load scRNA-seq data. scRNA-seq data should be prepared as an anndata.
        Preprocessing (cell and gene filtering, calculate DR and cluster, etc.) should be done before loading data.
        The method imports RAW GENE COUNTS because unscaled and uncentered gene expression data are required for the GRN inference and simulation.
        See tutorial notebook for the details about how to process scRNA-seq data.

        Args:
            adata (anndata): anndata object that stores scRNA-seq data.

            cluster_column_name (str): the column name about cluster info in anndata.obs.
                Clustering data suppose to be in anndata.obs.

            embedding_name (str): the key name about a dimensional reduction in anndata.obsm.
                Dimensional reduction (or 2D trajectory graph) should be in anndata.obsm.

            transform (str): The method for log-transformation. Chose one from "natural_log" or "log2".

        """
        if adata.X.min() < 0:
            raise ValueError("gene expression matrix (adata.X) does not seems to be raw_count because it contains negavive values.")

        # store data
        self.adata = adata.copy()
        self.cluster_column_name = cluster_column_name
        self.embedding_name = embedding_name

        #if hasattr(self.adata, "raw"):
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
        col_dict = _get_clustercolor_from_anndata(adata=self.adata,
                                                  cluster_name=self.cluster_column_name,
                                                  return_as="dict")
        self.colorandum = np.array([col_dict[i] for i in self.adata.obs[self.cluster_column_name]])

        # variable gene detection for the QC of simulation
        N = adata.shape[1]
        if N >= 3000:
            N = 3000
        self.score_cv_vs_mean(int(N/3)-1, plot=False, max_expr_avg=35)
        self.high_var_genes = self.cv_mean_selected_genes.copy()
        self.cv_mean_selected_genes = None

    def import_anndata_as_normalized_count(self, adata, cluster_column_name=None, embedding_name=None):
        """
        Load scRNA-seq data. scRNA-seq data should be prepared as an anndata.
        Preprocessing (cell and gene filtering, calculate DR and cluster, etc.) should be done before loading data.
        The method will import NORMALIZED and LOG TRANSFORMED but NOT SCALED and NOT CENTERED DATA.
        See the tutorial for the details for an example of how to process scRNA-seq data.

        Args:
            adata (anndata): anndata object that store scRNA-seq data.

            cluster_column_name (str): the column name about cluster info in anndata.obs.
                Clustering data suppose to be in anndata.obs.

            embedding_name (str): the key name about a dimensional reduction in anndata.obsm.
                Dimensional reduction (or 2D trajectory graph) should be in anndata.obsm.

            transform (str): The method for log-transformation. Chose one from "natural_log" or "log2".
        """
        if adata.X.min() < 0:
            raise ValueError("gene expression matrix (adata.X) contains negavive values. Please use UNSCALED and UNCENTERED data.")

        # store data
        self.adata = adata.copy()
        self.cluster_column_name = cluster_column_name
        self.embedding_name = embedding_name

        # store raw count data
        #self.adata.layers["raw_count"] = adata.X.copy()

        # normalization and log transformation
        self.adata.layers["normalized_count"] = self.adata.X.copy()

        # update color information
        col_dict = _get_clustercolor_from_anndata(adata=self.adata,
                                                  cluster_name=self.cluster_column_name,
                                                  return_as="dict")
        self.colorandum = np.array([col_dict[i] for i in self.adata.obs[self.cluster_column_name]])

        # variable gene detection for the QC of simulation
        N = adata.shape[1]
        if N >= 3000:
            N = 3000
        self.score_cv_vs_mean(int(N/3)-1, plot=False, max_expr_avg=35)
        self.high_var_genes = self.cv_mean_selected_genes.copy()
        self.cv_mean_selected_genes = None



    ####################################
    ### 2. Methods for GRN inference ###
    ####################################
    def fit_GRN_for_simulation(self, GRN_unit="cluster", alpha=1, use_cluster_specific_TFdict=False):
        """
        Do GRN inference.
        Please see the paper of CellOracle for details.

        GRN can be constructed at an arbitrary cell group.
        If you want to infer cluster-specific GRN, please set [GRN_unit="cluster"].
        GRN will be inferred for each cluster. You can select Cluster information when you import data (not when you run this method.).

        If you set [GRN_unit="whole"], GRN will be made using all cells.

        Args:
            GRN_unit (str): select "cluster" or "whole"

            alpha (float or int): the strength of regularization.
                If you set a lower value, the sensitivity increase, and you can detect a weak network connection, but it might get more noize.
                With a higher value of alpha may reduce the chance of overfitting.
        """
        # prepare data for GRN calculation
        gem_imputed = _adata_to_df(self.adata, "imputed_count")

        self.adata.layers["simulation_input"] = self.adata.layers["imputed_count"].copy()
        self.alpha_for_trajectory_GRN = alpha

        if use_cluster_specific_TFdict & (self.cluster_specific_TFdict is not None):
            self.coef_matrix_per_cluster = {}
            cluster_info = self.adata.obs[self.cluster_column_name]

            print(f"calculating GRN using cluster specicif TF dict...")
            for cluster in np.unique(cluster_info):
                print(f"calculating GRN in {cluster}")
                cells_in_the_cluster_bool = (cluster_info == cluster)
                gem_ = gem_imputed[cells_in_the_cluster_bool]
                self.coef_matrix_per_cluster[cluster] = _getCoefMatrix(gem=gem_,
                                                                       TFdict=self.cluster_specific_TFdict[cluster],
                                                                       alpha=alpha)


        else:
            if GRN_unit == "whole":
                self.coef_matrix = _getCoefMatrix(gem=gem_imputed, TFdict=self.TFdict, alpha=alpha)
            if GRN_unit == "cluster":
                self.coef_matrix_per_cluster = {}
                cluster_info = self.adata.obs[self.cluster_column_name]
                for cluster in np.unique(cluster_info):
                    print(f"calculating GRN in {cluster}")
                    cells_in_the_cluster_bool = (cluster_info == cluster)
                    gem_ = gem_imputed[cells_in_the_cluster_bool]
                    self.coef_matrix_per_cluster[cluster] = _getCoefMatrix(gem=gem_,
                                                                           TFdict=self.TFdict,
                                                                           alpha=alpha)



    #######################################################
    ### 3. Methods for simulation of signal propagation ###
    #######################################################

    def simulate_shift(self, perturb_condition=None, GRN_unit="whole",
                       n_propagation=3):
        """
        Simulate signal propagation with GRNs. Please see the paper of CellOracle for details.
        This function simulates a gene expression pattern in the near future.
        Simulated values will be stored in anndata.layers: ["simulated_count"]


        Three data below are used for the simulation.
        (1) GRN inference results (coef_matrix).
        (2) perturb_condition: You can set arbitrary perturbation condition.
        (3) gene expression matrix: simulation starts from imputed gene expression data.

        Args:
            perturb_condition (dictionary): condition for perturbation.
               if you want to simulate knockout for GeneX, please set [perturb_condition={"GeneX": 0.0}]
               Although you can set any non-negative values for the gene condition, avoid setting biologically unfeasible values for the perturb condition.
               It is strongly recommended to check actual gene expression values in your data before selecting perturb condition.

            n_propagation (int): Calculation will be performed iteratively to simulate signal propagation in GRN.
                you can set the number of steps for this calculation.
                With a higher number, the results may recapitulate signal propagation for many genes.
                However, a higher number of propagation may cause more error/noise.
        """

        # 0. Reset previous simulation results if it exist
        self.ixs_mcmc = None
        self.mcmc_transition_id = None
        self.corrcoef = None
        self.transition_prob = None
        self.tr = None


        # 1. prepare perturb information
        if not perturb_condition is None:

            for i in perturb_condition.keys():
                if i not in self.high_var_genes:
                    print(f"Variability score of Gene {i} is low. Simulation accuracy may be poor with this gene.")

            # reset simulation initiation point
            self.adata.layers["simulation_input"] = self.adata.layers["imputed_count"].copy()
            simulation_input = _adata_to_df(self.adata, "simulation_input")

            for i in perturb_condition.keys():
                if not i in simulation_input.columns:
                    raise ValueError(f"{i} is not in the data")
                else:
                    simulation_input[i] = perturb_condition[i]
        else:
            simulation_input = _adata_to_df(self.adata, "simulation_input")

        # 2. load gene expression matrix (initiation information for the simulation)
        gem_imputed = _adata_to_df(self.adata, "imputed_count")

        # 3. do simulation for signal propagation within GRNs
        if GRN_unit == "whole":
            gem_simulated = _do_simulation(self.coef_matrix,
                                           simulation_input,
                                           gem_imputed,
                                           n_propagation)

        elif GRN_unit == "cluster":
            simulated = []
            cluster_info = self.adata.obs[self.cluster_column_name]
            for cluster in np.unique(cluster_info):
                cells_in_the_cluster_bool = (cluster_info == cluster)

                simulation_input_ = simulation_input[cells_in_the_cluster_bool]
                gem_ = gem_imputed[cells_in_the_cluster_bool]

                simulated_in_the_cluster = _do_simulation(
                                             self.coef_matrix_per_cluster[cluster],
                                             simulation_input_,
                                             gem_,
                                             n_propagation)
                simulated.append(simulated_in_the_cluster)
            gem_simulated = pd.concat(simulated, axis=0)
            gem_simulated = gem_simulated.reindex(gem_imputed.index)

        else:
            raise ValueError("GRN_unit shold be either of 'whole' or 'cluster'")

        # 4. store simulation results
        #  simulated future gene expression matrix
        self.adata.layers["simulated_count"] = gem_simulated.values

        #  difference between simulated values and original values
        self.adata.layers["delta_X"] = self.adata.layers["simulated_count"] - self.adata.layers["imputed_count"]


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

        diag_step_dist = np.sqrt((meshes_tuple[0][0,0] - meshes_tuple[0][0,1])**2 + (meshes_tuple[1][0,0] - meshes_tuple[1][1,0])**2)
        min_dist = diag_step_dist / 2
        ixs = ixs[dist < min_dist]
        gridpoints_coordinates = gridpoints_coordinates[dist.flat[:]<min_dist,:]
        dist = dist[dist < min_dist]

        ixs = np.unique(ixs)
        self.ixs_mcmc = ixs

        if verbose:
            plt.scatter(self.embedding[ixs, 0], self.embedding[ixs, 1],
                        c=self.colorandum[ixs], alpha=1, s=30, lw=0.4,
                        edgecolor="0.4")

        self.prepare_markov(sigma_D=diag_step_dist, sigma_W=diag_step_dist/2.,
                       direction='forward', cells_ixs=ixs)


    def run_markov_chain_simulation(self, n_steps=500, n_duplication=5, seed=123):
        """
        Do Markov simlation to predict cell transition after perturbation.
        The transition probability between cells has been calculated
        based on simulated gene expression values in the signal propagation process.
        The cell state transition will be simulated based on the probability.
        You can simulate the process for multiple times to get a robust outcome.

        Args:
            n_steps (int): steps for Markov simulation. This value is equivalent to the time after perturbation.

            n_duplication (int): the number for multiple calculations.

        """
        np.random.seed(seed)
        _numba_random_seed(seed)

        self.prepare_markov_simulation()

        transition_prob = self.tr.toarray()
        n_cells = transition_prob.shape[0]

        start_cell_id_array = np.repeat(np.arange(n_cells), n_duplication)

        transition = _walk(start_cell_id_array, transition_prob, n_steps)
        transition = self.ixs_mcmc[transition]

        li = None

        ind = np.repeat(self.ixs_mcmc, n_duplication)
        self.mcmc_transition_id = pd.DataFrame(transition, ind)

    def summarize_mc_results_by_cluster(self, cluster_use):
        """
        This function summarizes the simulated cell state-transition by groping the results into each cluster.
        It returns sumarized results as a pandas.DataFrame.

        Args:
            cluster_use (str): cluster information name in anndata.obs.
               You can use any arbitrary cluster information in anndata.obs.
        """
        transition = self.mcmc_transition_id.values
        mcmc_transition_cluster = np.array(self.adata.obs[cluster_use])[transition]
        mcmc_transition_cluster = pd.DataFrame(mcmc_transition_cluster,
                                               index=self.mcmc_transition_id.index)
        return mcmc_transition_cluster


    def plot_mc_resutls_as_sankey(self, cluster_use, start=0, end=-1, order=None, font_size=10):
        """
        Plot the simulated cell state-transition as a Sankey-diagram after groping by the cluster.

        Args:
            cluster_use (str): cluster information name in anndata.obs.
               You can use any arbitrary cluster information in anndata.obs.

            start (int): The starting point of Sankey-diagram. Please select a  step in the Markov simulation.

            end (int): The end point of Sankey-diagram. Please select a  step in the Markov simulation.
                if you set [end=-1], the final step of Markov simulation will be used.

            order (list of str): The order of cluster name in sankey-diagram.

            font_size (int): Font size for cluster name label in Sankey diagram.

        """
        mcmc_transition_cluster = self.summarize_mc_results_by_cluster(cluster_use)
        mcmc_color_dict =  _adata_to_color_dict(self.adata, cluster_use)

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

        sankey(left=df['start'], right=df['end'],
               aspect=2, fontsize=font_size,
               colorDict=mcmc_color_dict,
               leftLabels=order_left, rightLabels=order_right)


    def plot_mc_result_as_kde(self, n_time, args={}):
        """
        Pick up one timepoint in the cell state-transition simulation and plot as a kde plot.

        Args:
            n_time (int): the number in Markov simulation

            args (dictionary): An argument for seaborn.kdeplot.
                See seaborn documentation for detail (https://seaborn.pydata.org/generated/seaborn.kdeplot.html#seaborn.kdeplot).

        """
        cell_ix = self.mcmc_transition_id.iloc[:, n_time].values

        x = self.embedding[cell_ix, 0]
        y = self.embedding[cell_ix, 1]

        sns.kdeplot(x, y, **args)

    def plot_mc_result_as_trajectory(self, cell_name, time_range, args={}):
        """
        Pick up several timepoints in the cell state-transition simulation and plot as a line plot.
        This function can be used to visualize how cell-state changes after perturbation focusing on a specific cell.

        Args:
            cell_name (str): cell name. chose from adata.obs.index

            time_range (list of int): the number in markov simulation

            args (dictionary): dictionary for the arguments for matplotlib.pyplit.plot.
                See matplotlib documentation for detail (https://matplotlib.org/api/_as_gen/matplotlib.pyplot.plot.html#matplotlib.pyplot.plot).

        """
        cell_ix = np.where(self.adata.obs.index == cell_name)[0][0]
        cell_ix_in_mcmctid = np.where(self.mcmc_transition_id.index == cell_ix)[0]

        # plot all cells in gray color
        plt.scatter(self.embedding[:,0], self.embedding[:,1], s=1, c="lightgray")


        for i in cell_ix_in_mcmctid:
            self._plot_one_trajectory(i, time_range, args)

        # plot cell of interest (initiation point of simulation) in red color
        plt.scatter(self.embedding[cell_ix,0], self.embedding[cell_ix,1], s=50, c="red")

    def _plot_one_trajectory(self, cell_ix_in_mcmctid, time_range, args={}):
        tt = self.mcmc_transition_id.iloc[cell_ix_in_mcmctid,:].values[time_range]
        plt.plot(self.embedding[:,0][tt], self.embedding[:,1][tt], **args)


    ###################################################
    ### 5. GRN inference for Network score analysis ###
    ###################################################
    def get_links(self, cluster_name_for_GRN_unit=None, alpha=10, bagging_number=20, verbose_level=1, test_mode=False):
        """
        Make GRN for each cluster and returns results as a Links object.
        Several preprocessing should be done before using this function.

        Args:
            cluster_name_for_GRN_unit (str): Cluster name for GRN calculation. The cluster information should be stored in Oracle.adata.obs.

            alpha (float or int): the strength of regularization.
                If you set a lower value, the sensitivity increase, and you can detect a weak network connection, but it might get more noize.
                With a higher value of alpha may reduce the chance of overfitting.

            bagging_number (int): The number for bagging calculation.


            verbose_level (int): if [verbose_level>1], most detailed progress information will be shown.
                if [verbose_level > 0], one progress bar will be shown.
                if [verbose_level == 0], no progress bar will be shown.

            test_mode (bool): If test_mode is True, GRN calculation will be done for only one cluster rather than all clusters.

        """
        links = get_links(oracle_object=self,
                          cluster_name_for_GRN_unit=cluster_name_for_GRN_unit,
                          alpha=alpha, bagging_number=bagging_number,
                          verbose_level=verbose_level, test_mode=test_mode)
        return links
