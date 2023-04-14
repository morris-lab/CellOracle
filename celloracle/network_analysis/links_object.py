# -*- coding: utf-8 -*-
"""

"""

###########################
### 0. Import libraries ###
###########################


# 0.1. libraries for fundamental data science and data processing

import warnings
import os
import pandas as pd
import numpy as np
from scipy import stats
from copy import deepcopy

from igraph import Graph

from ..utility.hdf5_processing import dump_hdf5, load_hdf5

from .use_r_scripts import (
    _get_network_score_by_Rscripts_inparallel,
    _check_R_libraries_installation,
)
from .network_structure_analysis import (
    plot_degree_distributions,
    plot_score_discributions,
    plot_network_entropy_distributions,
)

from .gene_analysis import (
    plot_scores_as_rank,
    plot_score_comparison_2D,
    plot_score_comparison_2D_with_plotly,
    plot_score_per_cluster,
    plot_cartography_scatter_per_cluster,
    plot_cartography_term,
)


def load_links(file_path):
    """
    Load links object saved as a hdf5 file.

    Args:
        file_path (str): file path.

    Returns:
        Links: loaded links object.

    """
    links = load_hdf5(filename=file_path, obj_class=Links)

    _update_links_object(links)

    return links


def _update_links_object(links):
    if hasattr(links, "thread_number"):  # if the links object has old version attribute
        if links.thread_number is not None:
            links.threshold_number = deepcopy(links.thread_number)
            delattr(links, "thread_number")


class Links:
    """
    This is a class for the processing and visualization of GRNs.
    Links object stores cluster-specific GRNs and metadata.
    Please use "get_links" function in Oracle object to generate Links object.

    Attributes:
        links_dict (dictionary): Dictionary that store unprocessed network data.
        filtered_links (dictionary): Dictionary that store filtered network data.
        merged_score (pandas.dataframe): Network scores.
        cluster (list of str): List of cluster name.
        name (str): Name of clustering unit.
        palette (pandas.dataframe): DataFrame that store color information.


    """

    def __init__(self, name, links_dict={}):
        """
        Initialize Links object.

        Args:
            name (str): Name for linkLists.
            links_dict (dictionary): Python dictionary that store linkLists.
                Keys are cluster name, values are linkLists.
        """

        self.name = name
        self.links_dict = links_dict.copy()
        self.cluster = list(self.links_dict.keys())
        self.entropy = None

    def to_hdf5(self, file_path):
        """
        Save object as hdf5.

        Args:
            file_path (str): file path to save file. Filename needs to end with '.celloracle.links'
        """
        if file_path.endswith(".celloracle.links"):
            pass
        else:
            raise ValueError("Filename needs to end with '.celloracle.links'")

        compression_opts = 7
        dump_hdf5(
            obj=self,
            filename=file_path,
            data_compression=compression_opts,
            chunks=(2048, 2048),
            noarray_compression=compression_opts,
            pickle_protocol=2,
        )

    def _pipeline(self):
        """
        Perform several processing methods below.
        [self.merge_links(), self.filter_links(), self.get_score()]

        """
        self.filter_links()
        self.get_score()

    def filter_links(
        self,
        p=0.001,
        weight="coef_abs",
        threshold_number=10000,
        genelist_source=None,
        genelist_target=None,
        thread_number=None,
    ):
        """

        Filter network edges.
        In most cases, inferred GRN has non-significant random edges.
        We have to remove these edges before analyzing the network structure.
        You can do the filtering in any of the following ways.

        (1) Filter based on the p-value of the network edge.
            Please enter p-value for thresholding.
        (2) Filter based on network edge number.
            If you set the number, network edges will be filtered based on the order of a network score. The top n-th network edges with network weight will remain, and the other edges will be removed.
            The network data has several types of network weight, so you have to select which network weight do you want to use.
        (3) Filter based on an arbitrary gene list. You can set a gene list for source nodes or target nodes.

        Args:
            p (float): threshold for p-value of the network edge.
            weight (str): Please select network weight name for the filtering
            genelist_source (list of str): gene list to remain in regulatory gene nodes. Default is None.
            genelist_target (list of str): gene list to remain in target gene nodes. Default is None.

        """
        if thread_number is not None:
            threshold_number = thread_number

        self.filtered_links = {}
        self.threshold_number = threshold_number
        for i in self.cluster:
            self.filtered_links[i] = _thresholding(
                linkList=self.links_dict[i],
                p=p,
                weight=weight,
                threshold_number=threshold_number,
                genelist_source=genelist_source,
                genelist_target=genelist_target,
            )

    def get_network_score(self):
        """
         Get several network sores using igraph library.
         The following scores are calculated: ['degree_all', 'degree_centrality_all', 'degree_in',
        'degree_centrality_in', 'degree_out', 'degree_centrality_out',
        'betweenness_centrality', 'eigenvector_centrality']

        """
        if hasattr(self, "filtered_links"):
            if isinstance(self.filtered_links, dict):
                pass
            else:
                raise ValueError(
                    "Filtered network was not found. Prease run 'filter_links' first."
                )
        else:
            raise ValueError(
                "Filtered network was not found. Prease run 'filter_links' first."
            )

        network_scores = []
        for key, val in self.filtered_links.items():
            df = _get_network_score(filtered_linklist_df=val)
            df["cluster"] = key
            network_scores.append(df)
        self.merged_score = pd.concat(network_scores, axis=0)

    def _get_network_score_if_not_exist(self):
        if hasattr(self, "merged_score"):
            if isinstance(self.merged_score, pd.core.frame.DataFrame):
                pass
            else:
                warnings.warn(
                    "\nNetwork score was not found. Celloracle will try to calculate network score first.",
                    UserWarning,
                )
                self.get_network_score()
        else:
            warnings.warn(
                "\nNetwork score was not found. Celloracle will try to calculate network score first.",
                UserWarning,
            )
            self.get_network_score()

    def get_score(self, test_mode=False, n_jobs=-1):
        """
        Get several network sores using R-igraph, linkcomm, and rnetcarto.
        This require R packages.
        """
        warnings.warn(
            "This is the deprecated function. Please use 'get_network_score' instead. This function will be removed in the future version.",
            DeprecationWarning,
        )

        if _check_R_libraries_installation():
            pass
        else:
            raise ValueError(
                "This function is deprecated. Please use 'get_network_score' instead. If you still use this function, you need to install R packages. Please make sure these R packages are installed: 'igraph', 'linkcomm', 'rnetcarto'."
            )

        li = list(self.filtered_links.keys())  # make list of cluster name

        # make dictionary. we make unique id for each cluster and use it for temporary file name.
        id_dict = {}
        for id_, i in enumerate(li):
            id_dict[i] = id_

        _get_network_score_by_Rscripts_inparallel(
            dict_links=self.filtered_links,
            id_dict=id_dict,
            output_folder="network_analysis",
            message=False,
            n_parallel=n_jobs,
        )

        network_scores = {}
        for i in li:
            network_scores[i] = _load_network_analysis_results(
                f"./network_analysis/{id_dict[i]}"
            )
        self.merged_score = _merge_df(network_scores)

        if not test_mode:
            os.system(f"rm -r ./network_analysis/")

        # print(f"the scores are saved in ./{self.name}/")

    def get_network_entropy(self, value="coef_abs"):
        """
        Calculate network entropy scores.

        Args:
            value (str): Default is "coef_abs".
        """
        li = []
        for i in self.cluster:
            mat = _link2mat(self.links_dict[i], value)
            df = _getNetworkEntropy(mat)
            df["cluster"] = i
            # df = df.dropna(axis=0)
            li.append(df)
        self.entropy = pd.concat(li, axis=0)

    #####################################################
    #### Methods for network structure visualization ####
    #####################################################
    def plot_degree_distributions(self, plot_model=False, save=None):
        """
        Plot the network degree distributions (the number of edge per gene).
        The network degree will be visualized in both linear scale and log scale.

        Args:
            links (Links): See network_analysis.Links class for detail.
            plot_model (bool): Whether to plot linear approximation line.
            save (str): Folder path to save plots. If the folder does not exist in the path, the function creates the folder.
               Plots will not be saved if [save=None]. Default is None.

        """
        if hasattr(self, "filtered_links"):
            if isinstance(self.filtered_links, dict):
                pass
            else:
                raise ValueError(
                    "Filtered network was not found. Prease run 'filter_links' first."
                )
        else:
            raise ValueError(
                "Filtered network was not found. Prease run 'filter_links' first."
            )

        plot_degree_distributions(links=self, plot_model=plot_model, save=save)

    def plot_score_discributions(self, values=None, method="boxplot", save=None):
        """
        Plot the distribution of network scores.
        An individual data point is a network edge (gene).

        Args:
            links (Links): See Links class for details.
            values (list of str): The list of score to visualize. If it is None, all of the network score will be used.
            method (str): Plotting method. Select either "boxplot" or "barplot".
            save (str): Folder path to save plots. If the folder does not exist in the path, the function creates the folder.
               Plots will not be saved if [save=None]. Default is None.
        """
        self._get_network_score_if_not_exist()
        plot_score_discributions(links=self, values=values, method=method, save=save)

    def plot_network_entropy_distributions(
        self, update_network_entropy=False, save=None
    ):
        """
        Plot the distribution for network entropy.
        See the CellOracle paper for more detail.

        Args:
            links (Links object): See network_analysis.Links class for detail.
            values (list of str): The list of score to visualize. If it is None, all network score (listed above) will be used.
            update_network_entropy (bool): Whether to recalculate network entropy.
            save (str): Folder path to save plots. If the folder does not exist in the path, the function creates the folder.
               Plots will not be saved if [save=None]. Default is None.
        """
        plot_network_entropy_distributions(
            links=self, update_network_entropy=update_network_entropy, save=save
        )

    ###################################
    #### Methods for gene analysis ####
    ###################################
    def plot_scores_as_rank(self, cluster, n_gene=50, save=None):
        """
        Pick up top n-th genes wich high-network scores and make plots.

        Args:
            links (Links): See network_analysis.Links class for detail.
            cluster (str): Cluster name to analyze.
            n_gene (int): Number of genes to plot. Default is 50.
            save (str): Folder path to save plots. If the folder does not exist in the path, the function creates the folder.
               Plots will not be saved if [save=None]. Default is None.
        """
        self._get_network_score_if_not_exist()

        plot_scores_as_rank(links=self, cluster=cluster, n_gene=n_gene, save=save)

    def plot_score_comparison_2D(
        self,
        value,
        cluster1,
        cluster2,
        percentile=99,
        annot_shifts=None,
        save=None,
        plt_show=True,
        interactive=False,
    ):
        """
        Make a scatter plot that compares specific network scores in two groups.

        Args:
            links (Links): See network_analysis.Links class for detail.
            value (srt): The network score type.
            cluster1 (str): Cluster name. Network scores in cluster1 will be visualized in the x-axis.
            cluster2 (str): Cluster name. Network scores in cluster2 will be visualized in the y-axis.
            percentile (float): Genes with a network score above the percentile will be shown with annotation. Default is 99.
            annot_shifts ((float, float)): Annotation visualization setting.
            save (str): Folder path to save plots. If the folder does not exist in the path, the function creates the folder.
               Plots will not be saved if [save=None]. Default is None.
        """
        self._get_network_score_if_not_exist()

        if interactive:
            return plot_score_comparison_2D_with_plotly(
                links=self,
                value=value,
                cluster1=cluster1,
                cluster2=cluster2,
                fillna_with_zero=True,
            )
        else:
            plot_score_comparison_2D(
                links=self,
                value=value,
                cluster1=cluster1,
                cluster2=cluster2,
                percentile=percentile,
                annot_shifts=annot_shifts,
                save=save,
                plt_show=plt_show,
            )

    def plot_score_per_cluster(self, goi, save=None, plt_show=True):
        """
        Plot network score for a gene.
        This function visualizes the network score for a specific gene between clusters to get an insight into the dynamics of the gene.

        Args:
            links (Links): See network_analysis.Links class for detail.
            goi (srt): Gene name.
            save (str): Folder path to save plots. If the folder does not exist in the path, the function creates the folder.
               Plots will not be saved if [save=None]. Default is None.
        """
        self._get_network_score_if_not_exist()

        plot_score_per_cluster(links=self, goi=goi, save=save, plt_show=plt_show)

    def plot_cartography_scatter_per_cluster(
        self,
        gois=None,
        clusters=None,
        scatter=True,
        kde=False,
        auto_gene_annot=False,
        percentile=98,
        args_dot={"n_levels": 105},
        args_line={"c": "gray"},
        args_annot={},
        save=None,
    ):
        """
        Make a gene network cartography plot.
        Please read the original paper describing gene network cartography for more information.
        https://www.nature.com/articles/nature03288

        Args:
            links (Links): See network_analysis.Links class for detail.
            gois (list of srt): List of gene name to highlight.
            clusters (list of str): List of cluster name to analyze. If None, all clusters in Links object will be analyzed.
            scatter (bool): Whether to make a scatter plot.
            auto_gene_annot (bool): Whether to pick up genes to make an annotation.
            percentile (float): Genes with a network score above the percentile will be shown with annotation. Default is 98.
            args_dot (dictionary): Arguments for scatter plot.
            args_line (dictionary): Arguments for lines in cartography plot.
            args_annot (dictionary): Arguments for annotation in plots.
            save (str): Folder path to save plots. If the folder does not exist in the path, the function creates the folder.
               Plots will not be saved if [save=None]. Default is None.

        """
        warnings.warn(
            "This is the deprecated function. This function will be removed in the future version.",
            DeprecationWarning,
        )

        if "role" not in self.merged_score.columns:
            print(
                "Cartography is not calculated yet. Please run 'get_score' first. This function require several R libraries."
            )
            return None

        plot_cartography_scatter_per_cluster(
            links=self,
            gois=gois,
            clusters=clusters,
            scatter=scatter,
            kde=kde,
            auto_gene_annot=auto_gene_annot,
            percentile=percentile,
            args_dot=args_dot,
            args_line=args_line,
            args_annot=args_annot,
            save=save,
        )

    def plot_cartography_term(self, goi, save=None, plt_show=True):
        """
        Plot the gene network cartography term like a heatmap.
        Please read the original paper of gene network cartography for the principle of gene network cartography.
        https://www.nature.com/articles/nature03288

        Args:
            links (Links): See network_analysis.Links class for detail.
            gois (list of srt): List of gene name to highlight.
            save (str): Folder path to save plots. If the folder does not exist in the path, the function creates the folder.
               Plots will not be saved if [save=None]. Default is None.
        """
        warnings.warn(
            "This is the deprecated function. This function will be removed in the future version.",
            DeprecationWarning,
        )

        if "role" not in self.merged_score.columns:
            print(
                "Cartography is not calculated yet. Please run 'get_score' first. This function require several R libraries."
            )
            return None

        plot_cartography_term(links=self, goi=goi, save=save, plt_show=plt_show)


def _link2mat(link, value="coef_abs", fillna=0):
    mat = pd.pivot(data=link, values=[value], index="target", columns="source")
    mat = mat.fillna(fillna)

    return mat


def _getNetworkEntropy(linkMat):
    tmp = linkMat.copy()
    k = (tmp != 0).sum(axis=1)
    tmp = linkMat.copy()
    ent = []
    ent_norm = []
    genes = []
    for i in tmp.index:
        if k[i] > 1:
            en = stats.entropy(tmp.loc[i])
            ent.append(en)
            ent_norm.append(en / np.log(k[i]))
            genes.append(i)

    df = pd.DataFrame({"entropy": ent, "entropy_norm": ent_norm}, index=genes).dropna(
        axis=0
    )
    return df


def _thresholding(
    linkList,
    p=None,
    weight="coef_abs",
    threshold_number=10000,
    genelist_source=None,
    genelist_target=None,
):
    li = linkList.copy()
    if not genelist_source is None:
        li = li[li.source.isin(genelist_source)]
    if not genelist_target is None:
        li = li[li.target.isin(genelist_target)]

    if not p is None:
        li = li[li["p"] <= p]
    li = li.sort_values(weight, ascending=False)

    if not threshold_number is None:
        li = li[:threshold_number]

    # li = li[["source", "target", weight]]

    return li


def _load_network_analysis_results(folder):
    files = os.listdir(folder)
    network_score = pd.read_csv(
        os.path.join(folder, "base_natwork_analysis.csv"), index_col=0
    )
    # overlapping_cluster = pd.read_csv(os.path.join(folder, "overlapping_cluster.csv"), index_col=0)
    # if "overlapping_cluster_GO.csv" in files:
    #   overlapping_cluster_GO = pd.read_csv(os.path.join(folder, "overlapping_cluster_GO.csv"), index_col=0)
    #    return network_score, overlapping_cluster, overlapping_cluster_GO
    # else:
    return network_score  # , overlapping_cluster


def _merge_df(link_dict):
    merged = []
    clusters = link_dict.keys()
    for i in clusters:
        tmp = link_dict[i].copy()
        tmp["cluster"] = i
        merged.append(tmp)
    merged = pd.concat(merged, axis=0)

    return merged


def _get_network_score(filtered_linklist_df):
    """
    This is the function to get network score for each node.
    This function is made to reproduce previous functions made with R-igraph.
    """
    # Make igraph object
    g = Graph.DataFrame(
        filtered_linklist_df[["source", "target"]], directed=True, use_vids=False
    )
    g.es["weight"] = filtered_linklist_df["coef_abs"].values.copy()

    # Make placeholder
    df = g.get_vertex_dataframe()

    # Calculate scores
    for i in ["all", "in", "out"]:
        df[f"degree_{i}"] = g.degree(mode=i)
        df[f"degree_centrality_{i}"] = df[f"degree_{i}"] / (df.shape[0] - 1)
    df["betweenness_centrality"] = g.betweenness(directed=True, weights="weight")
    df["eigenvector_centrality"] = g.eigenvector_centrality(
        directed=False, weights="weight"
    )

    # Use gene name for index
    df = df.set_index("name")
    df.index.name = None

    return df
