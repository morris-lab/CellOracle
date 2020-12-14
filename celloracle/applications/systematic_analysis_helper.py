# -*- coding: utf-8 -*-



import io
import logging
import os
import pickle
import subprocess
import sys
from copy import deepcopy

import pandas as pd
import numpy as np
from tqdm.notebook import tqdm
import matplotlib.pyplot as plt
#import h5py


from IPython.display import display, HTML
from ipywidgets import interactive, FloatText

from .development_module import Oracle_development_module
from ..visualizations.config import CONFIG


class Oracle_systematic_analysis_helper(Oracle_development_module):

    def __init__(self, hdf5_file_path, memory_save_mode=False):
        super().__init__()

        self.set_hdf_path(path=hdf5_file_path)
        self.hdf5_info = self.get_hdf5_info()
        self.estimated_scales_for_visulization = None
        self.negative_ip_sum = None
        self.positive_ip_sum = None
        self.memory_save_mode = memory_save_mode
        self._exemptions_when_del_attrs = \
            ["hdf5_info", "negative_ip_sum", "positive_ip_sum", "estimated_scales_for_visulization"]

    def get_negative_ip_sum_for_all_data(self, verbose=True, return_result=True):

        self.del_attrs()

        gene_misc_lists = self.hdf5_info["gene_misc_lists"]

        li = []
        if verbose:
            loop = tqdm(gene_misc_lists)
        else:
            loop = gene_misc_lists

        for gene, misc in loop:

            self.load_hdf5(gene=gene, misc=misc, specify_attributes=["inner_product_df"])
            df = self.get_sum_of_negative_ips()
            df["gene"] = gene
            df["misc"] = misc
            li.append(df)

        negative_ip_sum = pd.concat(li, axis=0).reset_index(drop=True)

        self.negative_ip_sum = negative_ip_sum

        # Clear memory
        self.del_attrs()

        if return_result:
            return self.negative_ip_sum


    def get_positive_ip_sum_for_all_data(self, verbose=True, return_result=True):

        self.del_attrs()

        gene_misc_lists = self.hdf5_info["gene_misc_lists"]

        li = []
        if verbose:
            loop = tqdm(gene_misc_lists)
        else:
            loop = gene_misc_lists

        self.del_attrs()

        for gene, misc in loop:
            self.load_hdf5(gene=gene, misc=misc, specify_attributes=["inner_product_df"])
            df = self.get_sum_of_positive_ips()
            df["gene"] = gene
            df["misc"] = misc
            li.append(df)

        positive_ip_sum = pd.concat(li, axis=0).reset_index(drop=True)

        self.positive_ip_sum = positive_ip_sum

        # Clear memory
        self.del_attrs()

        if return_result:
            return positive_ip_sum

    def sort_TFs_by_neagative_ip(self, misc, pseudotime="0,1,2,3,4,5,6,7,8,9"):
        if self.negative_ip_sum is None:
            self.get_negative_ip_sum_for_all_data(return_result=False, verbose=False)

        negative_ip_sum = self.negative_ip_sum

        # Focus on a misc
        df = negative_ip_sum[negative_ip_sum.misc == misc]

        # Focus on a specific pseudotime
        pseudotime = [i for i in pseudotime if i in list("0123456789")]
        df = df[df.pseudotime_id.isin(pseudotime)]

        # Get sum of negative ip values
        df = df[["gene", "score"]].groupby(["gene"]).sum()

        # Sore by sum score
        df = df.sort_values("score", ascending=True).reset_index(drop=False)

        return df

    def interactive_sort_TFs_by_neagative_ip(self):

        if self.hdf5_info is None:
            self.get_hdf5_info()

        def wrapper(misc, pseudotime, n_TFs=20):
            df = self.sort_TFs_by_neagative_ip(misc=misc, pseudotime=pseudotime)

            print(f"Top {n_TFs} in {misc}")
            display(HTML(df.iloc[:min(n_TFs, df.shape[0])].to_html()))

        interactive_table = interactive(wrapper,
                               #{'manual': True},
                               misc=self.hdf5_info["misc_list"],
                               pseudotime="0,1,2,3,4,5,6,7,8,9",
                               n_TFs=(5, 50, 1))

        return interactive_table


    def estimate_scale_for_visualization(self, return_result=True):

        if self.hdf5_info is None:
            self.get_hdf5_info()

        misc_list = self.hdf5_info["misc_list"]

        scales = np.array([self._estimate_scale(misc=misc) for misc in misc_list])
        scales = pd.DataFrame(scales, index=misc_list, columns=["vm", "simulation", "pseudotime"])

        self.estimated_scales_for_visulization = scales

        if return_result:
            return scales


    def _estimate_scale(self, misc):

        self.del_attrs()

        reference_scaling_dictionary = {"inner_product": 0.5,
                                        "flow": 5,
                                        "ref_flow": 18}

        scales = []
        for i in ["inner_product", "flow", "ref_flow"]:
            max_ = self.__check_max_of_max_norm_attr(misc=misc, attr_name=i)
            scales.append(max_ * reference_scaling_dictionary[i])

        # Clear memory
        self.del_attrs()

        # this is the lis of scales
        return scales

    def __check_max_of_max_norm_attr(self, misc, attr_name):

        """
        (1) Load attr and check max norm value.
        (2) Repeat step (1) for all genes in the misc.
        (3) Return max of all max of norm value list.

        if attr is 2D value, calculate max of l2 norm.

        """

        if self.hdf5_info is None:
            self.get_hdf5_info()

        gene_list = self.hdf5_info["misc_gene_dictionary"][misc]


        if attr_name == "ref_flow": # Don't need to check all gene because ref_flow is unique to misc.
            max_ = self.__check_max_norm_attr(misc=misc, gene=gene_list[0], attr_name=attr_name)

        else:
            maxs = []
            for gene in gene_list:
                abs_max = self.__check_max_norm_attr(misc=misc, gene=gene, attr_name=attr_name)
                maxs.append(abs_max)
            max_ = max(maxs)

        return max_

    ### For interactive systematic visualization
    def __check_max_norm_attr(self, misc, gene, attr_name):

        """
        Load attr and check max norm value.
        if attr is 2D value, calculate max of l2 norm.

        """
        self.load_hdf5(gene=gene, misc=misc, specify_attributes=[attr_name])
        attr = getattr(self, attr_name)

        if attr.ndim == 1:
            abs_max = np.absolute(attr).max()
        elif attr.ndim == 2:
            abs_max = np.linalg.norm(attr, axis=1).max()

        return abs_max

    def interactive_visualize_layout_0(self, scale_simulation=None, scale_pseudotime=None, vm=None, s=5, s_grid=30):

        # 1. Define a custom function
        def wrapper(gene, misc, scale_simulation, scale_pseudotime, vm, background):
            # Load data
            #self.del_attrs()

            self.load_hdf5(gene=gene, misc=misc)
            # Visualize
            self.visualize_development_module_layout_0(s=s,
                                                       scale_for_simulation=scale_simulation,
                                                       s_grid=s_grid,
                                                       scale_for_pseudotime=scale_pseudotime,
                                                       vm=vm,
                                                       show_background=background)
            print("Gene: ", gene)
            # Delete loaded data
            #self.del_attrs(exemptions=self._exemptions_when_del_attrs)


        # 2. Parameter setting
        # Use automatically estimated value for default scale parameters if it is not specified.

        if np.any([i is None for i in [scale_simulation, scale_pseudotime, vm]]):
            if self.estimated_scales_for_visulization is None:
                self.estimate_scale_for_visualization(return_result=False)
            scale = self.estimated_scales_for_visulization.max(axis=0)
            scale_simulation, scale_pseudotime, vm = scale[["simulation", "pseudotime", "vm"]]
            #scale_simulation, vm = [i * 0.5 for i in [scale_simulation, vm]]
            scale_simulation, scale_pseudotime, vm = [np.round(i, 2) for i in [scale_simulation, scale_pseudotime, vm]]

        ui_gene_list = self.hdf5_info["gene_list"]
        ui_misc_list = self.hdf5_info["misc_list"]

        interactive_viz = interactive(wrapper,
                                       #{'manual': True},
                                  scale_simulation=FloatText(value=scale_simulation),
                                  scale_pseudotime=FloatText(value=scale_pseudotime),
                                  vm=FloatText(value=vm),
                                  gene=ui_gene_list,
                                  misc=ui_misc_list,
                                  background=[True, False])

        return interactive_viz


    def get_all_ips(self, misc):

        misc = "Whole_cells"
        genes = self.hdf5_info["misc_gene_dictionary"][misc]

        # Load inner product for all genes
        ips = []
        for gene in genes:
            self.load_hdf5(gene=gene, misc=misc, specify_attributes=["inner_product"])
            ips.append(helper.inner_product)
        ips = np.stack(ips)
        ips = pd.DataFrame(ips, index=genes)

        # Clear memory
        self.del_attrs()

        return ips

    def get_corrcoef_ip(self, misc):

        ips = self.get_all_ips(misc=misc)

        # Calculate correlation
        corrcoef = np.corrcoef(ips)

        return corrcoef
