# -*- coding: utf-8 -*-
'''
This file contains custom functions for the analysis of ATAC-seq data.
Genomic activity information (peak of ATAC-seq) will be extracted first.
Then the peak DNA sequence will be subjected to TF motif scan.
Finally we will get list of TFs that potentially binds to a specific gene.

Codes were written by Kenji Kamimoto.


'''

###########################
### 0. Import libralies ###
###########################


# 0.1. libraries for fundamental data science and data processing

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from tqdm import tqdm_notebook as tqdm

##
from ..utility import exec_process
from .. import network as tn

#########################
### prepare chip data ###
#########################


@np.vectorize
def intersect_chip_bed_with_atac_bed(base_bed_file_path, chip_bed_file_path):

    file_name = chip_bed_file_path.split("/")[-1]
    GEO_name = file_name.split("_")[0]
    output_file_path = os.path.dirname(chip_bed_file_path) + f"/{GEO_name}_intersected.bed"

    command = f"bedtools intersect -a {base_bed_file_path} -b {chip_bed_file_path} -wa > {output_file_path}"
    exec_process(commands=command, message=True, wait_finished=True)

    print(f"processed: {file_name}")
    return output_file_path
intersect_chip_bed_with_atac_bed.excluded.add(0)


def load_and_merge_chip_bed_files(base_bed_file_path,
                                  chip_bed_file_paths,
                                  chip_type,
                                  save_path=None):


    print("1/3 loading base atac bed file...")
    merged = _read_peaks(base_bed_file_path)

    print("2/3 processing chip files...")
    for chip_bed_file_path in tqdm(chip_bed_file_paths):


        file_name = chip_bed_file_path.split("/")[-1]
        GEO_name = file_name.split("_")[0]
        sample_name = f"{GEO_name}_{chip_type}"
        print(f" processing {file_name}...")

        chip_bed_df = _read_peaks(chip_bed_file_path)
        chip_bed_df[sample_name] = 1
        if chip_bed_df.shape[0] == 0:
            print(" empty data")
        else:
            merged = pd.concat([merged, chip_bed_df[[sample_name]]],
                           axis=1, sort=False)

    print("3/3 merging files...")
    merged = merged.fillna(value=0)
    merged = merged.drop(["peak_chr", "peak_start", "peak_end"], axis=1)

    if save_path is None:
        save_path = os.path.dirname(chip_bed_file_paths[0]) + "/processed_chip_data.parquet"

    merged.to_parquet(save_path)
    print(f"processed file was saved in \n{save_path}")
    return merged

def _read_peaks(path):
    try:
        tmp = pd.read_csv(path, header=None, delimiter="\t")

    except Exception as e:
        print(e)
        tmp = pd.DataFrame(columns=["peak_chr", "peak_start", "peak_end"])
        return tmp

    tmp.columns = ["peak_chr", "peak_start", "peak_end"]
    ind = [tmp.iloc[i]["peak_chr"] + "_" +
           str(tmp.iloc[i]["peak_start"]) + "_" +
           str(tmp.iloc[i]["peak_end"]) for i in range(len(tmp))]
    tmp.index = ind

    # remove duplicated values
    tmp = tmp[~tmp.index.duplicated(keep="first")]

    return tmp

#####################################################
### integrate Chip data and GRN inference results ###
###  fucusing on a TF                             ###
#####################################################

def integrate_Chip_and_TN(tn_object, Chip_data, TF_of_interest,
                          mode="gene_mode"):

    # version190605


    # extract grn data
    grnDF_mean = tn.getDF_peakxTF(tn_object, "coef_abs")
    grnDF_ps = tn.getDF_peakxTF(tn_object, "-logp")

    if not TF_of_interest in grnDF_mean.columns:
        grnDF_mean[TF_of_interest] = 0
        grnDF_ps[TF_of_interest] = 0

    grnDF_mean = grnDF_mean[["gene_short_name", TF_of_interest]]
    grnDF_ps = grnDF_ps[[TF_of_interest]]
    grnDF = pd.concat([grnDF_mean, grnDF_ps], axis=1)

    grnDF.columns = ["gene_short_name", "coef_abs", "-logp"]

    if mode == "gene_mode":
        grnDF = grnDF.groupby("gene_short_name").max()
        promoter_enhancer_peak = None
    if mode == "peak_mode":
        grnDF = grnDF.drop("gene_short_name", axis=1)
        promoter_enhancer_peak = grnDF.index.values

    Chip_data = Chip_data.reindex(grnDF.index.values).fillna(0)

    integrated = pd.concat([Chip_data, grnDF], axis=1)



    ## change data structure

    chip_sample_names = list(Chip_data.columns)
    chip_sample_names = [i for i in chip_sample_names if not "gene" in i]

    li = []
    for i in chip_sample_names:

        tmp = integrated[[i, "coef_abs", "-logp"]]
        tmp = tmp.sort_values(by=i, ascending=False)
        tmp.index.name = "gene_short_name"
        tmp = tmp.reset_index(drop=False)
        tmp["chip_score_order"] = tmp.index.values

        tmp.index = tmp["gene_short_name"]
        tmp = tmp.drop('gene_short_name', axis=1)

        tmp["sample_id"] = i

        tmp = tmp.rename(columns={i : "chip_bind"})
        li.append(tmp)


    integrated = pd.concat(li, axis=0)


    return integrated, promoter_enhancer_peak


def loadAndProcessData(save_dir, TF_of_interest, chip_data, tn_paths):
    os.makedirs(save_dir, exist_ok=True)
    for tn_path in tqdm(tn_paths):
        # define save path
        save_path = os.path.join(save_dir,
                                 tn_path.split("/")[-1] + "_result.parquet")
        # load TransNetObject
        tn1 = tn.load_compressed_TN_object(tn_path)
        # integrate Chip seq data with TN results
        integrated, _ = integrate_Chip_and_TN(tn1, chip_data, TF_of_interest)
        # save integrated object
        integrated.to_parquet(os.path.join(save_path))

#####################################################
### integrate Chip data and GRN inference results ###
###  fucusing on a Histone Marker                 ###
#####################################################

def integrate_Chip_and_TN_for_HM_analysus(tn_object, Chip_data):

    # version190607

    # extract grn data
    grnDF = tn.getDF_peakxTF(tn_object, "coef_abs")

    # get sum of absolute coef value.
    # this values sppose to represent activity score for this gene locus
    grnDF = pd.DataFrame(grnDF.drop("gene_short_name", axis=1).abs().sum(axis=1),
                         columns=["coef_abs_sum"])
    #return grnDF
    grnDF = grnDF.groupby("peak_id").max()
    grnDF.index.name = None

    promoter_enhancer_peak = grnDF.index.values

    # integragte gene activity score with Chip data
    integrated = pd.concat([Chip_data.reindex(grnDF.index.values).fillna(0),
                            grnDF], axis=1)
    #return integrated
    # change data format
    integrated = integrated.reset_index(drop=False)
    integrated = integrated.melt(id_vars=["index", "coef_abs_sum"])
    integrated = integrated.set_index("index")
    integrated = integrated.rename(columns={"variable": "sample_id",
                                            "value": "chip_bind"})
    return integrated, promoter_enhancer_peak



def loadAndProcessData_for_HM_analysis(
        save_dir, chip_data, tn_paths):
    os.makedirs(save_dir, exist_ok=True)
    for tn_path in tqdm(tn_paths):
        # define save path
        save_path = os.path.join(save_dir, tn_path.split("/")[-1] + "_result.parquet")
        # load TransNetObject
        tn1 = tn.load_net_from_patquets(tn_path)
        # integrate Chip seq data with TN results
        integrated, _ = integrate_Chip_and_TN_for_HM_analysus(tn1, chip_data)
        # save integrated object
        integrated.to_parquet(os.path.join(save_path))
