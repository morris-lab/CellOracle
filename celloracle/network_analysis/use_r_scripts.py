# -*- coding: utf-8 -*-
'''
This is a series of custom functions for the inferring of GRN from single cell RNA-seq data.

Codes were written by Kenji Kamimoto.


'''

###########################
### 0. Import libralies ###
###########################


# 0.1. libraries for fundamental data science and data processing


import os
import pandas as pd
from joblib import dump, load

from multiprocessing import cpu_count

N_CPU = cpu_count()

from ..utility import exec_process
from ..network_analysis import __path__ as parent_path
#import seaborn as sns


def test_R_libraries_installation():
    """
    CellOracle.network_analysis use some R libraries for network analysis.
    This is a test function to check instalation of necessary R libraries.
    """

    r_libraries = ["igraph", "linkcomm", "rnetcarto"]#,"gProfileR"]
    for i in r_libraries:
        try:
            exec_process(f"Rscript {parent_path[0]}/rscripts_for_network_analysis/test_{i}.R", message=True)
            print(f"checking R library installation: {i} -> OK")
        except:
            print(f"checking_installation: {i} -> NG")
            print(f" R library, {i} is unavailable. Please check installation.")


def _get_network_score_by_Rscripts(linkList, name, output_folder="network_analysis",
             GO=True, message=False):
    folder = os.path.join(output_folder, name)
    os.makedirs(folder, exist_ok=True)

    link_path =folder + "/linkList.csv"
    linkList[["source", "target", "coef_abs"]].to_csv(link_path, index=None)
    if GO:
        command = f"Rscript {parent_path[0]}/rscripts_for_network_analysis/get_newtork_scores.R {folder}"
    else:
        command = f"Rscript {parent_path[0]}/rscripts_for_network_analysis/get_newtork_scores.R {folder} FALSE"

    exec_process(command, message)


def _get_network_score_by_Rscripts_inparallel(dict_links, output_folder="network_analysis",
                      GO=True, message=False, n_parallel=-1):

    if n_parallel == -1:
        n_parallel = N_CPU
    li = list(dict_links.keys())
    N = len(li)

    def internal(li_):
        process_list = []
        for i in li_:
            folder = os.path.join(output_folder, i)
            os.makedirs(folder, exist_ok=True)

            link_path =folder + "/linkList.csv"
            dict_links[i][["source", "target", "coef_abs"]].to_csv(link_path, index=None)

            if GO:
                command = f"Rscript {parent_path[0]}/rscripts_for_network_analysis/get_newtork_scores.R {folder}"
            else:
                command = f"Rscript {parent_path[0]}/rscripts_for_network_analysis/get_newtork_scores.R {folder} FALSE"
            process_list.append(exec_process(command, message=message,
                                             wait_finished=False, return_process=True))

        for process, name in zip(process_list, li_):

            process.wait()
            if process.returncode != 0:
                print(f'{name}: Build process aborts.')
            else:
                print(f'{name}: finished.')

    if N%n_parallel == 0:
        for k in range(N//n_parallel):
            print(f"processing... batch {k+1}/{N//n_parallel}")
            sub_list = li[n_parallel*k:n_parallel*(k+1)]
            internal(sub_list)

    else:
        for k in range(N//n_parallel):
            print(f"processing... batch {k+1}/{(N//n_parallel)+1}")
            sub_list = li[n_parallel*k:n_parallel*(k+1)]
            internal(sub_list)

        print(f"processing... batch {(N//n_parallel)+1}/{(N//n_parallel)+1}")
        sub_list = li[(N//n_parallel)*n_parallel: N]
        internal(sub_list)
