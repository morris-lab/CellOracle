import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


import seaborn as sns


import os, sys, shutil, importlib, glob
from tqdm.notebook import tqdm

from celloracle import motif_analysis as ma
import celloracle as co
print("co version: ", co.__version__)

plt.rcParams['figure.figsize'] = [6, 4.5]
plt.rcParams["savefig.dpi"] = 300

def main():
    test_tutorial_process_peak_data()

def test_tutorial_process_peak_data():
    # Start analysis
    import time
    start = time.time()
    
    # 1. Load data
    print("1. Load data")
    # 1.0 Download demo data
    print()
    from celloracle.utility.data_download_from_web import download_demo_data
    download_demo_data(file="all_peaks.csv")
    download_demo_data(file="cicero_connections.csv")

    # 1.1. Load data
    # Load scATAC-seq peak list.
    peaks = pd.read_csv("all_peaks.csv", index_col=0)
    peaks = peaks.x.values
    print(peaks)

    # Load Cicero coaccessibility scores.
    cicero_connections = pd.read_csv("cicero_connections.csv", index_col=0)
    cicero_connections.head()

    # 2. Annotate TSS
    print("2. Annotate TSS")
    print(ma.SUPPORTED_REF_GENOME)

    ##!! Please make sure to specify the correct reference genome here
    tss_annotated = ma.get_tss_info(peak_str_list=peaks, ref_genome="mm10")
    # Check results
    print(tss_annotated.tail())

    # 3. Integrate TSS info and cicero connections
    print("3. Integrate TSS info and cicero connections")
    integrated = ma.integrate_tss_peak_with_cicero(tss_peak=tss_annotated,
                                               cicero_connections=cicero_connections)
    print(integrated.shape)
    print(integrated.head())

    # 4. Filter peaks
    print("4. Filter peaks")
    peak = integrated[integrated.coaccess >= 0.8]
    peak = peak[["peak_id", "gene_short_name"]].reset_index(drop=True)

    print(peak.shape)
    print(peak.head())

    # 5. Save data
    print("5. Save data")
    peak.to_csv("processed_peak_file.csv")

    # 6. Remove files
    print("6. Remove files")
    for i in ["all_peaks.csv", "cicero_connections.csv", "processed_peak_file.csv"]:
        os.remove(i)
        #pass
    elapsed = time.time() - start

    print("Success")
    print(f"Total time: {elapsed/60:.3g} minutes")


if __name__ == "__main__":
    main()
