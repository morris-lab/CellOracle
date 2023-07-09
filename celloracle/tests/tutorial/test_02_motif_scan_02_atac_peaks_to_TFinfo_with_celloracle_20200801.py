import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


import seaborn as sns

import os, sys, shutil, importlib, glob
from tqdm.notebook import tqdm

import celloracle as co
from celloracle import motif_analysis as ma
from celloracle.utility import save_as_pickled_object
print("celloracle version: ", co.__version__)

plt.rcParams['figure.figsize'] = (15,7)
plt.rcParams["savefig.dpi"] = 600

def main():
    test_tutorial_motif_scan()

def test_tutorial_motif_scan():

    use_small_data = True

    # Start analysis
    import time
    start = time.time()


    # 1. Reference genome data preparation
    print("\n1. Reference genome data preparation")
    # 1.1 Check referenoce genome installation
    ref_genome = "mm10"
    genome_installation = ma.is_genome_installed(ref_genome=ref_genome, genomes_dir=None)
    print(ref_genome, "installation: ", genome_installation)
    ## 1.2. Install reference genome (if refgenome is not installed)
    if not genome_installation:
        import genomepy
        genomepy.install_genome(ref_genome, "UCSC", genomes_dir=None)
    else:
        print(ref_genome, "is installed.")


    # 2. Load data
    # 2.0 Download demo data
    print("\n2. Load data")
    from celloracle.utility.data_download_from_web import download_demo_data
    for i in ["processed_peak_file.csv"]:
        download_demo_data(file=i)

    ## 2.1. Load processed peak data
    # Load annotated peak data.
    peaks = pd.read_csv("processed_peak_file.csv", index_col=0)
    peaks.head()


    peaks = ma.check_peak_format(peaks_df=peaks, ref_genome=ref_genome, genomes_dir=None)


    # 3. Instantiate TFinfo object and search for TF binding motifs
    print("\n3. Instantiate TFinfo object and search for TF binding motifs")
    ## 3.1. Instantiate TFinfo object

    if use_small_data:
        peaks = peaks[:100]
    tfi = ma.TFinfo(peak_data_frame=peaks,
                    ref_genome=ref_genome,
                    genomes_dir=None)
    ## 3.2. Motif scan
    tfi.scan(fpr=0.02,
             motifs=None,  # If you enter None, default motifs will be loaded.
             verbose=True)

    # Save tfinfo object
    tfi.to_hdf5(file_path="test1.celloracle.tfinfo")
    # Check motif scan results
    tfi.scanned_df.head()


    # 4. Filtering motifs
    print("\n4. Filtering motifs")
    # Reset filtering
    tfi.reset_filtering()
    # Do filtering
    tfi.filter_motifs_by_score(threshold=10)
    # Format post-filtering results.
    tfi.make_TFinfo_dataframe_and_dictionary(verbose=True)


    # 5. Get final base GRN
    print("\n5. Get final base GRN")
    df = tfi.to_dataframe()
    df.head()

    # 6. Save result as a dataframe
    print("\n6. Save result as a dataframe")
    df = tfi.to_dataframe()
    df.to_parquet("base_GRN_dataframe.parquet")


    # Remove files
    print("Remove files")
    for i in ["processed_peak_file.csv", "test1.celloracle.tfinfo"]:
        os.remove(i)
        #pass
    elapsed = time.time() - start

    print("Success")
    print(f"Total time: {elapsed/60:.3g} minutes")


if __name__ == "__main__":
    main()
