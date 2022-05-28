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
    genome_installation = ma.is_genome_installed(ref_genome=ref_genome)
    print(ref_genome, "installation: ", genome_installation)
    ## 1.2. Install reference genome (if refgenome is not installed)
    if not genome_installation:
        import genomepy
        genomepy.install_genome(ref_genome, "UCSC")
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


    ## 2.2. Check data
    def decompose_chrstr(peak_str):
        """
        Args:
            peak_str (str): peak_str. e.g. 'chr1_3094484_3095479'

        Returns:
            tuple: chromosome name, start position, end position
        """

        *chr_, start, end = peak_str.split("_")
        chr_ = "_".join(chr_)
        return chr_, start, end

    from genomepy import Genome

    def check_peak_format(peaks_df, ref_genome):
        """
        Check peak format.
         (1) Check chromosome name.
         (2) Check peak size (length) and remove sort DNA sequences (<5bp)

        """

        df = peaks_df.copy()

        n_peaks_before = df.shape[0]

        # Decompose peaks and make df
        decomposed = [decompose_chrstr(peak_str) for peak_str in df["peak_id"]]
        df_decomposed = pd.DataFrame(np.array(decomposed), index=peaks_df.index)
        df_decomposed.columns = ["chr", "start", "end"]
        df_decomposed["start"] = df_decomposed["start"].astype(int)
        df_decomposed["end"] = df_decomposed["end"].astype(int)

        # Load genome data
        genome_data = Genome(ref_genome)
        all_chr_list = list(genome_data.keys())


        # DNA length check
        lengths = np.abs(df_decomposed["end"] - df_decomposed["start"])


        # Filter peaks with invalid chromosome name
        n_threshold = 5
        df = df[(lengths >= n_threshold) & df_decomposed.chr.isin(all_chr_list)]

        # DNA length check
        lengths = np.abs(df_decomposed["end"] - df_decomposed["start"])

        # Data counting
        n_invalid_length = len(lengths[lengths < n_threshold])
        n_peaks_invalid_chr = n_peaks_before - df_decomposed.chr.isin(all_chr_list).sum()
        n_peaks_after = df.shape[0]


        #
        print("Peaks before filtering: ", n_peaks_before)
        print("Peaks with invalid chr_name: ", n_peaks_invalid_chr)
        print("Peaks with invalid length: ", n_invalid_length)
        print("Peaks after filtering: ", n_peaks_after)

        return df


    peaks = check_peak_format(peaks, ref_genome)


    # 3. Instantiate TFinfo object and search for TF binding motifs
    print("\n3. Instantiate TFinfo object and search for TF binding motifs")
    ## 3.1. Instantiate TFinfo object

    if use_small_data:
        peaks = peaks[:100]
    tfi = ma.TFinfo(peak_data_frame=peaks,
                    ref_genome=ref_genome)
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
