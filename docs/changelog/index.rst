.. _changelog:

=========
Changelog
=========

* `0.18.0 <2024-03-31>`
 `-` We released celloracle version 0.18.0. Refactoring was done in the motif analysis to improve the stability of the function. 
 
 `-` Docker image was updated to the latest version of celloracle. From this version, we use python 3.10 in the docker image.

 `-` Bug fix in the unit test function. 

* `0.17.1 <2024-03-19>`
 `-` We updated the network function to support scikit-learn>=1.2. 

 `-` CellOracle build test pipeline has been updated to check for Python 3.10. We now support Python 3.10.

* `0.16.0 <2024-02-01>`
 `-` We fixed a bug in the systematic simulation function to support the latest update in pandas.

* `0.15.0 <2023-07-25>`
 `-` Refactoring in the Markov simulation function. 

* `0.14.0 <2023-07-08>`
 `-` Refactoring in some functions related to numpy. 

 `-` Bugfix and update in test functions.

* `0.13.0 <2023-07-08>`
 `-` We updated motif analysis function. We can now select the location of the reference genome.


* `0.10.15 <2023-02-02>`
 `-` We fixed a bug in motif data loading function.

* `0.10.14 <2023-01-03>`
 `-` We added a TSS information and promoter base GRN data for Xenopus Laevis reference genome: "Xenopus_laevis_v10.1"., to answer `user request #93 <https://github.com/morris-lab/CellOracle/issues/93>`_ .

 `-` Default Motif for Xenopus tropicalis was updated from "CisBP_ver2_Xenopus_tropicalis_and_Xenopus_laevis.pfm" to "CisBP_ver2_Xenopus_tropicalis.pfm". Because of this motif data update, default promoter base GRN for Xenopus tropicalis is also updated.

* `0.10.12 <2022-10-10>`
 `-` We added a function to check distribution range of simulated value.

* `0.10.11 <2022-09-27>`
 `-` We added a function to update cluster color information. See `celloracle GitHub issue #91 <https://github.com/morris-lab/CellOracle/issues/91>`_ for details.

* `0.10.10 <2022-09-09>`
 `-` Due to a change in the latest the python-igraph package (version 0.10.0), it has been reported that the get_network_score function is failing. This bug has been fixed. As a result, the latest version 0.10.0 of python-igraph has been temporarily made unsupported, since errors are expected to occur in the louvain function of scanpy.

 `-` Bug fix in the network score visualization plot to deal with recent update in seaborn package.

* `0.10.5 <2022-05-08>`
 `-` We refactored the data loading function. In the version with 0.10.5 or higher, the data can be automatically downloaded from GitHub if the data is not found in the local celloracle library directory.

 `-` We added `Travis CI page <https://app.travis-ci.com/github/morris-lab/CellOracle>`_ to automatically check celloracle build status.

 `-` We updated CellOracle GitHub repository README file to add celloracle package status quality check badges.

* `0.10.1 <2022-04-13>`
 `-` We refactored the code to simplify the installation process and calculation speed.

 `-` Links.get_score() function was deprecated. Please use Links.get_network_score() instead. Because the new network score calculation is implemented solely in Python code, celloracle no longer requires the R package.

 `-` Codes were updated to comply with the recent numpy updates. e.g. np.int was replaced with int.

* `0.10.0 <2022-04-13>`
 `-` Tutorial data was updated.

* `0.9.0 <2022-01-21>`
 `-` Change the default color scheme for the inner product score / perturbation score visualization in the in silico gene perturbation module.

* `0.8.5 <2022-01-11>`
 `-` Bug fix on plot saving function in the Network analysis module. The previous version may generate an error when the user try to save plots in network analysis. We fixed the bug.

* `0.8.4 <2021-12-29>`
 `-` Code refactoring in simulation results visualization function.

 `-` Code refactoring in motif analysis module.

 `-` Add TSS data for Guinea Pig reference genome, Cavpor3.0.


* `0.8.3 <2021-11-25>`
 `-` Fix typo in the Markov walk function.

* `0.8.2 <2021-10-31>`
 `-` Code refactoring in motif analysis module.

 * `0.8.1 <2021-10-30>`
  `-` Change requirements.

* `0.8.0 <2021-08-28>`
 `-` Change requirements. From this version, numba>=0.50.1 is required.

 `-` Update installation page in documentation.

* `0.7.5 <2021-07-28>`
 `-` Correct requirements.txt file name.

* `0.7.4 <2021-07-27>`
 `-` Update Arabidopsis promoter base GRN data.

* `0.7.3 <2021-07-25>`
 `-` Update Arabidopsis motif data.

* `0.7.0 <2021-07-16>`
 `-` Overhaul documentation.

* `0.7.0 <2021-07-11>`
 `-` Add pre-built promoter base GRNs.

* `0.7.1 <2021-07-15>`
 `-` Aad function for oracle transition probability calculation.

* `0.7.0 <2021-07-11>`
 `-` Add pre-built promoter base GRNs.

* `0.6.17 <2021-07-08>`
 `-` Add chicken and guinea pig motif

 `-` Update Arabidopsis ref genome name

* `0.6.12 <2021-06-11>`
 `-` Add functions to oracle object to check current data status.

* `0.6.11 <2021-06-09>`
 `-` Add data loading function. Demo oracle data and links data can be loaded using data loading functions.

* `0.6.9 <2021-05-14>`
 `-` Code refactoring in network visualization.

* `0.6.8 <2021-05-10>`
 `-` Update Seurat data conversion module.

* `0.6.8 <2021-05-08>`
 `-` Change requirements. From this version, numba=0.48.0  is required.

* `0.6.7 <2021-05-5>`
 `-` Add function to check status of installed dependent package version.

* `0.6.5 <2021-03-25>`
 `-` Minor bug fix in the installation process.

* `0.6.4 <2021-02-18>`
 `-` Minor change for oracle object. Metadata will be shown if you print oracle object.

 `-` Add new function to oracle class.

* `0.6.3 <2021-01-26>`
 `-` Big fix to solve [this issue](https://github.com/morris-lab/CellOracle/issues/42).

 `-` Bug fix. Anndata>=0.7.5 is required.

* `0.6.2 <2021-12-16>`
 `-` Big fix. h5py>=3.1.0 is required.

* `0.6.0 <2021-12-14>`
 `-` Add new modules: dev_modules and analysis_helper.

* `0.5.1 <2020-08-4>`
 `-` Add new promoter-TSS reference data for several reference genomes; (1)"Xenopus": ["xenTro2", "xenTro3"], (2)"Rat": ["rn4", "rn5", "rn6"], (3)"Drosophila": ["dm3", "dm6"], (4)"C.elegans": ["ce6", "ce10"], (5)"Arabidopsis": ["tair10"].

 `-` Add new motif data for several species: "Xenopus", "Rat", "Drosophila", "C.elegans" and "Arabidopsis".

* `0.5.0 <2020-08-3>`
 `-` Add new functions for custom motifs. You can select motifs from several options. Also, we updated our web tutorial to introduce how to load / make different motif data.

 `-` Change default motifs for S.cerevisiae and zebrafish.

 `-` Change requirements for dependent package: gimmemotifs and geomepy. Celloracle codes were updated to support new version of gimmemotifs (0.14.4) and genomepy (0.8.4).


* `0.4.2 <2020-07-14>`
 `-` Add promoter-TSS information for zebrafish reference genome (danRer7, danRer10 and danRer11).

* `0.4.1 <2020-07-02>`
 `-` Add promoter-TSS information for S.cerevisiae reference genome (sacCer2 and sacCer3).

* `0.4.0 <2020-06-28>`
 `-` Change requirements.

 `-` From this version, pandas version 1.0.3 or later is required.

 `-` From this version, scanpy version 1.5.3 or later is required.

* `0.3.7 <2020-06-12>`
 `-` Delete GO function from r-script

 `-` Update some functions for network visualization

* `0.3.6 <2020-06-08>`
 `-` Fix a bug on the transition probability calculation in Markov simulation

 `-` Add new function "count_cells_in_mc_results" to oracle class

* `0.3.5 <2020-05-09>`
 `-` Fix a bug on the function for gene cortography visualization

 `-` Change some settings for installation

 `-` Update data conversion module

* `0.3.4 <2020-04-29>`
 `-` Change pandas version restriction

 `-` Fix a bug on the function for gene cortography visualization

 `-` Add new functions for R-path configuration

* `0.3.3 <2020-04-24>`
 `-` Add promoter-TSS information for hg19 and hg38 reference genome

* `0.3.1 <2020-03-23>`
 `-` Fix an error when try to save file larger than 4GB file

* `0.3.0 <2020-2-17>`
 `-` Release beta version
