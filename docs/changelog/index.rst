.. _changelog:

=========
Changelog
=========

* `0.5.1 <2020-08-4>`
 `-` Add new promoter-TSS reference data for several reference genomes; (1)"Xenopus": ["xenTro2", "xenTro3"],n(2)"Rat": ["rn4", "rn5", "rn6"], (3)"Drosophila": ["dm3", "dm6"], (4)"C.elegans": ["ce6", "ce10"], (5)"Arabidopsis": ["tair10"].

 `-` Add new motif data for several species: "Xenopus", "Rat", "Drosophila", "C.elegans" and "Arabidopsis".


* `0.5.0 <2020-08-3>`
 `-` Add now functions for custom motifs. You can select motifs from several options. Also, we updated our web tutorial to introduce how to load / make a different motif data.

 `-` Change default motifs for S.cerevisiae and Zebrafish.

 `-` Change requirements for dependent package: gimmemotifs and geomepy. Celloracle codes were updated to support new version of gimmemotifs (0.14.4) and genomepy (0.8.4).


* `0.4.2 <2020-07-14>`
 `-` Add promoter-TSS information for Zebrafish reference genome (danRer7, danRer10 and danRer11).

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
