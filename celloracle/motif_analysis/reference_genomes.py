# -*- coding: utf-8 -*-
'''

'''


import pandas as pd

SUPPORTED_REF_GENOME = \
    pd.DataFrame([["Human", "hg38", "UCSC"],
                  ["Human", "hg19", "UCSC"],
                  ["Mouse", 'mm39', "UCSC"],
                  ["Mouse", 'mm10', "UCSC"],
                  ["Mouse", 'mm9', "UCSC"],
                  ["S.cerevisiae", "sacCer2", "UCSC"],
                  ["S.cerevisiae", "sacCer3", "UCSC"],
                  ["Zebrafish", "danRer7", "UCSC"],
                  ["Zebrafish", "danRer10", "UCSC"],
                  ["Zebrafish", "danRer11", "UCSC"],
                  ["Xenopus", "xenTro2", "UCSC"],
                  ["Xenopus", "xenTro3", "UCSC"],
                  ["Rat", "rn4", "UCSC"],
                  ["Rat", "rn5", "UCSC"],
                  ["Rat", "rn6", "UCSC"],
                  ["Drosophila", "dm3", "UCSC"],
                  ["Drosophila", "dm6", "UCSC"],
                  ["C.elegans", "ce6", "UCSC"],
                  ["C.elegans", "ce10", "UCSC"],
                  ["Arabidopsis", "TAIR10", "Ensembl"],
                  ["Chicken", "galGal4", "UCSC"],
                  ["Chicken", "galGal5", "UCSC"],
                  ["Chicken", "galGal6", "UCSC"],
                  ["Guinea_Pig", "Cavpor3.0", "Ensembl"],
                  ["Axolotl", "AmexG_v6.0-DD", "Axolotl-omics.org"]
                  ],
                 columns=["species", "ref_genome", "provider"])
