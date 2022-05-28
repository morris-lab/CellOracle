# -*- coding: utf-8 -*-
import pandas as pd

import celloracle as co

# CisBP_ver2_Arabidopsis_thaliana
def test_load_motifs_CisBP_ver2_Arabidopsis_thaliana():
    motifs_name = 'CisBP_ver2_Arabidopsis_thaliana.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=False)
    assert isinstance(motifs, list)

def test_load_motifs_dl_CisBP_ver2_Arabidopsis_thaliana():
    motifs_name = 'CisBP_ver2_Arabidopsis_thaliana.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=True)
    assert isinstance(motifs, list)

# CisBP_ver2_Arabidopsis_thaliana_GENE_ID
def test_load_motifs_CisBP_ver2_Arabidopsis_thaliana_GENE_ID():
    motifs_name = 'CisBP_ver2_Arabidopsis_thaliana_GENE_ID.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=False)
    assert isinstance(motifs, list)

def test_load_motifs_dl_CisBP_ver2_Arabidopsis_thaliana_GENE_ID():
    motifs_name = 'CisBP_ver2_Arabidopsis_thaliana_GENE_ID.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=True)
    assert isinstance(motifs, list)

# CisBP_ver2_Caenorhabditis_elegans
def test_load_motifs_CisBP_ver2_Caenorhabditis_elegans():
    motifs_name = 'CisBP_ver2_Caenorhabditis_elegans.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=False)
    assert isinstance(motifs, list)

def test_load_motifs_dl_CisBP_ver2_Caenorhabditis_elegans():
    motifs_name = 'CisBP_ver2_Caenorhabditis_elegans.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=True)
    assert isinstance(motifs, list)

# CisBP_ver2_Cavia_porcellus
def test_load_motifs_CisBP_ver2_Cavia_porcellus():
    motifs_name = 'CisBP_ver2_Cavia_porcellus.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=False)
    assert isinstance(motifs, list)

def test_load_motifs_dl_CisBP_ver2_Cavia_porcellus():
    motifs_name = 'CisBP_ver2_Cavia_porcellus.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=True)
    assert isinstance(motifs, list)

# CisBP_ver2_Danio_rerio
def test_load_motifs_CisBP_ver2_Danio_rerio():
    motifs_name = 'CisBP_ver2_Danio_rerio.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=False)
    assert isinstance(motifs, list)

def test_load_motifs_dl_CisBP_ver2_Danio_rerio():
    motifs_name = 'CisBP_ver2_Danio_rerio.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=True)
    assert isinstance(motifs, list)

# CisBP_ver2_Drosophila_ananassae
def test_load_motifs_CisBP_ver2_Drosophila_ananassae():
    motifs_name = 'CisBP_ver2_Drosophila_ananassae.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=False)
    assert isinstance(motifs, list)

def test_load_motifs_dl_CisBP_ver2_Drosophila_ananassae():
    motifs_name = 'CisBP_ver2_Drosophila_ananassae.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=True)
    assert isinstance(motifs, list)

# CisBP_ver2_Drosophila_erecta
def test_load_motifs_CisBP_ver2_Drosophila_erecta():
    motifs_name = 'CisBP_ver2_Drosophila_erecta.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=False)
    assert isinstance(motifs, list)

def test_load_motifs_dl_CisBP_ver2_Drosophila_erecta():
    motifs_name = 'CisBP_ver2_Drosophila_erecta.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=True)
    assert isinstance(motifs, list)

# CisBP_ver2_Drosophila_grimshawi
def test_load_motifs_CisBP_ver2_Drosophila_grimshawi():
    motifs_name = 'CisBP_ver2_Drosophila_grimshawi.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=False)
    assert isinstance(motifs, list)

def test_load_motifs_dl_CisBP_ver2_Drosophila_grimshawi():
    motifs_name = 'CisBP_ver2_Drosophila_grimshawi.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=True)
    assert isinstance(motifs, list)

# CisBP_ver2_Drosophila_melanogaster
def test_load_motifs_CisBP_ver2_Drosophila_melanogaster():
    motifs_name = 'CisBP_ver2_Drosophila_melanogaster.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=False)
    assert isinstance(motifs, list)

def test_load_motifs_dl_CisBP_ver2_Drosophila_melanogaster():
    motifs_name = 'CisBP_ver2_Drosophila_melanogaster.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=True)
    assert isinstance(motifs, list)

# CisBP_ver2_Drosophila_mix
def test_load_motifs_CisBP_ver2_Drosophila_mix():
    motifs_name = 'CisBP_ver2_Drosophila_mix.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=False)
    assert isinstance(motifs, list)

def test_load_motifs_dl_CisBP_ver2_Drosophila_mix():
    motifs_name = 'CisBP_ver2_Drosophila_mix.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=True)
    assert isinstance(motifs, list)

# CisBP_ver2_Drosophila_mojavensis
def test_load_motifs_CisBP_ver2_Drosophila_mojavensis():
    motifs_name = 'CisBP_ver2_Drosophila_mojavensis.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=False)
    assert isinstance(motifs, list)

def test_load_motifs_dl_CisBP_ver2_Drosophila_mojavensis():
    motifs_name = 'CisBP_ver2_Drosophila_mojavensis.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=True)
    assert isinstance(motifs, list)

# CisBP_ver2_Drosophila_persimilis
def test_load_motifs_CisBP_ver2_Drosophila_persimilis():
    motifs_name = 'CisBP_ver2_Drosophila_persimilis.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=False)
    assert isinstance(motifs, list)

def test_load_motifs_dl_CisBP_ver2_Drosophila_persimilis():
    motifs_name = 'CisBP_ver2_Drosophila_persimilis.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=True)
    assert isinstance(motifs, list)

# CisBP_ver2_Drosophila_pseudoobscura
def test_load_motifs_CisBP_ver2_Drosophila_pseudoobscura():
    motifs_name = 'CisBP_ver2_Drosophila_pseudoobscura.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=False)
    assert isinstance(motifs, list)

def test_load_motifs_dl_CisBP_ver2_Drosophila_pseudoobscura():
    motifs_name = 'CisBP_ver2_Drosophila_pseudoobscura.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=True)
    assert isinstance(motifs, list)

# CisBP_ver2_Drosophila_sechellia
def test_load_motifs_CisBP_ver2_Drosophila_sechellia():
    motifs_name = 'CisBP_ver2_Drosophila_sechellia.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=False)
    assert isinstance(motifs, list)

def test_load_motifs_dl_CisBP_ver2_Drosophila_sechellia():
    motifs_name = 'CisBP_ver2_Drosophila_sechellia.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=True)
    assert isinstance(motifs, list)

# CisBP_ver2_Drosophila_simulans
def test_load_motifs_CisBP_ver2_Drosophila_simulans():
    motifs_name = 'CisBP_ver2_Drosophila_simulans.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=False)
    assert isinstance(motifs, list)

def test_load_motifs_dl_CisBP_ver2_Drosophila_simulans():
    motifs_name = 'CisBP_ver2_Drosophila_simulans.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=True)
    assert isinstance(motifs, list)

# CisBP_ver2_Drosophila_virilis
def test_load_motifs_CisBP_ver2_Drosophila_virilis():
    motifs_name = 'CisBP_ver2_Drosophila_virilis.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=False)
    assert isinstance(motifs, list)

def test_load_motifs_dl_CisBP_ver2_Drosophila_virilis():
    motifs_name = 'CisBP_ver2_Drosophila_virilis.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=True)
    assert isinstance(motifs, list)

# CisBP_ver2_Drosophila_willistoni
def test_load_motifs_CisBP_ver2_Drosophila_willistoni():
    motifs_name = 'CisBP_ver2_Drosophila_willistoni.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=False)
    assert isinstance(motifs, list)

def test_load_motifs_dl_CisBP_ver2_Drosophila_willistoni():
    motifs_name = 'CisBP_ver2_Drosophila_willistoni.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=True)
    assert isinstance(motifs, list)

# CisBP_ver2_Drosophila_yakuba
def test_load_motifs_CisBP_ver2_Drosophila_yakuba():
    motifs_name = 'CisBP_ver2_Drosophila_yakuba.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=False)
    assert isinstance(motifs, list)

def test_load_motifs_dl_CisBP_ver2_Drosophila_yakuba():
    motifs_name = 'CisBP_ver2_Drosophila_yakuba.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=True)
    assert isinstance(motifs, list)

# CisBP_ver2_Gallus_gallus
def test_load_motifs_CisBP_ver2_Gallus_gallus():
    motifs_name = 'CisBP_ver2_Gallus_gallus.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=False)
    assert isinstance(motifs, list)

def test_load_motifs_dl_CisBP_ver2_Gallus_gallus():
    motifs_name = 'CisBP_ver2_Gallus_gallus.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=True)
    assert isinstance(motifs, list)

# CisBP_ver2_Homo_sapiens
def test_load_motifs_CisBP_ver2_Homo_sapiens():
    motifs_name = 'CisBP_ver2_Homo_sapiens.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=False)
    assert isinstance(motifs, list)

def test_load_motifs_dl_CisBP_ver2_Homo_sapiens():
    motifs_name = 'CisBP_ver2_Homo_sapiens.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=True)
    assert isinstance(motifs, list)

# CisBP_ver2_Mus_musculus
def test_load_motifs_CisBP_ver2_Mus_musculus():
    motifs_name = 'CisBP_ver2_Mus_musculus.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=False)
    assert isinstance(motifs, list)

def test_load_motifs_dl_CisBP_ver2_Mus_musculus():
    motifs_name = 'CisBP_ver2_Mus_musculus.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=True)
    assert isinstance(motifs, list)

# CisBP_ver2_Rattus_norvegicus
def test_load_motifs_CisBP_ver2_Rattus_norvegicus():
    motifs_name = 'CisBP_ver2_Rattus_norvegicus.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=False)
    assert isinstance(motifs, list)

def test_load_motifs_dl_CisBP_ver2_Rattus_norvegicus():
    motifs_name = 'CisBP_ver2_Rattus_norvegicus.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=True)
    assert isinstance(motifs, list)

# CisBP_ver2_Saccharomyces_cerevisiae
def test_load_motifs_CisBP_ver2_Saccharomyces_cerevisiae():
    motifs_name = 'CisBP_ver2_Saccharomyces_cerevisiae.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=False)
    assert isinstance(motifs, list)

def test_load_motifs_dl_CisBP_ver2_Saccharomyces_cerevisiae():
    motifs_name = 'CisBP_ver2_Saccharomyces_cerevisiae.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=True)
    assert isinstance(motifs, list)

# CisBP_ver2_Xenopus_laevis
def test_load_motifs_CisBP_ver2_Xenopus_laevis():
    motifs_name = 'CisBP_ver2_Xenopus_laevis.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=False)
    assert isinstance(motifs, list)

def test_load_motifs_dl_CisBP_ver2_Xenopus_laevis():
    motifs_name = 'CisBP_ver2_Xenopus_laevis.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=True)
    assert isinstance(motifs, list)

# CisBP_ver2_Xenopus_tropicalis
def test_load_motifs_CisBP_ver2_Xenopus_tropicalis():
    motifs_name = 'CisBP_ver2_Xenopus_tropicalis.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=False)
    assert isinstance(motifs, list)

def test_load_motifs_dl_CisBP_ver2_Xenopus_tropicalis():
    motifs_name = 'CisBP_ver2_Xenopus_tropicalis.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=True)
    assert isinstance(motifs, list)

# CisBP_ver2_Xenopus_tropicalis_and_Xenopus_laevis
def test_load_motifs_CisBP_ver2_Xenopus_tropicalis_and_Xenopus_laevis():
    motifs_name = 'CisBP_ver2_Xenopus_tropicalis_and_Xenopus_laevis.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=False)
    assert isinstance(motifs, list)

def test_load_motifs_dl_CisBP_ver2_Xenopus_tropicalis_and_Xenopus_laevis():
    motifs_name = 'CisBP_ver2_Xenopus_tropicalis_and_Xenopus_laevis.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=True)
    assert isinstance(motifs, list)
