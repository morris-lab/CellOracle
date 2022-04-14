 # -*- coding: utf-8 -*-
'''


'''
import pandas as pd
import numpy as np

from urllib import request
import sys
import os


from goatools.go_enrichment import GOEnrichmentStudy
from goatools.obo_parser import GODag
from goatools.associations import read_ncbi_gene2go

go_folder = os.path.join(os.path.dirname(__file__), "data")
Xtable_human=pd.read_csv(os.path.join(go_folder,'hg19_xref.txt'), sep='\t')
Xtable_mouse=pd.read_csv(os.path.join(go_folder,'biomart_xref.mm10.txt'), sep='\t')

def _check_data_and_download_if_necessary(data_folder):
    files = os.listdir(data_folder)

    if not 'gene2go.txt' in files:
        print("gene2go file was not found in the PC. Downloading gene2go file from ncbi ...")
        path = os.path.join(data_folder, "gene2go.txt.gz")
        #os.system(f"wget -O {path} https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz")
        request.urlretrieve("https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz", path)
        os.system(f"gunzip {path}")

    if not 'go-basic.obo' in files:
        print("go-basic file was not found in the PC. Downloading gene2go file from geneontology.org ...")
        path = os.path.join(data_folder, "go-basic.obo")
        os.system(f"wget -O {path} http://geneontology.org/ontology/go-basic.obo")


def geneSymbol2ID(symbols, species="mouse"):
    """
    Convert gene symbol into Entrez gene id.

    Args:
        symbols (array of str): gene symbol

        species (str): Select species. Either "mouse" or "human"

    Returns:
        list of str: Entrez gene id

    """

    if (species=='human'):
        Xtable=pd.read_csv(os.path.join(go_folder,'hg19_xref.txt'), sep='\t')

    elif(species=='mouse'):
        Xtable=pd.read_csv(os.path.join(go_folder,'biomart_xref.mm10.txt'), sep='\t')

    Xtable=Xtable[['Associated Gene Name','EntrezGene ID']].dropna()
    Xtable.index=Xtable['Associated Gene Name']
    GOIs_entrez=[x for x in np.unique(Xtable.loc[symbols].dropna()['EntrezGene ID'])]

    return GOIs_entrez


def geneID2Symbol(IDs, species="mouse"):
    """
    Convert Entrez gene id into gene symbol.

    Args:
        IDs (array of str): Entrez gene id.

        species (str): Select species. Either "mouse" or "human".

    Returns:
        list of str: Gene symbol

    """

    if (species=='human'):
        Xtable=Xtable_human

    elif(species=='mouse'):
        Xtable=Xtable_mouse

    Xtable=Xtable[['Associated Gene Name','EntrezGene ID']].dropna()
    Xtable.index=Xtable['EntrezGene ID']
    symbols=[x for x in np.unique(Xtable.loc[IDs].dropna()['Associated Gene Name'])]

    return symbols

def _ids2symbols(study_ids, species):

    if study_ids == "":
        return []
    ids = study_ids.replace(" ", "").split(",")
    ids = [int(i) for i in ids]
    genes = geneID2Symbol(ids, species)

    return genes


def get_GO(gene_query, species='mouse'):
    """
    Get Gene Ontologies (GOs).

    Args:
        gene_query (array of str): gene list.

        species (str): Select species. Either "mouse" or "human"

    Returns:
        pandas.dataframe: GO analysis results as dataframe.
    """

    sig_thresh = 3
    num_genes = None
    GOIs = gene_query

    # prepare files
    # check files
    _check_data_and_download_if_necessary(go_folder)


    obodag = GODag(os.path.join(go_folder, "go-basic.obo"))

    #go analysis

    if (species=='human'):

        geneid2gos = read_ncbi_gene2go(os.path.join(go_folder,"gene2go.txt"), taxids=[9606])
        print("{N:,} annotated genes".format(N=len(geneid2gos)))

        Xtable=pd.read_csv(os.path.join(go_folder,'hg19_xref.txt'), sep='\t')
        Xtable.index=Xtable['Approved Symbol']
        GOIs_entrez=[int(x) for x in np.unique(Xtable.loc[GOIs].dropna()['EntrezGene ID'])]

    elif (species=='mouse'):

        geneid2gos = read_ncbi_gene2go(os.path.join(go_folder,"gene2go.txt"), taxids=[10090])
        print("{N:,} annotated genes".format(N=len(geneid2gos)))

        from goatools.test_data.genes_NCBI_10090_ProteinCoding import GENEID2NT as GeneID2nt_mus

        Xtable=pd.read_csv(os.path.join(go_folder,'biomart_xref.mm10.txt'), sep='\t')
        Xtable=Xtable[['Associated Gene Name','EntrezGene ID']].dropna()
        Xtable.index=Xtable['Associated Gene Name']
        GOIs_entrez=[int(x) for x in np.unique(Xtable.loc[GOIs].dropna()['EntrezGene ID'])]


    print("processing " + str(len(GOIs)) + " genes ...")


    goeaobj = GOEnrichmentStudy(
        GeneID2nt_mus.keys(), # List of mouse protein-coding genes
        geneid2gos, # geneid/GO associations
        obodag, # Ontologies
        propagate_counts = False,
        alpha = 0.05, # default significance cut-off
        methods = ['fdr_bh']) # defult multipletest correction method


    goea_results = goeaobj.run_study(GOIs_entrez)

    li=[]
    names=[]

    go_default_output = goea_results[0].get_prtflds_default()


    for i in goea_results:
        li.append(i.get_field_values(go_default_output))
        names.append(i.name)

    df_GO = pd.DataFrame(li)

    if len(li) != 0:
        df_GO.columns = go_default_output
        df_GO["genes"] = df_GO.study_items.apply(lambda x: _ids2symbols(x, species))
    else:
        print("Found No GO with significant p-value")

    return df_GO
