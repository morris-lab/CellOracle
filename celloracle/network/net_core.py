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

from copy import deepcopy
from datetime import datetime
from time import ctime

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from joblib import Parallel, delayed
from tqdm.notebook import tqdm

from scipy.stats import ttest_1samp, norm

from ..utility.hdf5_processing import dump_hdf5, load_hdf5

# 0.2. custom libraries
from ..utility import  save_as_pickled_object, standard, intersect
from .regression_models import get_bagging_ridge_coefs as _get_bagging_ridge_coefs
from .regression_models import get_bayesian_ridge_coefs as _get_bayesian_ridge_coefs


RIDGE_SOLVER = "auto"

##############################
### 1.  define main class  ###
##############################

class Net():
    '''
    Net is a custom class for inferring sample-specific GRN from scRNA-seq data.
    This class is used inside the Oracle class for GRN inference.
    This class requires two types of information below.

    (1) Single-cell RNA-seq data:
        The Net class needs processed scRNA-seq data.
        Gene and cell filtering, quality check, normalization, log-transformation (but not scaling and centering) have to be done before starting the GRN calculation with this class.
        You can also use any arbitrary metadata (i.e., mRNA count, cell-cycle phase) for GRN input.

    (2) Potential regulatory connection (or base GRN):
        This method uses the list of potential regulatory TFs as input.
        This information can be calculated from ATAC-seq data using the motif-analysis module.
        If sample-specific ATAC-seq data is not available,
        you can use general TF-binding info derived from public ATAC-seq dataset of various tissue/cell type.

    Attributes:
        linkList (pandas.DataFrame): The results of the GRN inference.
        all_genes (numpy.array): An array of all genes that exist in the input gene expression matrix
        embedding_name (str): The key name name in adata.obsm containing dimensional reduction coordinates
        annotation(dictionary): Annotation. you can add custom annotation.
        coefs_dict (dictionary): Coefs of linear regression.
        stats_dict (dictionary): Statistic values about coefs.
        fitted_genes (list of str): List of genes where the regression model was successfully calculated.
        failed_genes (list of str): List of genes that were not assigned coefs
        cellstate (pandas.DataFrame): A metadata for GRN input
        TFinfo (pandas.DataFrame): Information about potential regulatory TFs.
        gem (pandas.DataFrame): Merged matrix made with gene_expression_matrix and cellstate matrix.
        gem_standerdized (pandas.DataFrame): Almost the same as gem, but the gene_expression_matrix was standardized.
        library_last_update_date (str): Last update date of this code. This info is for code development. It can be deprecated in the future
        object_initiation_date (str): The date when this object was made.

    '''


    #######################################
    ### 1.1. Initialize transNet object ###
    #######################################
    def __init__(self, gene_expression_matrix, gem_standerdized=None, TFinfo_matrix=None, cellstate=None, TFinfo_dic=None, annotation=None, verbose=True):
        '''
        Instantiate Net object

        Args:
            gene_expression_matrix (pandas.DataFrame): scRNA gene expression matrix.
                Column names should be genes. Row names should be cell names.
                Data preprocessing with Seurat or scanpy is recommended.
                See tutorials for more information.

            TFinfo (pandas.DataFrame): information about potential regulatory TFs.

            cellstate (pandas.DataFrame): optional input data.
                Metadata that was acquired during data preprocessing phase can be used for GRN inference.

        '''
        ## 1. Initiate attributes
        if verbose:
            print("initiating Net object ...")
        self.gem = None
        self.gem_standerdized = None
        self.cellstate = None
        self.all_genes = None
        self.TFinfo = None
        self.TFdict = {}
        self.linkList = None
        self.annotation = annotation
        self.fitted_genes = []
        self.failed_genes = []
        self.coefs_dict = {}
        self.stats_dict = {}
        self.library_last_update_date = ctime(os.path.getmtime(__file__))
        self.object_initiation_date = datetime.now().ctime()
        if self.annotation is None:
            self.annotation = {}

        ## 2. Process attributes
        # 2.1. gene expression matrix, standerdized_gem, cellstate
        self.gem = gene_expression_matrix.copy()
        self.gem.index.name = None
        self.gem.columns.name = None

        if gem_standerdized is None:
            self.gem_standerdized = standard(self.gem)
        else:
            self.gem_standerdized = gem_standerdized.copy()
            self.gem_standerdized.index.name = None
            self.gem_standerdized.columns.name = None

        if not cellstate is None:
            self.cellstate = cellstate.copy()
            self.cellstate.index.name = None
            self.cellstate.columns.name = None
            print(f"cellstate shape: {self.cellstate.shape}")

            if not self.gem.shape[0] == self.cellstate.shape[0]:
                raise ValueError("cellnumber in GEM and cellstate shold be same.")
            else:
                self.gem = pd.concat([self.gem, self.cellstate], axis=1)
                self.gem_standerdized = pd.concat([self.gem_standerdized,
                                                   self.cellstate], axis=1)
                self.annotation["cellstate"] = list(self.cellstate.columns)
        self.gem_standerdized = self.gem_standerdized.astype("float32")
        self.gem = self.gem.astype("float32")


        if verbose:
            print(f"gem_shape: {self.gem.shape}")

        # 2.2. Gene list
        self.all_genes = np.unique(self.gem.columns)

        # 2.3. TF_information
        if not TFinfo_matrix is None:
            self.TFinfo = TFinfo_matrix.copy()
            self.TFinfo.index.name = None
            self.TFinfo.columns.name = None
            self.TFinfo = self.TFinfo.reset_index(drop=True)


            tmp = self.TFinfo.copy()
            tmp = tmp.drop(["peak_id"], axis=1)
            tmp = tmp.groupby(by="gene_short_name").sum()
            self.TFdict = dict(tmp.apply(lambda x: x[x > 0].index.values, axis=1))
            if verbose:
                print(f"TF info shape: {self.TFinfo.shape}")


        if not TFinfo_dic is None:
            self.TFdict.update(TFinfo_dic)

        if verbose:
            print("initiation completed.")

    def copy(self):
        """
        Deepcopy itself
        """
        return deepcopy(self)

    def addTFinfo_matrix(self, TFinfo_matrix):
        '''
        Load TF info dataframe.

        Args:
            TFinfo (pandas.DataFrame): information about potential regulatory TFs.

        '''
        self.TFinfo = TFinfo_matrix.copy()
        tmp = TFinfo_matrix.copy()
        tmp = tmp.drop(["peak_id"], axis=1)
        tmp = tmp.groupby(by="gene_short_name").sum()
        self.TFdict = dict(tmp.apply(lambda x: x[x > 0].index.values, axis=1))

    def updateTFinfo_dictionary(self, TFdict):
        '''
        Update TF info matrix

        Args:
            TFdict (dictionary): A python dictionary in which a key is Target gene, value are potential regulatory genes for the target gene.

        '''
        self.TFdict.update(TFdict)

    def addTFinfo_dictionary(self, TFdict):
        """
        Add a new TF info to pre-exiting TFdict.

        Args:
            TFdict (dictionary): python dictionary of TF info.
        """

        for tf in TFdict:
            if tf in self.TFdict.keys():
                targets = self.TFdict[tf]
                targets = list(TFdict[tf]) + list(targets)
                targets = np.unique(targets)
                self.TFdict.update({tf: targets})
            else:
                self.TFdict.update({tf: TFdict[tf]})



    def addAnnotation(self, annotation_dictionary):
        '''
        Add a new annotation.

        Args:
            annotation_dictionary (dictionary): e.g. {"sample_name": "NIH 3T3 cell"}
        '''

        self.annotation.update(annotation_dictionary)



    ##############################################################
    ## 2.1. Make and fit Regularized Linear model to get Coefs  ##
    ##############################################################


    def fit_All_genes_parallel(self, bagging_number=200, scaling=True, log=None, verbose=10):
        """
        IMPORTANT: this function being debugged and is currently unavailable.

        Make ML models for all genes.
        The calculation will be performed in parallel using joblib parallel module.

        Args:
            bagging_number (int): The number of estimators for bagging.
            scaling (bool): Whether or not to scale regulatory gene expression values.
            log (logging object): log object to output log
            verbose (int): verbose for joblib parallel
        """

        # prepare for parallel calculation
        def process(target_gene):

            coefs = _get_bagging_ridge_coefs(target_gene=target_gene,
                                             gem=self.gem,
                                             gem_scaled=self.gem_standerdized,
                                             TFdict=self.TFdict,
                                             cellstate=self.cellstate,
                                             bagging_number=bagging_number,
                                             scaling=scaling,
                                             n_jobs=1,
                                             alpha=alpha)
            if type(coefs) is np.int:
                fitted_gene = "na"
                failed_gene = target_gene
            else:
                fitted_gene = target_gene
                failed_gene = "na"
            return _get_stats_df_bagging_ridge(coefs), fitted_gene, failed_gene

        # do parallel calculation with joblib
        genes = np.array(intersect(self.all_genes, self.TFdict.keys()))
        results = Parallel(n_jobs=-1,
                           backend="threading",
                           verbose=verbose)([delayed(process)(i) for i in genes])

        # order and save results
        results = np.array(results)
        values, fitted_genes, failed_genes = results[:, 0], results[:, 1], results[:, 2]
        self.stats_dict = dict(zip(fitteg_genes, values))
        self.stats_dict.pop("na")
        fitted_genes = np.delete(fitted_genes, np.where(fitted_genes == "na"))
        failed_genes = np.delete(failed_genes, np.where(failed_genes == "na"))
        self.fitted_genes = fitted_genes
        self.failed_genes = failed_genes

    def fit_All_genes(self, bagging_number=200, scaling=True, model_method="bagging_ridge",
                      command_line_mode=False, log=None, alpha=1, verbose=True):
        """
        Make ML models for all genes.
        The calculation will be performed in parallel using scikit-learn bagging function.
        You can select a modeling method (bagging_ridge or bayesian_ridge).  This calculation usually takes a long time.

        Args:
            bagging_number (int): The number of estimators for bagging.
            scaling (bool): Whether or not to scale regulatory gene expression values.
            model_method (str): ML model name. Please select either "bagging_ridge" or "bayesian_ridge"
            command_line_mode (bool): Please select False if the calculation is performed on jupyter notebook.
            log (logging object): log object to output log
            alpha (int) : Strength of regularization.
            verbose (bool): Whether or not to show a progress bar.
        """
        self.fit_genes(target_genes=self.all_genes,
                       bagging_number=bagging_number,
                       scaling=scaling,
                       model_method=model_method,
                       save_coefs=False,
                       command_line_mode=command_line_mode,
                       log=log,
                       alpha=alpha,
                       verbose=verbose)

    def fit_genes(self, target_genes, bagging_number=200, scaling=True, model_method="bagging_ridge",
                  save_coefs=False, command_line_mode=False, log=None, alpha=1, verbose=True):
        """
        Make ML models for genes of interest.
        This calculation will be performed in parallel using scikit-learn's bagging function.
        You can select a modeling method; Please chose either bagging_ridge or bayesian_ridge.

        Args:
            target_genes (list of str): gene list
            bagging_number (int): The number of estimators for bagging.
            scaling (bool): Whether or not to scale regulatory gene expression values.
            model_method (str): ML model name. Please select either "bagging_ridge" or "bayesian_ridge"
            save_coefs (bool): Whether or not to store details of coef values in bagging model.
            command_line_mode (bool): Please select False if the calculation is performed on jupyter notebook.
            log (logging object): log object to output log
            alpha (int) : Strength of regularization.
            verbose (bool): Whether or not to show a progress bar.

        """
        genes = np.array(intersect(target_genes, self.TFdict.keys()))
        genes = np.array(intersect(genes, self.all_genes))
        if verbose:
            #print(f"method: {model_method}")
            if model_method == "bagging_ridge":
                #print(f"alpha: {alpha}")
                pass

        if command_line_mode:

            N = len(genes)
            log_step = 10

            if model_method == "bagging_ridge":

                for i, target_gene in enumerate(genes):
                    coefs = _get_bagging_ridge_coefs(target_gene=target_gene,
                                                     gem=self.gem,
                                                     gem_scaled=self.gem_standerdized,
                                                     TFdict=self.TFdict,
                                                     cellstate=self.cellstate,
                                                     bagging_number=bagging_number,
                                                     scaling=scaling,
                                                     alpha=alpha,
                                                     solver=RIDGE_SOLVER)

                    if isinstance(coefs, np.int):
                        self.failed_genes.append(target_gene)

                    else:
                        self.fitted_genes.append(target_gene)
                        self.stats_dict[target_gene] = _get_stats_df_bagging_ridge(coefs)
                        if save_coefs:
                            self.coefs_dict[target_gene] = coefs

                    if i/N*100 >= log_step:
                        #print(f"{datetime.now().ctime()}: {log_step}% processed..")
                        log.info(f"{log_step}% processed..")
                        log_step += 10

            elif model_method == "bayesian_ridge":

                for i, target_gene in enumerate(genes):
                    coef_mean, coef_variance, coef_names = \
                        _get_bayesian_ridge_coefs(target_gene=target_gene,
                                                  gem=self.gem,
                                                  gem_scaled=self.gem_standerdized,
                                                  TFdict=self.TFdict,
                                                  cellstate=self.cellstate,
                                                  scaling=True)

                    if isinstance(coef_mean, np.int):
                        self.failed_genes.append(target_gene)

                    else:
                        self.fitted_genes.append(target_gene)
                        stats_df = _get_stats_df_from_bayesian_ridge(coef_mean=coef_mean,
                                                                     coef_variance=coef_variance,
                                                                     coef_names=coef_names)
                        self.stats_dict[target_gene] = stats_df

                    if i/N*100 >= log_step:
                        #print(f"{datetime.now().ctime()}: {log_step}% processed..")
                        log.info(f"{log_step}% processed..")
                        log_step += 10

            log.info("100% processed..")
            #print(f"{datetime.now().ctime()}: 100% processed..")

        else:

            if model_method == "bagging_ridge":
                if verbose:
                    loop = tqdm(genes)
                else:
                    loop = genes

                for target_gene in loop:
                    coefs = _get_bagging_ridge_coefs(target_gene=target_gene,
                                                     gem=self.gem,
                                                     gem_scaled=self.gem_standerdized,
                                                     TFdict=self.TFdict,
                                                     cellstate=self.cellstate,
                                                     bagging_number=bagging_number,
                                                     scaling=scaling,
                                                     alpha=alpha,
                                                     solver=RIDGE_SOLVER)

                    if isinstance(coefs, np.int):
                        self.failed_genes.append(target_gene)

                    else:
                        self.fitted_genes.append(target_gene)
                        self.stats_dict[target_gene] = _get_stats_df_bagging_ridge(coefs)

            elif model_method == "bayesian_ridge":

                for target_gene in tqdm(genes):
                    coef_mean, coef_variance, coef_names = \
                        _get_bayesian_ridge_coefs(target_gene=target_gene,
                                                  gem=self.gem,
                                                  gem_scaled=self.gem_standerdized,
                                                  TFdict=self.TFdict,
                                                  cellstate=self.cellstate,
                                                  scaling=True)

                    if isinstance(coef_mean, np.int):
                        self.failed_genes.append(target_gene)

                    else:
                        self.fitted_genes.append(target_gene)
                        stats_df = _get_stats_df_from_bayesian_ridge(coef_mean=coef_mean,
                                                                     coef_variance=coef_variance,
                                                                     coef_names=coef_names)
                        self.stats_dict[target_gene] = stats_df

    #########################
    ### 3.1. Plot results ###
    #########################

    def plotCoefs(self, target_gene, sort=True, threshold_p=None):
        """
        Plot the distribution of Coef values (network edge weights).

        Args:
            target_gene (str): gene name
            sort (bool): Whether or not to sort genes by its strength
            bagging_number (int): The number of estimators for bagging.
            threshold_p (float): the threshold for p-values. TFs will be filtered based on the p-value.
                if None, no filtering is applied.
        """
        coefs = self.coefs_dict[target_gene].copy()
        stat = self.stats_dict[target_gene].copy()

        if not threshold_p is None:
            sig_TFs = (stat.percentile_of_zero > (100-threshold_p*100)) |\
                (stat.percentile_of_zero < (threshold_p*100))
            sig_TFs = stat.index.values[sig_TFs]
            coefs = coefs[sig_TFs]
            stat = stat.loc[sig_TFs]

        if sort:
            sorted_columns = stat.sort_values(by="positive_score",
                                              ascending=False).index.values
            coefs = coefs[sorted_columns]

        melted = _get_melted_df(coefs)
        sns.boxplot(data=melted, x="variable", y="value")
        plt.title(target_gene)
        plt.axhline(0)
        plt.xticks(range(len(coefs.columns)), coefs.columns, rotation="vertical")
    #############################
    ### 4.1. Process results  ###
    #############################

    # get linkList
    def updateLinkList(self, verbose=True):
        """
        Update LinkList.
        LinkList is a data frame that store information about inferred GRNs.

        Args:
            verbose (bool): Whether or not to show a progress bar

        """
        if not self.fitted_genes: # if the sequence is empty
            print("No model found. Do fit first.")

        linkList = []
        if verbose:
            loop = tqdm(np.unique(self.fitted_genes))
        else:
            loop = np.unique(self.fitted_genes)

        for i in loop:
            linkList.append(_stats2LinkList(i, stat_df=self.stats_dict[i]))

        linkList = pd.concat(linkList, axis=0)
        linkList = linkList.reset_index(drop=True)

        self.linkList = linkList





    ####################################
    ### 5.1 Save and Output results  ###
    ####################################


    def _save_as_pickle(self, folder):

        """
        Save itself as a pickled data.
        The pickled data may be large.
        Please use "load_pickled_object" function in utility module to load saved file.

        Args:
            foder (str): path to save the Net object.
                The folder will be created if the folder does not exist.

        """

        os.makedirs(folder, exist_ok=True)
        save_as_pickled_object(self, os.path.join(folder, "transNet.pickle"))

        print(f"file saved at: {folder}")


    def _save_as_parquet(self, folder=None):
        """
        Save itself. DataFrame in this object will be saved as a parquet file.
        Please use "load_compressed_TN_object" function in the network module to load the saved file.

        Args:
            folder (str): path to save the Net object.
                The folder will be created if the folder does not exist.

        """

        os.makedirs(folder, exist_ok=True)

        tmp_self = self.copy()

        tmp_self.gem.to_parquet(folder + "/gem.parquet")
        tmp_self.gem = None

        tmp_self.gem_standerdized.to_parquet(folder + "/gem_standerdized.parquet")
        tmp_self.gem_standerdized = None

        if not self.linkList is None:
            tmp_self.linkList.to_parquet(folder + "/linkList.parquet")
            tmp_self.linkList = None

        if not self.TFinfo is None:
            tmp_self.TFinfo.to_parquet(folder + "/TFinfo.parquet")
            tmp_self.TFinfo = None

        if not self.cellstate is None:
            tmp_self.cellstate.to_parquet(folder + "/cellstate.parquet")
            tmp_self.cellstate = None


        save_as_pickled_object(tmp_self, folder + "/transNet.pickle")

        print(f"file saved at: {folder}")


    def to_hdf5(self, file_path):
        """
        Save object as hdf5.

        Args:
            file_path (str): file path to save file. Filename needs to end with '.celloracle.net'
        """
        if file_path.endswith(".celloracle.net"):
            pass
        else:
            raise ValueError("Filename needs to end with '.celloracle.net'")

        compression_opts = 7
        dump_hdf5(obj=self, filename=file_path,
                  data_compression=compression_opts,  chunks=(2048, 2048),
                  noarray_compression=compression_opts, pickle_protocol=4)

####################################################
### 2. Define functions for transNet calculation ###
####################################################



# this is a function to convert DataFrame to visualize itself by seaborn plot function.
def _get_melted_df(df):
    d_ = df.copy()
    d_["index_name"] = df.index
    melted = pd.melt(d_, id_vars=['index_name'], var_name='variable', value_name='value')
    return melted

# this function process coefs to get several statistical values
def _get_stats_df_bagging_ridge(df):

    if isinstance(df, np.int):
        return 0

    mean = df.mean()
    p = df.apply(lambda x: ttest_1samp(x.dropna(), 0)[1])
    neg_log_p = -np.log10(p.fillna(1))

    result = pd.concat([mean, mean.abs(),
                        p, neg_log_p, #positive_score, negative_score
                        ], axis=1, sort=False)
    result.columns = ["coef_mean", "coef_abs", "p", "-logp",
                      #"positive_score", "negative_score"
                      ]

    return result

def _get_stats_df_from_bayesian_ridge(coef_mean, coef_variance, coef_names):

    coef_abs = np.abs(coef_mean)
    p = norm.cdf(x=0, loc=coef_abs, scale=np.sqrt(coef_variance))*2
    neg_log_p = -np.log(p)

    stats_df = pd.DataFrame({"coef_mean": coef_mean,
                             "coef_abs": coef_abs,
                             "coef_variance": coef_variance,
                             "p": p,
                             "-logp": neg_log_p},
                             index = coef_names)
    return stats_df

def _stats2LinkList(tg, stat_df):

    links = pd.DataFrame({"source": stat_df.index.values,
                          "target": np.repeat(tg, len(stat_df))})
    linkList = pd.concat([links, stat_df.reset_index(drop=True)], axis=1)

    return linkList
