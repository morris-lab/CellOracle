# -*- coding: utf-8 -*-

from setuptools import setup, find_packages


with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

"""with open('requirements.txt') as f:
    requirements = f.read()
"""
setup(
    name='celloracle',
    version='0.6.3',
    description='GRN analysis with single cell data',
    long_description=readme,
    python_requires='>=3.6',
    classifiers=[# How mature is this project? Common values are
                #   3 - Alpha
                #   4 - Beta
                #   5 - Production/Stable
                'Development Status :: 4 - Beta',

                # Indicate who your project is intended for
                'Intended Audience :: Developers',
                'Topic :: Software Development :: Build Tools',

                # Pick your license as you wish (should match "license" above)
                # 'License :: OSI Approved :: MIT License',

                # Specify the Python versions you support here. In particular, ensure
                # that you indicate whether you support Python 2, Python 3 or both.
                'Programming Language :: Python :: 3.6',
                'Programming Language :: Python :: 3.7'
            ],
    install_requires=['numpy',
                      'scipy',
                      'cython',
                      'numba',
                      'matplotlib',#'matplotlib==3.0.*',
                      'seaborn',
                      'scikit-learn',
                      'h5py>=3.1.0',
                      'pandas>=1.0.3', #'pandas==0.25.*',
                      'velocyto>=0.17',
                      'umap-learn', #'umap-learn==0.3.10',
                      'pyarrow>=0.17',
                      'jupyter',
                      'tqdm>=4.45',
                      "python-igraph",
                      "louvain", #"louvain==0.6.1",
                      "fa2",
                      "anndata>=0.7.5",
                      'scanpy>=1.6', #'scanpy==1.4.4',
                      'joblib',
                      'ipyvuetify',
                      'goatools',
                      'genomepy==0.8.4', #0.5.5
                      'gimmemotifs==0.14.4'], # 0.13.1
    author='Kenji Kamimoto at Samantha Morris Lab',
    author_email='kamimoto@wustl.edu',
    url='https://github.com/morris-lab/CellOracle',
    license=license,
    package_data={"celloracle": ["go_analysis/data/*.txt", "go_analysis/data/*.obo",
                                 "data_conversion/*.R",
                                 "motif_analysis/tss_ref_data/*.bed",
                                 "data/TFinfo_data/*.txt", "data/TFinfo_data/*.parquet", "data/motif_data/*.txt", "data/motif_data/*.pfm",
                                 "network_analysis/rscripts_for_network_analysis/*.R"]},
    packages=["celloracle", "celloracle.data_conversion", "celloracle.network", "celloracle.trajectory",
              "celloracle.data", "celloracle.go_analysis", "celloracle.oracle_utility",
              "celloracle.motif_analysis", "celloracle.network_analysis", "celloracle.utility", "celloracle.applications", "celloracle.visualizations"],
    entry_points={'console_scripts':['seuratToAnndata = celloracle.data_conversion.process_seurat_object:main']}

)
