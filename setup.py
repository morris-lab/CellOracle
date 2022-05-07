# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
from os import path
import re

with open("requirements.txt") as f:
    required = f.read().splitlines()

with open('README.rst') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

root_dir = path.abspath(path.dirname(__file__))
package_name = "celloracle"

with open(path.join(root_dir, package_name, '__init__.py')) as f:
    init_text = f.read()
    license = re.search(r'__license__\s*=\s*[\'\"](.+?)[\'\"]', init_text).group(1)
    author = re.search(r'__author__\s*=\s*[\'\"](.+?)[\'\"]', init_text).group(1)
    author_email = re.search(r'__author_email__\s*=\s*[\'\"](.+?)[\'\"]', init_text).group(1)
    url = re.search(r'__url__\s*=\s*[\'\"](.+?)[\'\"]', init_text).group(1)


with open(path.join(root_dir, package_name, 'version.py')) as f:
    init_text = f.read()
    version = re.search(r'__version__\s*=\s*[\'\"](.+?)[\'\"]', init_text).group(1)

# Start install process
setup(
    name=package_name,
    version=version, #
    description='in silico gene perturbation analysis and GRN analysis with single cell data',
    long_description=readme,
    keywords='scRNA-seq, GRN, simulation, gene perturbation',
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
                'Programming Language :: Python :: 3.7',
                'Programming Language :: Python :: 3.8'
            ],
    install_requires=required,
    author=author,
    author_email=author_email,
    url=url,
    license=license,
    package_data={"celloracle": ["go_analysis/data/*.txt", "go_analysis/data/*.obo",
                                 "data_conversion/*.R",
                                 "motif_analysis/tss_ref_data/*.bed",
                                 #"data/TFinfo_data/*.txt", "data/TFinfo_data/*.parquet",
                                 #"data/motif_data/*.txt", "data/motif_data/*.pfm",
                                 #"data/anndata/*.h5ad",
                                 #"data/tutorial_data/*.celloracle.oracle", "data/tutorial_data/*.celloracle.links",
                                 #"data/promoter_base_GRN/*.parquet",
                                 "network_analysis/rscripts_for_network_analysis/*.R",
                                 "utility/requirements.txt"]},
    packages=["celloracle", "celloracle.data_conversion", "celloracle.network", "celloracle.trajectory",
              "celloracle.data", "celloracle.go_analysis", "celloracle.oracle_utility",
              "celloracle.motif_analysis", "celloracle.network_analysis", "celloracle.utility", "celloracle.applications", "celloracle.visualizations"],
    entry_points={'console_scripts':['seuratToAnndata = celloracle.data_conversion.process_seurat_object:main']}
)
