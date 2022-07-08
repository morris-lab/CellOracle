FROM ubuntu:20.04

LABEL version=""
LABEL description="celloracle docker installation"
LABEL maintainer="kamimoto@wustl.edu"

# Setup ubuntu basic softwares
RUN apt-get update \
 && apt-get install -y wget git nano gcc g++ libz-dev bedtools \
 && rm -rf /var/lib/apt/lists/*

# Clone celloracle
RUN cd \
 && git clone https://github.com/morris-lab/CellOracle.git

# Install miniconda 
RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh \
 && bash miniconda.sh -b -p $HOME/miniconda \
 && rm miniconda.sh

# Configure miniconda and install packages
RUN . "$HOME/miniconda/etc/profile.d/conda.sh" \
 && cd \
 && hash -r \
 && export PATH="$HOME/miniconda/bin:${PATH}" \
 && conda config --set always_yes yes --set changeps1 no \
 && conda update -q conda \
 && conda create -q -n celloracle_env python=3.8 \
 && conda activate celloracle_env \
 && conda install cython numpy pytest \
 && wget https://anaconda.org/bioconda/gimmemotifs/0.17.1/download/linux-64/gimmemotifs-0.17.1-py38h8ded8fe_1.tar.bz2 \
 && conda install --offline gimmemotifs-0.17.1-py38h8ded8fe_1.tar.bz2 \
 && rm gimmemotifs-0.17.1-py38h8ded8fe_1.tar.bz2 \
 && cd $HOME/CellOracle \
 && pip install . --default-timeout=100 \
 && pytest \
 && cd  \
 && rm -r CellOracle \
 && rm -r $HOME/celloracle_data \
 && conda clean --all \
 && conda init bash \
 && echo "conda activate celloracle_env" >> $HOME/.bashrc

