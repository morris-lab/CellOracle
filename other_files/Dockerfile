FROM ubuntu:20.04

LABEL version="0.18.0"
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
 && bash miniconda.sh -b -p /opt/miniconda \
 && rm miniconda.sh

# Configure miniconda and install packages.
# Here, re-installing gimmemotifs using conda for safety.
RUN . "/opt/miniconda/etc/profile.d/conda.sh" \
 && cd \
 && hash -r \
 && export PATH="/opt/miniconda/bin:${PATH}" \
 && conda config --set always_yes yes --set changeps1 no \
 && conda create -q -n celloracle_env python=3.10 \
 && conda activate celloracle_env \
 && conda install cython numpy pytest -c conda-forge \
 && conda install gimmemotifs==0.17.2 -c bioconda \
 && cd $HOME/CellOracle \
 && pip install . --default-timeout=100 \
 && pytest \
 && cd  \
 && rm -r CellOracle \
 && rm -r $HOME/celloracle_data \
 && conda clean --all \
 && conda init bash \
 && echo "conda activate celloracle_env" >> $HOME/.bashrc


ENV PATH /opt/miniconda/envs/celloracle_env/bin:$PATH
ENV CONDA_DEFAULT_ENV celloracle_env