# Get base image
FROM condaforge/mambaforge:24.1.2-0

# environment
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH $PATH:/opt/conda/bin
ENV NUMBA_CACHE_DIR=/tmp/numba_cache


RUN  apt-get update \
  && apt-get install -y wget \
  && apt-get -y install sudo \
  && sudo apt update \
  && sudo apt install -y zip \
  && sudo apt install -y unzip

COPY environment.yaml .
RUN /opt/conda/bin/conda config --add channels defaults && \
    /opt/conda/bin/conda config --add channels bioconda && \
    /opt/conda/bin/conda config --add channels conda-forge && \
    /opt/conda/bin/mamba install -n base python=3.12 && \
    /opt/conda/bin/mamba env update -n base -f environment.yaml && \
    /opt/conda/bin/mamba install --force-reinstall java-jdk && \
    /opt/conda/bin/mamba clean -afy