# Get base image
FROM bioconductor/bioconductor_docker:RELEASE_3_16

# environment
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH $PATH:/opt/conda/bin

# get conda installed
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean --all -f -y && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    apt-get update

COPY environment.yaml .
RUN /opt/conda/bin/conda config --add channels defaults && \
    /opt/conda/bin/conda config --add channels bioconda && \
    /opt/conda/bin/conda config --add channels conda-forge && \
    /opt/conda/bin/conda install -c conda-forge mamba && \
    mamba env update -n base -f environment.yaml && \
    conda clean -afy


# Install R packages
WORKDIR /home/rstudio
COPY r-pkgs.R r-pkgs.R
RUN R -e 'source("r-pkgs.R")'