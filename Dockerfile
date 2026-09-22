# Base image: RStudio + R (Linux)
FROM rocker/rstudio:4.3.1

RUN apt-get update && apt-get install -y \
    wget \
    curl \
    git \
    build-essential \
    cmake \
    python3 \
    python3-pip \
    libncurses-dev \
    libdeflate-dev \
    libssl-dev \
    libbz2-dev \
    zlib1g-dev \
    liblzma-dev \
    clustalo \
    libxml2-dev \
    libblas-dev \
    liblapack-dev \
    muscle \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /opt/tools

##############################
# Install samtools from source
##############################
RUN curl -L -O https://github.com/samtools/samtools/releases/download/1.22/samtools-1.22.tar.bz2 \
    && tar -jxf samtools-1.22.tar.bz2 \
    && cd samtools-1.22 \
    && ./configure --prefix=/usr/local \
    && make -j$(nproc) \
    && make install \
    && cd /opt/tools \
    && rm -rf samtools-1.22 samtools-1.22.tar.bz2

##############################
# Install bcftools from source
##############################
RUN curl -L -O https://github.com/samtools/bcftools/releases/download/1.22/bcftools-1.22.tar.bz2 \
    && tar -jxf bcftools-1.22.tar.bz2 \
    && cd bcftools-1.22 \
    && ./configure --prefix=/usr/local \
    && make -j$(nproc) \
    && make install \
    && cd /opt/tools \
    && rm -rf bcftools-1.22 bcftools-1.22.tar.bz2

##############################
# Install R packages
##############################
RUN R -e "install.packages(c('shiny','data.table','plotly','DT'), repos='https://cran.rstudio.com/')"
RUN R -e "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager'); BiocManager::install(c('Rsamtools','Biostrings'))"
# CRAN packages
RUN R -e "install.packages(c('ggplot2','data.table','shiny','plotly','DT'), repos='https://cran.rstudio.com/')"

# Bioconductor packages
RUN R -e "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager'); BiocManager::install(c('Biostrings','msa'))"

# Bio3D package
RUN R -e "install.packages('bio3d', repos='https://cran.rstudio.com/')"

# r3dmol from CRAN
RUN R -e "install.packages('r3dmol', repos='https://cran.rstudio.com/')"

# Expose RStudio and Shiny ports
EXPOSE 8787 3838


# Default command
CMD ["/init"]
