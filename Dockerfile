FROM rocker/r-ver:4.5.3 AS builder
LABEL author="zdwyer@gmail.com" description="Builder image for TCGA-database"

ENV DEBIAN_FRONTEND=noninteractive
ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8

RUN apt-get update && apt-get install -y libcurl4-openssl-dev libssl-dev libxml2-dev zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Install required R packages
RUN R -e "install.packages(c('dplyr', 'tidyr', 'duckdb', 'BiocManager'), repos='https://packagemanager.posit.co/cran/__linux__/noble/2026-04-01')" && \
    R -e "BiocManager::install(c('TCGAbiolinks', 'SummarizedExperiment'), ask=FALSE, update=FALSE)"


RUN R -q -e "library(dplyr); library(tidyr); library(duckdb); library(TCGAbiolinks); library(SummarizedExperiment);"

FROM rocker/r-ver:4.5.3
LABEL author="zdwyer@gmail.com" description="Runtime image for TCGA-database"

ENV DEBIAN_FRONTEND=noninteractive
ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8

RUN apt-get update && apt-get install -y \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /usr/local/lib/R/site-library /usr/local/lib/R/site-library