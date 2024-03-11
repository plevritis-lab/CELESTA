# Base image 
FROM rocker/r-ubuntu:latest

##install packagaes from CRAN
RUN apt-get update && apt-get install -y \
    libudunits2-dev \
    libgdal-dev \
 && rm -rf /var/lib/apt/lists/*

RUN R -e "options(warn=2); install.packages('spdep')"

RUN install2.r --error \
    Rmixmod \
    ggplot2 \
    reshape2 \
    zeallot \
    optparse \
    remotes 
  
## Create directories
RUN mkdir -p /output \
    mkdir -p /data \
    mkdir -p /code

##install additional package CELESTA from GitHub
RUN installGithub.r plevritis/CELESTA


