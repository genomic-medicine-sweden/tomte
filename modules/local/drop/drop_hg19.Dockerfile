################## BASE IMAGE ######################
FROM condaforge/mambaforge:24.9.2-0

################## METADATA ######################
LABEL software="drop"
LABEL software.version="1.4.0"

# Create a non-root user and set up the home directory
ARG USERNAME=dropuser
ARG USER_UID=1000
ARG USER_GID=1000

RUN addgroup --gid ${USER_GID} ${USERNAME} && \
    adduser --disabled-password \
    --gecos "" \
    --uid ${USER_UID} \
    --gid ${USER_GID} \
    --home /home/${USERNAME} \
    ${USERNAME}

# Install the "drop" module, R packages, and clean up in a single layer
RUN mamba install -y -c conda-forge -c bioconda drop=1.4.0 pyreadr && \
    R -e 'system("defaults write org.R-project.R force.LANG C.UTF-8")' && \
    R -e 'BiocManager::install(c("BSgenome.Hsapiens.UCSC.hg19"))' && \
    R -e 'BiocManager::install(c("MafDb.gnomAD.r2.1.hs37d5"))' && \
    R -e "install.packages('DT', repos='https://cloud.r-project.org/', dependencies=TRUE)" && \
    conda clean -ya && \
    rm -rf /opt/conda/pkgs/* /tmp/* /var/tmp/*

# Set user and home directory
USER ${USERNAME}
WORKDIR /home/${USERNAME}

# Create a data directory for the module and ensure ownership
RUN mkdir /home/${USERNAME}/data && \
    chown -R ${USERNAME}:${USERNAME} /home/${USERNAME}

# Set the default user when running the container
USER ${USERNAME}

# Specify the data directory as the working directory
WORKDIR /home/${USERNAME}/data
