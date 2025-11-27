# --- HEADER ---
FROM ubuntu:22.04

# --- LABELS (%labels) ---
LABEL author="Beatriz Rodrigues-Estevam <bia.estevam.25@gmail.com>"
LABEL version="1.0.0 (Nov-2025)"
LABEL description="Docker image for the Nextflow pipeline SLseek (Spliced Leader sequence discovery)."
LABEL Nextflow-Pipeline="SLseek"

# Set a non-interactive environment for apt
ENV DEBIAN_FRONTEND=noninteractive

# --- ENVIRONMENT (%environment) ---
ENV MINICONDA_HOME=/opt/miniconda3
ENV CONDA_ENV_PATH=/opt/miniconda3/envs/SLseek_env
ENV PATH=$CONDA_ENV_PATH/bin:$MINICONDA_HOME/bin:$PATH
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
ENV PYTHONUNBUFFERED=1
ENV PIPELINE_CONTAINER_VERSION=1.0.0

# --- INSTALLATION (%post) ---
# Use one large RUN command to minimize Docker layers and image size.
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        wget \
        build-essential \
        git \
        bzip2 \
        ca-certificates \
        liblzma-dev \
        libbz2-dev \
        zlib1g-dev \
        libncurses5-dev \
        libncursesw5-dev \
        libssl-dev \
        libcurl4-gnutls-dev \
        libxml2-dev \
        libffi-dev \
        zlib1g \
    # Install Miniconda
    && wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh \
    && bash miniconda.sh -b -p $MINICONDA_HOME \
    && rm miniconda.sh \
    \
    # Initialize Conda and setup environment
    && . $MINICONDA_HOME/etc/profile.d/conda.sh \
    && conda config --set auto_activate_base false \
    && CONDA_CHANNELS="-c bioconda -c conda-forge -c defaults" \
    \
    # Create Conda Environment and Install ALL Packages
    && conda create -y --name SLseek_env $CONDA_CHANNELS \
        python=3.10 \
        nextflow \
        jellyfish \
        cd-hit \
        fastk \
        biopython=1.85 \
        numpy=2.2.6 \
        pandas=2.3.3 \
        scipy=1.15.3 \
        pysam=0.23.3 \
        matplotlib=3.10.7 \
        pillow=12.0.0 \
        six regex dateutil pytz packaging kiwisolver cycler fonttools pyparsing patsy statsmodels tzdata \
    \
    # Activate environment for pip install and install pip packages
    && conda activate SLseek_env \
    && pip install \
        pyahocorasick==2.2.0 \
        mizani==0.14.3 \
        plotnine==0.15.1 \
        ahocorasick \
        argparse \
    \
    # Final Cleanup
    && conda clean -a -y \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*
    
# --- ENTRYPOINT / RUNSCRIPT (%runscript) ---
COPY slseek_entrypoint.sh /usr/local/bin/slseek_entrypoint.sh
RUN chmod +x /usr/local/bin/slseek_entrypoint.sh
ENTRYPOINT ["/usr/local/bin/slseek_entrypoint.sh"]
CMD ["/bin/bash"]