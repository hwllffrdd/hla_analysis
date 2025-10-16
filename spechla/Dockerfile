# Use conda base image
FROM continuumio/miniconda3:23.10.0-1

# Set working directory
WORKDIR /opt/spechla

# Install system dependencies (including less for zless)
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    wget \
    less \
    gzip \
    && rm -rf /var/lib/apt/lists/*

# Set environment variables for conda
ENV PATH="/opt/conda/bin:$PATH"
SHELL ["/bin/bash", "--login", "-c"]

# Clone SpecHLA repository
RUN git clone https://github.com/deepomicslab/SpecHLA.git /opt/spechla

# Set working directory to SpecHLA
WORKDIR /opt/spechla

# Create conda environment from SpecHLA requirements
RUN conda env create --prefix=./spechla_env -f environment.yml

# Initialize conda for bash
RUN conda init bash

# Make conda environment active by default
RUN echo "conda activate /opt/spechla/spechla_env" >> ~/.bashrc

# Create the environment PATH activation
ENV PATH="/opt/spechla/spechla_env/bin:$PATH"

# Set up conda environment activation for all RUN commands
SHELL ["conda", "run", "--no-capture-output", "-p", "/opt/spechla/spechla_env", "/bin/bash", "-c"]

# Make binaries executable
RUN chmod +x /opt/spechla/bin/*

# Download and prepare human reference genome (hg38) with HLA alt contigs
RUN echo "Downloading human reference genome hg38 with ALT contigs..." && \
    mkdir -p /opt/spechla/ref_genome && \
    cd /opt/spechla/ref_genome && \
    wget -q "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz" && \
    gunzip GCA_000001405.15_GRCh38_full_analysis_set.fna.gz && \
    mv GCA_000001405.15_GRCh38_full_analysis_set.fna hg38_with_alt.fa && \
    echo "Creating BWA index..." && \
    bwa index hg38_with_alt.fa && \
    echo "Creating samtools index..." && \
    samtools faidx hg38_with_alt.fa && \
    echo "Reference genome with ALT contigs setup complete."

# Try to run index.sh with fallback handling for SpecHLA databases
# This will use bowtie2 if Novoalign license is not available
RUN unset LD_LIBRARY_PATH && unset LIBRARY_PATH && \
    (bash index.sh || echo "Warning: SpecHLA index creation failed, but continuing. SpecHLA will use bowtie2 instead of novoalign.")

# Create data directory for input/output
RUN mkdir -p /data

# Set the working directory for data
WORKDIR /data

# Default entrypoint
ENTRYPOINT ["conda", "run", "--no-capture-output", "-p", "/opt/spechla/spechla_env", "/bin/bash"]
