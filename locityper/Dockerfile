# Use conda base image
FROM continuumio/miniconda3:23.10.0-1

# Set working directory
WORKDIR /opt/locityper

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    wget \
    curl \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    && rm -rf /var/lib/apt/lists/*

# Install Locityper and its dependencies via conda
RUN conda install -c bioconda -c conda-forge \
    locityper \
    samtools \
    && conda clean -a -y

# Create directories for references and data
RUN mkdir -p /opt/locityper/references \
    /opt/locityper/loci \
    /data

# Download a small helper script for reference setup (optional, can be run manually)
RUN echo '#!/bin/bash\n\
echo "Locityper Docker Container"\n\
echo "==========================="\n\
echo ""\n\
echo "Container is ready. Mount your data with:"\n\
echo "  -v /path/to/fastq:/input_data"\n\
echo "  -v /path/to/results:/output"\n\
echo ""\n\
echo "For reference genome setup, see /opt/locityper/setup_references.sh"\n\
echo ""\n' > /opt/locityper/welcome.sh && \
    chmod +x /opt/locityper/welcome.sh

# Create a reference setup script that users can run if needed
RUN echo '#!/bin/bash\n\
echo "Setting up Locityper references..."\n\
echo "This will download ~6GB of data and may take 20-30 minutes"\n\
echo ""\n\
cd /opt/locityper/references\n\
\n\
# Download reference genome if not present\n\
if [ ! -f "GRCh38_full_analysis_set_plus_decoy_hla.fa" ]; then\n\
  echo "Downloading GRCh38 reference..."\n\
  wget -q --show-progress https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa\n\
  wget -q https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai\n\
fi\n\
\n\
# Download pangenome VCF if not present\n\
if [ ! -f "../pangenome/hprc-v1.1-mc-grch38.raw.vcf.gz" ]; then\n\
  echo "Downloading HPRC pangenome VCF..."\n\
  mkdir -p ../pangenome\n\
  cd ../pangenome\n\
  wget -q --show-progress https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.raw.vcf.gz\n\
  wget -q https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.raw.vcf.gz.tbi\n\
  cd ../references\n\
fi\n\
\n\
# Build k-mer database if not present (uses less memory with smaller -s parameter)\n\
if [ ! -f "GRCh38_k25.jf" ]; then\n\
  echo "Building k-mer database (this may take 15-20 minutes)..."\n\
  jellyfish count -m 25 -s 1G -t 4 -C GRCh38_full_analysis_set_plus_decoy_hla.fa -o GRCh38_k25.jf\n\
fi\n\
\n\
echo ""\n\
echo "Reference setup complete!"\n\
echo "You can now run Locityper genotyping."\n' > /opt/locityper/setup_references.sh && \
    chmod +x /opt/locityper/setup_references.sh

# Set the working directory for data
WORKDIR /data

# Add conda to PATH
ENV PATH="/opt/conda/bin:$PATH"

# Run welcome message and default to bash
CMD ["/bin/bash", "-c", "/opt/locityper/welcome.sh && /bin/bash"]
