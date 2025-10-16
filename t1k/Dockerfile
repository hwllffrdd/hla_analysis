# Use conda base image
FROM continuumio/miniconda3:23.10.0-1

# Set working directory
WORKDIR /opt/t1k

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    wget \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Install samtools via conda (T1K dependency)
RUN conda install -c bioconda samtools

# Clone T1K repository
RUN git clone https://github.com/mourisl/T1K.git /opt/t1k

# Set working directory to T1K
WORKDIR /opt/t1k

# Compile T1K
RUN make

# Download and build HLA reference database
RUN perl t1k-build.pl -o hlaidx --download IPD-IMGT/HLA

# Download and build KIR reference database  
RUN perl t1k-build.pl -o kiridx --download IPD-KIR

# Create data directory for input/output
RUN mkdir -p /data

# Set the working directory for data
WORKDIR /data

# Set PATH to include T1K
ENV PATH="/opt/t1k:$PATH"

# Default entrypoint
ENTRYPOINT ["/bin/bash"]
