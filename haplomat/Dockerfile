FROM ubuntu:22.04

# Avoid prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    libgsl-dev \
    python3 \
    python3-pip \
    vim \
    wget \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Clone Hapl-o-Mat
WORKDIR /opt
RUN git clone https://github.com/DKMS/Hapl-o-Mat.git

# Build Hapl-o-Mat
WORKDIR /opt/Hapl-o-Mat
RUN make

# Create a work directory for data
WORKDIR /data

# Add Hapl-o-Mat to PATH
ENV PATH="/opt/Hapl-o-Mat:${PATH}"

# Keep container running with bash
CMD ["/bin/bash"]
