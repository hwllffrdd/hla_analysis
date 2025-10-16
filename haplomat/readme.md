# Hapl-o-Mat HLA Haplotype Analysis Guide

This guide outlines how to use Hapl-o-Mat for HLA haplotype frequency estimation inside a Docker container.

## 1. Docker Setup

Create a `Dockerfile`:

```dockerfile
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

# Create data directory for output
WORKDIR /data

# Add Hapl-o-Mat to PATH
ENV PATH="/opt/Hapl-o-Mat:${PATH}"

# Keep container running with bash
CMD ["/bin/bash"]
```

Build and run:

```bash
# Build the Docker image
docker build -t haplomat .

# Create local directory for data sharing
mkdir -p ~/hla_analysis/data

# Run the container with mounted volume
docker run -it --name haplomat_container -v ~/hla_analysis/data:/data haplomat
```

## 2. Data Preparation

Inside the container:

```bash
# Prepare reference data (essential first step)
cd /opt/Hapl-o-Mat/prepareData
python BuildData.py

# This creates necessary HLA reference files in /opt/Hapl-o-Mat/data
```

## 3. Input Data Conversion

Create a Python script to convert your HLA typing data to MAC format:

```python
#!/usr/bin/env python3

import os

# Input and output file paths
input_file = "/data/your_input_file.tsv"  # Change to your actual file
output_file = "/data/haplomat_input.txt"

# Create a proper MAC header line
with open(output_file, 'w') as outfile:
    # MAC format header
    outfile.write("ID\tA\tA\tB\tB\tC\tC\n")
    
    # Process the input data
    with open(input_file, 'r') as infile:
        # Skip the header line
        next(infile)
        
        for line in infile:
            fields = line.strip().split('\t')
            
            if len(fields) < 7:
                continue  # Skip malformed lines
                
            sample_id = fields[0]
            a1 = fields[1].replace("A*", "")
            a2 = fields[2].replace("A*", "")
            b1 = fields[3].replace("B*", "")
            b2 = fields[4].replace("B*", "")
            c1 = fields[5].replace("C*", "")
            c2 = fields[6].replace("C*", "")
            
            # Format MAC format with tabs
            haplomat_line = f"{sample_id}\t{a1}\t{a2}\t{b1}\t{b2}\t{c1}\t{c2}"
            outfile.write(haplomat_line + "\n")

print(f"Conversion complete! Output file: {output_file}")
```

Run the script:

```bash
python3 /data/convert_format.py
```

## 4. Create Parameter File

Create a file named `parametersMAC` in the container:

```bash
cat > /data/parametersMAC << EOF
#file names
FILENAME_INPUT=/data/haplomat_input.txt
FILENAME_HAPLOTYPES=/data/results/haplotypes.dat
FILENAME_GENOTYPES=/data/results/genotypes.dat
FILENAME_HAPLOTYPEFREQUENCIES=/data/results/hfs.dat
FILENAME_EPSILON_LOGL=/data/results/epsilon.dat
#reports
LOCI_AND_RESOLUTIONS=A:g,B:g,C:g
MINIMAL_FREQUENCY_GENOTYPES=1e-5
DO_AMBIGUITYFILTER=false
EXPAND_LINES_AMBIGUITYFILTER=false
WRITE_GENOTYPES=true
DO_ANALYTICS=false
#EM-algorithm
INITIALIZATION_HAPLOTYPEFREQUENCIES=perturbation
EPSILON=1e-6
CUT_HAPLOTYPEFREQUENCIES=1e-6
RENORMALIZE_HAPLOTYPEFREQUENCIES=true
SEED=0
EOF
```

## 5. Create Output Directory

```bash
mkdir -p /data/results
```

## 6. Run Hapl-o-Mat

From the Hapl-o-Mat directory:

```bash
cd /opt/Hapl-o-Mat
./haplomat MAC -c /data/parametersMAC
```

## 7. Output Files

Three output files will be created in `/data/results/`:

1. **hfs.dat**: Haplotype frequencies sorted in descending order
   - Format: `A*01:01g~B*08:01g~C*07:01g 0.07028014620864`
   - Each line is a haplotype with its estimated population frequency

2. **epsilon.dat**: EM algorithm convergence information
   - First column: Change in frequency estimates between iterations
   - Second column: Log-likelihood values

3. **genotypes.dat**: Processed genotype data
   - Individual ID, ambiguity codes, weight, and processed genotype

## Notes

- For different HLA loci combinations, adjust the `LOCI_AND_RESOLUTIONS` parameter
- Resolution options: g (g-group), G (G-group), P (P-group), 1f (first field), etc.
- Common issues: incorrect input format, missing reference data, non-existent output directories

To exit and remove the Docker container:
```bash
exit
docker rm haplomat_container
```
