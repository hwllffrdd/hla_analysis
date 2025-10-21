# HLA Analysis Pipeline Comparison

A comprehensive toolkit for comparing HLA typing tools and calculating HLA Evolutionary Divergence (HED) scores for oncogenetics research.

## Overview

This repository contains Dockerfiles, scripts, and tools for:
- Running and comparing HLA typing tools (T1K, SpecHLA, HLA-HD)
- Calculating HLA Evolutionary Divergence (HED) scores
- Analyzing peptide-MHC binding promiscuity with MHCflurry
- Haplotype analysis

## Project Structure
```
hla_analysis/
├── t1k/                          # T1K HLA typing tool
│   ├── Dockerfile
│   └── readme
├── spechla/                      # SpecHLA tool with LOH detection
│   ├── Dockerfile
│   ├── readme.md
│   └── results.md
├── hlahd/                        # HLA-HD typing tool
│   ├── Dockerfile
│   ├── hlahd.1.7.1.tar.gz
│   └── readme.md
├── mhcflurry/                    # MHC-peptide binding analysis
│   ├── peptides.txt
│   ├── mhcflurry_setup.sh
│   ├── mhcflurry_fix.sh
│   └── mhcflurry_promiscuity.py
├── hed_input_files/              # HED score calculation tools
│   └── convert_format.py
├── batch_hlahd.sh                # Batch processing for HLA-HD
├── batch_spechla.sh              # Batch processing for SpecHLA
├── enhanced_hla_comparison.py    # Compare outputs from all tools
├── extract_for_hed2.py           # Extract data for HED calculation
├── t1k_results_extraction.py     # Extract T1K results
└── wsl_venv_setup.sh             # WSL2 environment setup
```

## Prerequisites

- Docker (for containerized tools)
- Python 3.x
- WSL2 (if running on Windows)
- Git

## Setup

### 1. Clone this repository
```bash
git clone https://github.com/hwllffrdd/hla_analysis.git
cd hla_analysis
```

### 2. Set up HLA-HED (separate repository)
```bash
git clone https://github.com/CaronLab/HLA-HED.git
```
See `HLA-HED_SETUP.md` for detailed instructions.

### 3. Set up Python environment (WSL2)
```bash
bash wsl_venv_setup.sh
```

### 4. Build Docker containers
```bash
# Build T1K container
cd t1k
docker build -t t1k:latest .

# Build SpecHLA container
cd ../spechla
docker build -t spechla:latest .

# Build HLA-HD container
cd ../hlahd
docker build -t hlahd:latest .
```

## Usage

### Running HLA Typing Tools

#### Batch processing with HLA-HD
```bash
bash batch_hlahd.sh <input_directory> <output_directory>
```

#### Batch processing with SpecHLA
```bash
bash batch_spechla.sh <input_directory> <output_directory>
```

#### T1K analysis
See `t1k/readme` for specific usage instructions.

### Comparing Results
```bash
python enhanced_hla_comparison.py \
    --t1k t1k_results/ \
    --spechla spechla_results/ \
    --hlahd hlahd_results/ \
    --output comparison_report.tsv
```

### Calculating HED Scores
```bash
# Extract HLA types for HED calculation
python extract_for_hed2.py --input hla_types.txt --output hed_input.txt

# Convert to HLA-HED format
python hed_input_files/convert_format.py hed_input.txt

# Run HLA-HED (see HLA-HED_SETUP.md)
```

### MHCflurry Analysis
```bash
# Setup MHCflurry
bash mhcflurry/mhcflurry_setup.sh

# Apply fixes if needed
bash mhcflurry/mhcflurry_fix.sh

# Run promiscuity analysis
python mhcflurry/mhcflurry_promiscuity.py --peptides mhcflurry/peptides.txt
```

## Input Data Requirements

- **FASTQ files**: Paired-end sequencing data (for fair comparison across all tools)
- **Format**: `sample_1.fastq` and `sample_2.fastq` for each sample
- **Quality**: Standard Illumina quality scores (Phred+33)

## Tools Comparison Features

- **T1K**: Fast, supports both HLA and KIR typing, excellent for RNA-seq
- **SpecHLA**: Full-resolution typing with LOH detection capabilities
- **HLA-HD**: Comprehensive 29-loci typing, fastest runtime

## Citation

If you use this pipeline, please cite the original tools:

- **T1K**: [Citation from T1K paper]
- **SpecHLA**: [Citation from SpecHLA paper]
- **HLA-HD**: Kawaguchi et al. (2017) Human Mutation
- **HLA-HED**: Chowell et al. (2019) Nature Medicine

## License

Please refer to individual tool licenses in their respective directories.

## Contributing

Contributions are welcome! Please open an issue or submit a pull request.

## Contact

For questions or issues, please open a GitHub issue.
