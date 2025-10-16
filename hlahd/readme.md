# HLA-HD Docker Setup and Usage Guide

## Prerequisites
- Docker installed and running
- HLA-HD source code downloaded manually from https://w3.genome.med.kyoto-u.ac.jp/HLA-HD/
- FASTQ files (paired-end) in `data/` directory with naming pattern: `*_R1.fastq.gz` and `*_R2.fastq.gz`

## Directory Structure
```
~/hla_analysis/
├── hlahd/
│   ├── Dockerfile
│   └── hlahd.1.7.1.tar.gz        # Downloaded manually
├── data/
│   ├── sample1_R1.fastq.gz
│   ├── sample1_R2.fastq.gz
│   └── ...
└── results_hlahd/                # Will be created
```

## Dockerfile
```dockerfile
FROM ubuntu:20.04
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    build-essential \
    bowtie2 \
    g++ \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt

COPY hlahd.*.tar.gz .
RUN tar -zxf hlahd.*.tar.gz && \
    rm hlahd.*.tar.gz && \
    cd hlahd.* && \
    chmod +x install.sh && \
    ./install.sh

WORKDIR /data
ENV PATH="/opt/hlahd.1.7.1/bin:$PATH"
```

## Build and Run Commands

### 1. Build Docker Image
```bash
cd ~/hla_analysis/hlahd
docker build -t hlahd .
```

### 2. Create Output Directory
```bash
cd ~/hla_analysis
mkdir -p results_hlahd
```

### 3. Run Container
```bash
docker run -it --rm \
  -v $(pwd)/data:/input_data \
  -v $(pwd)/results_hlahd:/output \
  hlahd bash
```

### 4. Single Sample Test (Inside Container)
```bash
# Set sample name (adjust as needed)
sample="007_CZEPRS2xPKM10165_run133"

# Create output directory
cd /output
mkdir -p ${sample}_hlahd

# Run HLA-HD
/opt/hlahd.1.7.1/bin/hlahd.sh \
  -t 4 \
  -m 100 \
  -c 0.95 \
  -f /opt/hlahd.1.7.1/freq_data/ \
  /input_data/${sample}_R1.fastq.gz \
  /input_data/${sample}_R2.fastq.gz \
  /opt/hlahd.1.7.1/HLA_gene.split.txt \
  /opt/hlahd.1.7.1/dictionary/ \
  ${sample} \
  ${sample}_hlahd
```

### 5. Batch Processing (Inside Container)
```bash
cd /output

for r1 in /input_data/*_R1.fastq.gz; do
  sample=$(basename "$r1" _R1.fastq.gz)
  echo "Processing sample: $sample"
  
  mkdir -p ${sample}_hlahd
  
  /opt/hlahd.1.7.1/bin/hlahd.sh \
    -t 4 \
    -m 100 \
    -c 0.95 \
    -f /opt/hlahd.1.7.1/freq_data/ \
    /input_data/${sample}_R1.fastq.gz \
    /input_data/${sample}_R2.fastq.gz \
    /opt/hlahd.1.7.1/HLA_gene.split.txt \
    /opt/hlahd.1.7.1/dictionary/ \
    ${sample} \
    ${sample}_hlahd
    
  echo "Completed: $sample"
done
```

## Understanding Results

### Output Directory Structure
```
results_hlahd/
└── {sample}_hlahd/
    └── {sample}/
        ├── result/
        │   ├── {sample}_final.result.txt    # MAIN SUMMARY FILE
        │   ├── {sample}_A.est.txt           # HLA-A typing results
        │   ├── {sample}_B.est.txt           # HLA-B typing results
        │   ├── {sample}_C.est.txt           # HLA-C typing results
        │   ├── {sample}_DRB1.est.txt        # HLA-DRB1 typing results
        │   └── ... (other HLA genes)
        ├── log/                             # Log files
        └── ... (intermediate files)
```

### Key Result Files
- **`*_final.result.txt`**: Complete summary of all HLA typing results
- **`*_{gene}.est.txt`**: Individual gene typing results (main results)
- **`*_{gene}.read.txt`**: Detailed read mapping information (for validation)

### Example Result Format
```
# From *_A.est.txt
#Pair count	1
#Best allele pair	1
HLA-A*01:01:01:01,HLA-A*01:01:01:02N	-	exon2:58.2889:comp.0,exon3:42.1449:comp.0	-
```

## HLA-HD Parameters Explained
- `-t 4`: Use 4 threads
- `-m 100`: Minimum read length
- `-c 0.95`: Trimming rate (95%)
- `-f freq_data/`: Population frequency data directory
- Input: R1 and R2 FASTQ files
- Dictionary and gene split files: HLA reference data
- Sample name and output directory

## Troubleshooting
- **Permission denied**: Make sure `install.sh` has execute permissions (handled in Dockerfile)
- **Directory not found**: Create output directories before running HLA-HD
- **No input data visible**: Check volume mount paths are correct
- **Out of memory**: Reduce thread count (`-t` parameter) or process fewer samples simultaneously

## Comparison with Other Tools
- **Input**: Uses same FASTQ files as T1K and SpecHLA
- **Output**: Provides comprehensive HLA typing for all classical HLA genes
- **Special feature**: No LOH detection (unlike SpecHLA)
