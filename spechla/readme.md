# build the docker image
docker build -t spechla:latest .

# BAM creation
SAMPLE_NAME="007_CZEPRS2xBRCA10386_run133"

docker run --rm \
  -v ~/hla_analysis/data:/data \
  -v ~/hla_analysis/results:/results \
  spechla:latest \
  -c "
  echo 'Creating BAM file for ${SAMPLE_NAME}...' && \
  bwa mem -t 4 /opt/spechla/ref_genome/hg38_with_alt.fa \
    /data/${SAMPLE_NAME}_R1.fastq.gz \
    /data/${SAMPLE_NAME}_R2.fastq.gz | \
  samtools view -bS - | \
  samtools sort -o /data/${SAMPLE_NAME}.sorted.bam && \
  samtools index /data/${SAMPLE_NAME}.sorted.bam && \
  echo 'BAM file created successfully!'
  "

# HLA region extraction
docker run --rm \
  -v ~/hla_analysis/data:/data \
  -v ~/hla_analysis/results:/results \
  spechla:latest \
  -c "
  echo 'Extracting HLA reads...' && \
  bash /opt/spechla/script/ExtractHLAread.sh \
    -s ${SAMPLE_NAME} \
    -b /data/${SAMPLE_NAME}.sorted.bam \
    -r hg38 \
    -o /data/hla_extracted && \
  echo 'HLA extraction completed!'
  "

# SpecHLA typing
docker run --rm \
  -v ~/hla_analysis/data:/data \
  -v ~/hla_analysis/results:/results \
  spechla:latest \
  -c "
  echo 'Running SpecHLA typing...' && \
  bash /opt/spechla/script/whole/SpecHLA.sh \
    -n ${SAMPLE_NAME} \
    -1 /data/hla_extracted/${SAMPLE_NAME}_extract_1.fq.gz \
    -2 /data/hla_extracted/${SAMPLE_NAME}_extract_2.fq.gz \
    -o /results && \
  echo 'SpecHLA typing completed!'
  "
