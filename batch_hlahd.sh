#!/bin/bash
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
