Go to t1k directory and run

docker build -t t1k .

Start the Docker container with your volume mounts as usual:

docker run --platform linux/amd64 -it --rm -v $(pwd)/data:/input_data -v $(pwd)/results_t1k:/output t1k bash

Inside the container, run a shell loop over all R1 files and derive the sample name and corresponding R2 file:

bash
for r1 in /input_data/*_R1.fastq.gz; do
  sample=$(basename "$r1" _R1.fastq.gz)
  r2="/input_data/${sample}_R2.fastq.gz"

  /opt/t1k/run-t1k \
    -1 "$r1" \
    -2 "$r2" \
    -f /opt/t1k/hlaidx/hlaidx_dna_seq.fa \
    --preset hla-wes \
    -o "/output/${sample}_t1k_result"
done

This loop will run the T1K analysis on all pairs of samples in the /input_data directory, writing each result to the specified output directory with a folder named after each sample.

If you want to run this directly from your host machine without entering the container interactively, you can pass the loop as a command to Docker:

bash
docker run --rm -v $(pwd)/data:/input_data -v $(pwd)/results_t1k:/output t1k bash -c '\
for r1 in /input_data/*_R1.fastq.gz; do \
  sample=$(basename "$r1" _R1.fastq.gz); \
  r2="/input_data/${sample}_R2.fastq.gz"; \
  /opt/t1k/run-t1k -1 "$r1" -2 "$r2" -f /opt/t1k/hlaidx/hlaidx_dna_seq.fa --preset hla-wes -o "/output/${sample}_t1k_result"; \
done'
