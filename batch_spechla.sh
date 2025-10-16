#!/bin/bash

# SpecHLA Batch Processing Script
# This script processes all FASTQ pairs in the data directory through the complete SpecHLA pipeline

# Note: Removed 'set -e' to prevent script termination on cleanup permission errors

# Configuration
DATA_DIR="$HOME/hla_analysis/data"
RESULTS_DIR="$HOME/hla_analysis/results"
DOCKER_IMAGE="spechla:latest"
THREADS=4
REFERENCE="hg38"

# Create results directory if it doesn't exist
mkdir -p "$RESULTS_DIR"

# Function to log messages with timestamp
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# Function to process a single sample
process_sample() {
    local sample_name="$1"
    local r1_file="$2"
    local r2_file="$3"
    
    log_message "Starting processing for sample: $sample_name"
    
    # Create sample-specific result directory
    mkdir -p "$RESULTS_DIR/$sample_name"
    
    # Step 1: BAM creation
    log_message "[$sample_name] Step 1/3: Creating BAM file..."
    if ! docker run --rm \
        -v "$DATA_DIR":/data \
        -v "$RESULTS_DIR":/results \
        "$DOCKER_IMAGE" \
        -c "
        echo 'Creating BAM file for ${sample_name}...' && \
        bwa mem -t $THREADS /opt/spechla/ref_genome/hg38_with_alt.fa \
            /data/$r1_file \
            /data/$r2_file | \
        samtools view -bS - | \
        samtools sort -o /data/${sample_name}.sorted.bam && \
        samtools index /data/${sample_name}.sorted.bam && \
        echo 'BAM file created successfully for ${sample_name}!'
        "; then
        log_message "[$sample_name] ERROR: BAM creation failed"
        return 1
    fi
    
    # Step 2: HLA region extraction
    log_message "[$sample_name] Step 2/3: Extracting HLA reads..."
    if ! docker run --rm \
        -v "$DATA_DIR":/data \
        -v "$RESULTS_DIR":/results \
        "$DOCKER_IMAGE" \
        -c "
        echo 'Extracting HLA reads for ${sample_name}...' && \
        bash /opt/spechla/script/ExtractHLAread.sh \
            -s $sample_name \
            -b /data/${sample_name}.sorted.bam \
            -r $REFERENCE \
            -o /data/hla_extracted && \
        echo 'HLA extraction completed for ${sample_name}!'
        "; then
        log_message "[$sample_name] ERROR: HLA extraction failed"
        return 1
    fi
    
    # Step 3: SpecHLA typing
    log_message "[$sample_name] Step 3/3: Running SpecHLA typing..."
    if ! docker run --rm \
        -v "$DATA_DIR":/data \
        -v "$RESULTS_DIR":/results \
        "$DOCKER_IMAGE" \
        -c "
        echo 'Running SpecHLA typing for ${sample_name}...' && \
        bash /opt/spechla/script/whole/SpecHLA.sh \
            -n $sample_name \
            -1 /data/hla_extracted/${sample_name}_extract_1.fq.gz \
            -2 /data/hla_extracted/${sample_name}_extract_2.fq.gz \
            -o /results && \
        echo 'SpecHLA typing completed for ${sample_name}!'
        "; then
        log_message "[$sample_name] ERROR: SpecHLA typing failed"
        return 1
    fi
    
    # Cleanup intermediate files (optional) - handle Docker permission issues
    log_message "[$sample_name] Cleaning up intermediate files..."
    docker run --rm \
        -v "$DATA_DIR":/data \
        "$DOCKER_IMAGE" \
        -c "
        rm -f /data/${sample_name}.sorted.bam /data/${sample_name}.sorted.bam.bai 2>/dev/null || true
        rm -f /data/hla_extracted/${sample_name}_extract_1.fq.gz /data/hla_extracted/${sample_name}_extract_2.fq.gz 2>/dev/null || true
        " || log_message "[$sample_name] Warning: Some cleanup files could not be removed"
    
    log_message "[$sample_name] Processing completed successfully!"
    return 0
}

# Function to find and pair FASTQ files
find_fastq_pairs() {
    log_message "Scanning for FASTQ pairs in $DATA_DIR..."
    
    # Arrays to store sample information
    declare -A samples
    declare -a sample_names
    
    # Find all R1 files (assuming _R1 pattern)
    for r1_file in "$DATA_DIR"/*_R1.fastq* "$DATA_DIR"/*_1.fastq*; do
        if [ -f "$r1_file" ]; then
            # Extract sample name (remove path and R1/1 suffix)
            sample_name=$(basename "$r1_file" | sed -E 's/_R1\.fastq.*$|_1\.fastq.*$//')
            
            # Find corresponding R2 file
            r2_file=""
            for potential_r2 in "$DATA_DIR/${sample_name}_R2.fastq"* "$DATA_DIR/${sample_name}_2.fastq"*; do
                if [ -f "$potential_r2" ]; then
                    r2_file=$(basename "$potential_r2")
                    break
                fi
            done
            
            if [ -n "$r2_file" ]; then
                samples["$sample_name"]="$(basename "$r1_file") $r2_file"
                sample_names+=("$sample_name")
                log_message "Found pair: $sample_name -> $(basename "$r1_file"), $r2_file"
            else
                log_message "WARNING: No R2 file found for $sample_name ($(basename "$r1_file"))"
            fi
        fi
    done
    
    if [ ${#sample_names[@]} -eq 0 ]; then
        log_message "ERROR: No FASTQ pairs found in $DATA_DIR"
        log_message "Expected naming patterns: *_R1.fastq*, *_R2.fastq* or *_1.fastq*, *_2.fastq*"
        exit 1
    fi
    
    log_message "Found ${#sample_names[@]} FASTQ pairs to process"
    
    # Process each sample
    local success_count=0
    local failure_count=0
    local failed_samples=()
    
    for sample_name in "${sample_names[@]}"; do
        read -r r1_file r2_file <<< "${samples[$sample_name]}"
        
        if process_sample "$sample_name" "$r1_file" "$r2_file"; then
            ((success_count++))
        else
            ((failure_count++))
            failed_samples+=("$sample_name")
        fi
        
        log_message "Progress: $((success_count + failure_count))/${#sample_names[@]} samples processed"
    done
    
    # Final summary
    log_message "=== BATCH PROCESSING SUMMARY ==="
    log_message "Total samples: ${#sample_names[@]}"
    log_message "Successful: $success_count"
    log_message "Failed: $failure_count"
    
    if [ $failure_count -gt 0 ]; then
        log_message "Failed samples: ${failed_samples[*]}"
        exit 1
    fi
    
    log_message "All samples processed successfully!"
}

# Function to generate summary report
generate_summary_report() {
    local report_file="$RESULTS_DIR/batch_summary_report.txt"
    
    log_message "Generating summary report: $report_file"
    
    {
        echo "SpecHLA Batch Processing Summary"
        echo "Generated: $(date)"
        echo "========================================"
        echo
        
        # Process each result directory
        for result_dir in "$RESULTS_DIR"/*; do
            if [ -d "$result_dir" ] && [ "$(basename "$result_dir")" != "batch_summary_report.txt" ]; then
                sample_name=$(basename "$result_dir")
                echo "Sample: $sample_name"
                echo "------------------------"
                
                # Main HLA results
                if [ -f "$result_dir/hla.result.txt" ]; then
                    echo "HLA Typing Results:"
                    sed 's/^/  /' "$result_dir/hla.result.txt"
                else
                    echo "  ERROR: hla.result.txt not found"
                fi
                
                # G-group results
                if [ -f "$result_dir/hla.result.g.group.txt" ]; then
                    echo "G-Group Results:"
                    sed 's/^/  /' "$result_dir/hla.result.g.group.txt"
                fi
                
                # Check for potential LOH
                echo "LOH Analysis:"
                local loh_detected=false
                for locus in A B C DPA1 DPB1 DQA1 DQB1 DRB1; do
                    local allele1_file="$result_dir/hla.allele.1.HLA_${locus}.fasta"
                    local allele2_file="$result_dir/hla.allele.2.HLA_${locus}.fasta"
                    
                    if [ ! -f "$allele1_file" ] || [ ! -f "$allele2_file" ]; then
                        echo "  HLA-$locus: POTENTIAL LOH DETECTED"
                        loh_detected=true
                    else
                        local size1=$(wc -c < "$allele1_file" 2>/dev/null || echo "0")
                        local size2=$(wc -c < "$allele2_file" 2>/dev/null || echo "0")
                        if [ $size1 -eq 0 ] || [ $size2 -eq 0 ]; then
                            echo "  HLA-$locus: POTENTIAL LOH DETECTED"
                            loh_detected=true
                        fi
                    fi
                done
                
                if [ "$loh_detected" = false ]; then
                    echo "  No LOH detected"
                fi
                
                echo
                echo
            fi
        done
    } > "$report_file"
    
    log_message "Summary report generated: $report_file"
}

# Main execution
main() {
    log_message "Starting SpecHLA batch processing..."
    log_message "Data directory: $DATA_DIR"
    log_message "Results directory: $RESULTS_DIR"
    log_message "Docker image: $DOCKER_IMAGE"
    log_message "Threads: $THREADS"
    
    # Check if Docker image exists
    if ! docker image inspect "$DOCKER_IMAGE" &> /dev/null; then
        log_message "ERROR: Docker image '$DOCKER_IMAGE' not found"
        log_message "Please build the image first with: docker build -t $DOCKER_IMAGE ."
        exit 1
    fi
    
    # Check if data directory exists
    if [ ! -d "$DATA_DIR" ]; then
        log_message "ERROR: Data directory '$DATA_DIR' not found"
        exit 1
    fi
    
    # Process all samples
    find_fastq_pairs
    
    # Generate summary report
    generate_summary_report
    
    log_message "Batch processing completed successfully!"
}

# Run main function
main "$@"
