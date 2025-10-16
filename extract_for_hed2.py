#!/usr/bin/env python3
"""
Script to extract HLA typing data from comparison CSV and format for HLA-HED calculation.
Focuses on T1K results with handling for:
- Multiple ambiguous allele calls (takes first/best)
- LOH cases (duplicates single allele to represent homozygosity)
"""

import pandas as pd
import os
import sys
from pathlib import Path

def clean_and_select_allele(allele_string, sample_id=None, locus=None):
    """
    Clean allele string and handle multiple comma-separated alleles.
    Takes the first allele if multiple are present (most confident).
    Examples: 
      'B*15:01:01,HLA-B*15:01:60,HLA-B*15:456' -> 'B*15:01:01'
      'A*02:01' -> 'A*02:01'
    """
    if pd.isna(allele_string) or allele_string == '':
        return None
    
    allele_clean = str(allele_string).strip()
    
    # Handle comma-separated multiple alleles - take the first one
    if ',' in allele_clean:
        first_allele = allele_clean.split(',')[0].strip()
        if sample_id and locus:
            print(f"  {sample_id} - {locus}: Multiple alleles detected, using first: {first_allele}")
        return first_allele
    
    return allele_clean

def extract_locus_from_allele(allele_string):
    """
    Extract HLA locus (A, B, or C) from allele string.
    Examples: 'A*02:01' -> 'A', 'HLA-B*07:02' -> 'B'
    """
    if pd.isna(allele_string) or allele_string == '':
        return None
    
    allele_clean = str(allele_string).strip().upper()
    
    # Handle different formats
    if allele_clean.startswith('HLA-'):
        allele_clean = allele_clean[4:]  # Remove 'HLA-' prefix
    
    if '*' in allele_clean:
        locus = allele_clean.split('*')[0]
        return locus
    
    return None

def format_allele(allele_string, sample_id=None, locus=None):
    """
    Format allele to standard 2-field format (Locus*XX:XX).
    Handles multiple comma-separated alleles by taking the first.
    """
    if pd.isna(allele_string) or allele_string == '':
        return None
    
    # First, handle multiple alleles
    allele_clean = clean_and_select_allele(allele_string, sample_id, locus)
    if allele_clean is None:
        return None
    
    allele_clean = allele_clean.strip().upper()
    
    # Remove HLA- prefix if present
    if allele_clean.startswith('HLA-'):
        allele_clean = allele_clean[4:]
    
    # Ensure 2-field format (truncate if longer)
    if '*' in allele_clean and ':' in allele_clean:
        parts = allele_clean.split('*')
        if len(parts) == 2:
            locus = parts[0]
            field_parts = parts[1].split(':')
            if len(field_parts) >= 2:
                # Keep only first two fields
                formatted = f"{locus}*{field_parts[0]}:{field_parts[1]}"
                return formatted
    
    return allele_clean

def organize_by_locus(sample_data, tool_prefix, sample_id=None):
    """
    Organize alleles by locus (A, B, C) for a specific tool.
    Handles LOH by duplicating single allele when second is missing.
    Returns dict with loci as keys and [allele1, allele2] as values.
    """
    allele1_col = f"{tool_prefix}_Allele1"
    allele2_col = f"{tool_prefix}_Allele2"
    
    if allele1_col not in sample_data or allele2_col not in sample_data:
        return None
    
    allele1_raw = sample_data[allele1_col]
    allele2_raw = sample_data[allele2_col]
    
    # Extract locus first to use in logging
    temp_locus = extract_locus_from_allele(allele1_raw) or extract_locus_from_allele(allele2_raw)
    
    # Format both alleles (with logging for multiple alleles)
    allele1 = format_allele(allele1_raw, sample_id, temp_locus)
    allele2 = format_allele(allele2_raw, sample_id, temp_locus)
    
    # Handle LOH cases - only log for A, B, C loci
    if allele1 is None and allele2 is None:
        # Both missing - skip this locus
        return None
    elif allele1 is not None and allele2 is None:
        # LOH detected - only allele1 present, duplicate it
        locus = extract_locus_from_allele(allele1)
        if locus in ['A', 'B', 'C'] and sample_id:
            print(f"  {sample_id} - {locus}: LOH detected, duplicating {allele1}")
        allele2 = allele1
    elif allele1 is None and allele2 is not None:
        # LOH detected - only allele2 present, use it for both
        locus = extract_locus_from_allele(allele2)
        if locus in ['A', 'B', 'C'] and sample_id:
            print(f"  {sample_id} - {locus}: LOH detected, duplicating {allele2}")
        allele1 = allele2
    
    # Extract locus from first allele
    locus1 = extract_locus_from_allele(allele1)
    locus2 = extract_locus_from_allele(allele2)
    
    # Verify loci match (they should if LOH handling worked correctly)
    if locus1 != locus2 or locus1 is None:
        if sample_id:
            print(f"  {sample_id}: Warning - Inconsistent loci: {allele1}, {allele2}")
        return None
    
    return {
        'locus': locus1,
        'allele1': allele1,
        'allele2': allele2
    }

def extract_t1k_data(input_csv_path, output_dir):
    """
    Extract and format T1K HLA data for HLA-HED.
    """
    # Read the comparison CSV
    try:
        df = pd.read_csv(input_csv_path)
        print(f"Loaded {len(df)} rows from {input_csv_path}")
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        return
    
    # Create output directory if it doesn't exist
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    tool = 'T1K'
    print(f"\nProcessing {tool} results...")
    print("=" * 60)
    print("Reporting irregularities only (LOH in A/B/C, multiple alleles)...")
    print("=" * 60)
    
    # Group by sample and organize by locus
    samples_data = {}
    
    for idx, row in df.iterrows():
        sample_id = row['Sample_ID']
        
        if sample_id not in samples_data:
            samples_data[sample_id] = {'A': {}, 'B': {}, 'C': {}}
        
        # Process this row for T1K (logging happens inside organize_by_locus)
        locus_data = organize_by_locus(row, tool, sample_id)
        
        if locus_data:
            locus = locus_data['locus']
            if locus in ['A', 'B', 'C']:
                samples_data[sample_id][locus] = {
                    'allele1': locus_data['allele1'],
                    'allele2': locus_data['allele2']
                }
    
    # Create HLA-HED input format
    hed_data = []
    skipped_samples = []
    
    print("\n" + "=" * 60)
    print("Creating HLA-HED input file...")
    print("=" * 60)
    
    for sample_id, loci in samples_data.items():
        # Check if we have all three loci
        missing_loci = [locus for locus in ['A', 'B', 'C'] if locus not in loci or not loci[locus]]
        
        if missing_loci:
            print(f"Skipping {sample_id} - missing loci: {', '.join(missing_loci)}")
            skipped_samples.append(sample_id)
            continue
        
        hed_row = {
            'Sample': sample_id,
            'A1': loci['A']['allele1'],
            'A2': loci['A']['allele2'],
            'B1': loci['B']['allele1'],
            'B2': loci['B']['allele2'],
            'C1': loci['C']['allele1'],
            'C2': loci['C']['allele2']
        }
        hed_data.append(hed_row)
    
    # Write to tab-delimited file
    if hed_data:
        output_file = os.path.join(output_dir, "t1k_hed_input.txt")
        
        hed_df = pd.DataFrame(hed_data)
        hed_df.to_csv(output_file, sep='\t', index=False)
        
        print(f"\nSuccess! Created {output_file}")
        print(f"Total samples processed: {len(hed_data)}")
        print(f"Samples skipped: {len(skipped_samples)}")
        
        # Show first few rows for verification
        print(f"\nFirst 5 rows of output:")
        print("=" * 60)
        print(hed_df.head(5).to_string(index=False))
        
        # Show summary statistics
        print(f"\n" + "=" * 60)
        print("Summary:")
        print(f"  - Successfully formatted samples: {len(hed_data)}")
        print(f"  - Skipped samples (incomplete data): {len(skipped_samples)}")
        
        if skipped_samples:
            print(f"\nSkipped samples list:")
            for sample in skipped_samples[:10]:  # Show first 10
                print(f"  - {sample}")
            if len(skipped_samples) > 10:
                print(f"  ... and {len(skipped_samples) - 10} more")
    else:
        print(f"Error: No valid data found for {tool}")

def main():
    # Configuration - using Unix/WSL paths
    input_csv = "./enhanced_comparison_results/classical_hla_comparison.csv"
    output_directory = "./hed_input_files"
    
    print("=" * 60)
    print("T1K HLA-HED Data Extraction Script")
    print("=" * 60)
    print(f"Input file: {input_csv}")
    print(f"Output directory: {output_directory}")
    print("\nFeatures:")
    print("  - Handles multiple ambiguous allele calls")
    print("  - Handles LOH (duplicates single allele)")
    print("  - 2-field resolution format")
    print("=" * 60)
    
    # Check if input file exists
    if not os.path.exists(input_csv):
        print(f"Error: Input file not found: {input_csv}")
        return
    
    # Extract the data
    extract_t1k_data(input_csv, output_directory)
    
    print("\n" + "=" * 60)
    print("Extraction complete!")
    print("=" * 60)
    print(f"\nOutput file: {output_directory}/t1k_hed_input.txt")
    print("\nNext step - Run HLA-HED:")
    print("  cd /mnt/c/hla_analysis/HLA-HED")
    print("  python3 hla_hed.py -d database/grantham_matrix.txt -f database/ABC_prot.fa \\")
    print("         -i ../hed_input_files/t1k_hed_input.txt \\")
    print("         -o ../hed_input_files/t1k_hed_output.txt")

if __name__ == "__main__":
    main()
