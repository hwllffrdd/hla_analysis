#!/usr/bin/env python3
"""
T1K KIR Results Extraction Script
Extracts KIR typing results from T1K output files (*_t1k_kir_result_genotype.tsv)
and creates a consolidated summary table.

Input: Directory containing *_t1k_kir_result_genotype.tsv files
Output: CSV file with consolidated KIR typing results across all samples

Columns in input TSV:
gene_name num_diff_alleles allele_1 abundance_1 quality_1 allele_2 abundance_2 quality_2 secondary_alleles
"""

import os
import sys
import pandas as pd
import glob
from pathlib import Path
import argparse

def parse_t1k_kir_file(filepath):
    """
    Parse a single T1K KIR genotype file.
    
    Args:
        filepath: Path to *_t1k_kir_result_genotype.tsv file
        
    Returns:
        Dictionary with sample name and KIR typing results
    """
    # Extract sample name from filename
    filename = os.path.basename(filepath)
    # Remove the suffix '_t1k_kir_result_genotype.tsv'
    sample_name = filename.replace('_t1k_kir_result_genotype.tsv', '')
    
    # Read the TSV file
    df = pd.read_csv(filepath, sep='\t', header=None, 
                     names=['gene_name', 'num_diff_alleles', 'allele_1', 'abundance_1', 
                            'quality_1', 'allele_2', 'abundance_2', 'quality_2', 'secondary_alleles'])
    
    # Initialize results dictionary
    results = {'Sample': sample_name}
    
    # Process each KIR gene
    for _, row in df.iterrows():
        gene = row['gene_name']
        num_alleles = row['num_diff_alleles']
        
        # Store allele 1 information
        allele1 = row['allele_1']
        abundance1 = row['abundance_1']
        quality1 = row['quality_1']
        
        # Store allele 2 information (if present)
        allele2 = row['allele_2']
        abundance2 = row['abundance_2']
        quality2 = row['quality_2']
        
        # Store secondary alleles if present
        secondary = row['secondary_alleles']
        
        # Create columns for this gene
        results[f'{gene}_num_alleles'] = num_alleles
        results[f'{gene}_allele1'] = allele1
        results[f'{gene}_allele1_abundance'] = abundance1
        results[f'{gene}_allele1_quality'] = quality1
        
        if pd.notna(allele2) and allele2 != '':
            results[f'{gene}_allele2'] = allele2
            results[f'{gene}_allele2_abundance'] = abundance2
            results[f'{gene}_allele2_quality'] = quality2
        else:
            results[f'{gene}_allele2'] = ''
            results[f'{gene}_allele2_abundance'] = ''
            results[f'{gene}_allele2_quality'] = ''
        
        # Store genotype as combined alleles
        if pd.notna(allele2) and allele2 != '':
            results[f'{gene}_genotype'] = f'{allele1}/{allele2}'
        else:
            results[f'{gene}_genotype'] = allele1
            
        # Store secondary alleles if present
        if pd.notna(secondary) and secondary != '':
            results[f'{gene}_secondary'] = secondary
        else:
            results[f'{gene}_secondary'] = ''
    
    return results

def extract_all_kir_results(input_dir, output_file=None, verbose=False):
    """
    Extract KIR results from all T1K output files in a directory.
    
    Args:
        input_dir: Directory containing *_t1k_kir_result_genotype.tsv files
        output_file: Output CSV filename (optional)
        verbose: Print progress messages
        
    Returns:
        DataFrame with consolidated results
    """
    # Find all KIR result files
    pattern = os.path.join(input_dir, '*_t1k_kir_result_genotype.tsv')
    kir_files = glob.glob(pattern)
    
    if not kir_files:
        print(f"ERROR: No KIR result files found matching pattern: {pattern}")
        sys.exit(1)
    
    if verbose:
        print(f"Found {len(kir_files)} KIR result files")
    
    # Parse all files
    all_results = []
    for filepath in sorted(kir_files):
        if verbose:
            print(f"Processing: {os.path.basename(filepath)}")
        
        try:
            results = parse_t1k_kir_file(filepath)
            all_results.append(results)
        except Exception as e:
            print(f"ERROR processing {filepath}: {e}")
            continue
    
    # Create DataFrame
    df = pd.DataFrame(all_results)
    
    # Reorder columns: Sample first, then sorted by gene name
    cols = ['Sample']
    other_cols = [col for col in df.columns if col != 'Sample']
    
    # Sort other columns by gene name
    other_cols.sort()
    cols.extend(other_cols)
    
    df = df[cols]
    
    # Save to file if specified
    if output_file:
        df.to_csv(output_file, index=False)
        if verbose:
            print(f"\nResults saved to: {output_file}")
            print(f"Total samples: {len(df)}")
            print(f"Total columns: {len(df.columns)}")
    
    return df

def create_simple_summary(df, output_file=None, verbose=False):
    """
    Create a simplified summary with just genotypes for each KIR gene.
    
    Args:
        df: DataFrame from extract_all_kir_results
        output_file: Output CSV filename
        verbose: Print progress messages
        
    Returns:
        Simplified DataFrame
    """
    # Extract just Sample and genotype columns
    sample_col = ['Sample']
    genotype_cols = [col for col in df.columns if col.endswith('_genotype')]
    
    # Create simplified dataframe
    simple_df = df[sample_col + sorted(genotype_cols)].copy()
    
    # Rename columns to remove '_genotype' suffix
    rename_dict = {col: col.replace('_genotype', '') for col in genotype_cols}
    simple_df.rename(columns=rename_dict, inplace=True)
    
    if output_file:
        simple_df.to_csv(output_file, index=False)
        if verbose:
            print(f"Simple summary saved to: {output_file}")
    
    return simple_df

def create_presence_absence_matrix(df, output_file=None, verbose=False):
    """
    Create a presence/absence matrix for KIR genes across samples.
    
    Args:
        df: DataFrame from extract_all_kir_results
        output_file: Output CSV filename
        verbose: Print progress messages
        
    Returns:
        Presence/absence DataFrame
    """
    # Extract genotype columns
    sample_col = ['Sample']
    genotype_cols = [col for col in df.columns if col.endswith('_genotype')]
    
    # Create presence/absence dataframe
    pa_df = df[sample_col + sorted(genotype_cols)].copy()
    
    # Convert to presence (1) / absence (0)
    for col in genotype_cols:
        gene_name = col.replace('_genotype', '')
        pa_df[gene_name] = pa_df[col].apply(lambda x: 0 if pd.isna(x) or x == '' else 1)
        pa_df.drop(col, axis=1, inplace=True)
    
    if output_file:
        pa_df.to_csv(output_file, index=False)
        if verbose:
            print(f"Presence/absence matrix saved to: {output_file}")
    
    return pa_df

def main():
    parser = argparse.ArgumentParser(
        description='Extract and consolidate T1K KIR typing results',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Basic usage
  python t1k_kir_results_extraction.py /path/to/results_t1k
  
  # Specify output file
  python t1k_kir_results_extraction.py /path/to/results_t1k -o kir_results.csv
  
  # Create all output formats
  python t1k_kir_results_extraction.py /path/to/results_t1k --all-formats -v
        '''
    )
    
    parser.add_argument('input_dir', 
                        help='Directory containing *_t1k_kir_result_genotype.tsv files')
    parser.add_argument('-o', '--output', 
                        default='t1k_kir_consolidated_results.csv',
                        help='Output CSV file for detailed results (default: t1k_kir_consolidated_results.csv)')
    parser.add_argument('-s', '--simple', 
                        default='t1k_kir_simple_summary.csv',
                        help='Output CSV file for simple genotype summary (default: t1k_kir_simple_summary.csv)')
    parser.add_argument('-p', '--presence-absence',
                        default='t1k_kir_presence_absence.csv',
                        help='Output CSV file for presence/absence matrix (default: t1k_kir_presence_absence.csv)')
    parser             help='Create all output formats (detailed, simple, presence/absence)')
    parser.add_argument('-v', '--verbose', 
                        action='store_true',
                        help='Print verbose progress messages')
    
    args = parser.parse_args()
    
    # Check if input directory exists
    if not os.path.isdir(args.input_dir):
        print(f"ERROR: Input directory does not exist: {args.input_dir}")
        sys.exit(1)
    
    print("="*60)
    print("T1K KIR Results Extraction")
    print("="*60)
    
    # Extract detailed results
    print(f"\nProcessing files in: {args.input_dir}")
    df = extract_all_kir_results(args.input_dir, args.output, args.verbose)
    
    # Create simple summary
    if args.all_formats or args.simple:
        print("\nCreating simple genotype summary...")
        simple_df = create_simple_summary(df, args.simple, args.verbose)
        
        if args.verbose:
            print("\nFirst few rows of simple summary:")
            print(simple_df.head())
    
    # Create presence/absence matrix
    if args.all_formats or args.presence_absence:
        print("\nCreating presence/absence matrix...")
        pa_df = create_presence_absence_matrix(df, args.presence_absence, args.verbose)
        
        if args.verbose:
            print("\nFirst few rows of presence/absence matrix:")
            print(pa_df.head())
    
    print("\n" + "="*60)
    print("Processing complete!")
    print("="*60)
    
    # Print summary statistics
    print(f"\nSummary:")
    print(f"  Samples processed: {len(df)}")
    
    # Count KIR genes detected
    genotype_cols = [col for col in df.columns if col.endswith('_genotype')]
    kir_genes = set([col.replace('_genotype', '') for col in genotype_cols])
    print(f"  KIR genes detected: {len(kir_genes)}")
    print(f"  Genes: {', '.join(sorted(kir_genes))}")
    
    print(f"\nOutput files created:")
    print(f"  Detailed results: {args.output}")
    if args.all_formats or args.simple:
        print(f"  Simple summary: {args.simple}")
    if args.all_formats or args.presence_absence:
        print(f"  Presence/absence: {args.presence_absence}")

if __name__ == '__main__':
    main()
