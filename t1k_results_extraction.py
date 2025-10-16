#!/usr/bin/env python3
"""
T1K HLA Results Extraction Script
Extracts comprehensive HLA typing data from T1K results including:
- Classical HLA genes (HLA-A, B, C, DRB1, DQB1, DPB1, DQA1, DPA1)
- Non-classical HLA genes (HLA-E, F, G, H, J, K, L)
- Quality scores and read counts
- Alternative alleles
- Novel variant detection from VCF files

Supports both single sample and batch analysis.
"""

import pandas as pd
import os
from pathlib import Path
from typing import Dict, List

class T1KExtractor:
    def __init__(self, base_dir: str):
        self.base_dir = Path(base_dir)
        self.t1k_dir = self.base_dir
        
        # Classical HLA genes
        self.classical_genes = [
            'HLA-A', 'HLA-B', 'HLA-C',           # Class I
            'HLA-DRB1', 'HLA-DQB1', 'HLA-DPB1',  # Class II major
            'HLA-DQA1', 'HLA-DPA1'               # Class II alpha
        ]
        
        # Non-classical HLA genes of oncological interest
        self.nonclassical_genes = [
            'HLA-E', 'HLA-F', 'HLA-G',  # Immune evasion genes
            'HLA-H', 'HLA-J', 'HLA-K', 'HLA-L'  # Potential biomarkers
        ]
        
        self.all_genes = self.classical_genes + self.nonclassical_genes
    
    def extract_t1k_results(self, sample_id: str) -> Dict:
        """Extract comprehensive T1K results with extended gene coverage"""
        results = {
            'sample_id': sample_id,
            'tool': 'T1K',
            'resolution': '4-field',
            'classical_hla': {},
            'nonclassical_hla': {},
            'other_genes': {},
            'quality_scores': {},
            'read_counts': {},
            'alternatives': {},
            'novel_variants': {},
            'errors': []
        }
        
        try:
            # Main genotype results
            genotype_file = self.t1k_dir / f"{sample_id}_t1k_result_genotype.tsv"
            if not genotype_file.exists():
                results['errors'].append(f"Genotype file not found: {genotype_file}")
                return results
            
            with open(genotype_file, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 8:
                        gene = parts[0]
                        num_alleles = int(parts[1])
                        
                        gene_result = {
                            'num_alleles': num_alleles,
                            'allele_1': parts[2] if parts[2] != '.' else None,
                            'score_1': float(parts[3]) if parts[3] != '.' else 0,
                            'reads_1': int(parts[4]) if parts[4] != '-1' else 0,
                            'allele_2': parts[5] if parts[5] != '.' else None,
                            'score_2': float(parts[6]) if parts[6] != '.' else 0,
                            'reads_2': int(parts[7]) if parts[7] != '-1' else 0
                        }
                        
                        # Handle alternatives if present
                        if len(parts) > 8 and parts[8]:
                            alternatives = []
                            alt_parts = parts[8].split(';')
                            for i in range(0, len(alt_parts), 3):
                                if i+2 < len(alt_parts):
                                    alternatives.append({
                                        'allele': alt_parts[i],
                                        'score': float(alt_parts[i+1]),
                                        'reads': int(alt_parts[i+2])
                                    })
                            gene_result['alternatives'] = alternatives
                        
                        # Classify genes by type
                        if gene in self.classical_genes:
                            results['classical_hla'][gene] = gene_result
                        elif gene in self.nonclassical_genes:
                            results['nonclassical_hla'][gene] = gene_result
                        else:
                            results['other_genes'][gene] = gene_result
            
            # Extract novel variants from VCF file
            vcf_file = self.t1k_dir / f"{sample_id}_t1k_result_allele.vcf"
            if vcf_file.exists():
                with open(vcf_file, 'r') as f:
                    variant_count = 0
                    for line in f:
                        if not line.startswith('#') and line.strip():
                            parts = line.strip().split()
                            if len(parts) >= 5:
                                variant_count += 1
                                full_allele = parts[0]
                                
                                # Extract gene name
                                if '*' in full_allele:
                                    gene_part = full_allele.split('*')[0]
                                    gene = gene_part if gene_part.startswith('HLA-') else f"HLA-{gene_part}"
                                else:
                                    gene = full_allele if full_allele.startswith('HLA-') else f"HLA-{full_allele}"
                                
                                pos = parts[1]
                                ref = parts[3] if len(parts) > 3 else 'N'
                                alt = parts[4] if len(parts) > 4 else 'N'
                                
                                if gene not in results['novel_variants']:
                                    results['novel_variants'][gene] = []
                                
                                results['novel_variants'][gene].append({
                                    'position': pos,
                                    'ref_allele': ref,
                                    'alt_allele': alt,
                                    'variant_type': 'SNP' if len(ref) == len(alt) == 1 else 'INDEL',
                                    'reference_allele': full_allele
                                })
                
                    results['variant_summary'] = {
                        'total_variants': variant_count,
                        'genes_with_variants': len(results['novel_variants'])
                    }
        
        except Exception as e:
            results['errors'].append(f"T1K extraction error: {str(e)}")
        
        return results
    
    def create_results_table(self, sample_list: List[str]) -> pd.DataFrame:
        """Create comprehensive results table for all samples"""
        table_data = []
        
        for sample_id in sample_list:
            print(f"Processing sample: {sample_id}")
            
            # Extract T1K results
            t1k_results = self.extract_t1k_results(sample_id)
            
            # Process all genes (classical + non-classical)
            all_genes_data = {**t1k_results['classical_hla'], 
                            **t1k_results['nonclassical_hla'],
                            **t1k_results['other_genes']}
            
            for gene in self.all_genes:
                row = self.create_gene_row(sample_id, gene, t1k_results, all_genes_data)
                table_data.append(row)
        
        return pd.DataFrame(table_data)
    
    def create_gene_row(self, sample_id: str, gene: str, t1k_results: Dict, all_genes_data: Dict) -> Dict:
        """Create a row for a single gene"""
        
        # Determine gene classification
        if gene in self.classical_genes:
            gene_type = 'Classical'
        elif gene in self.nonclassical_genes:
            gene_type = 'Non-classical'
        else:
            gene_type = 'Other'
        
        row = {
            'Sample_ID': sample_id,
            'Gene': gene,
            'Gene_Type': gene_type,
        }
        
        # T1K results
        if gene in all_genes_data:
            t1k_data = all_genes_data[gene]
            row.update({
                'Allele_1': t1k_data.get('allele_1'),
                'Allele_2': t1k_data.get('allele_2'),
                'Score_1': round(t1k_data.get('score_1', 0), 2),
                'Score_2': round(t1k_data.get('score_2', 0), 2),
                'Reads_1': t1k_data.get('reads_1'),
                'Reads_2': t1k_data.get('reads_2'),
                'Num_Alleles': t1k_data.get('num_alleles'),
                'Total_Reads': t1k_data.get('reads_1', 0) + t1k_data.get('reads_2', 0),
            })
            
            # Alternative alleles information
            if 'alternatives' in t1k_data:
                alt_alleles = [alt['allele'] for alt in t1k_data['alternatives']]
                row['Alternative_Alleles'] = ';'.join(alt_alleles[:3])  # Top 3 alternatives
                row['Alternative_Count'] = len(t1k_data['alternatives'])
            else:
                row['Alternative_Alleles'] = None
                row['Alternative_Count'] = 0
            
            # Typing quality assessment
            avg_score = (t1k_data.get('score_1', 0) + t1k_data.get('score_2', 0)) / 2
            total_reads = t1k_data.get('reads_1', 0) + t1k_data.get('reads_2', 0)
            
            if avg_score > 0.9 and total_reads > 50:
                quality = "High"
            elif avg_score > 0.7 and total_reads > 20:
                quality = "Medium"
            else:
                quality = "Low"
            
            row['Typing_Quality'] = quality
        else:
            row.update({
                'Allele_1': None,
                'Allele_2': None,
                'Score_1': None,
                'Score_2': None,
                'Reads_1': None,
                'Reads_2': None,
                'Num_Alleles': None,
                'Total_Reads': None,
                'Alternative_Alleles': None,
                'Alternative_Count': 0,
                'Typing_Quality': 'Not Typed'
            })
        
        # Variant information
        if gene in t1k_results.get('novel_variants', {}):
            variants = t1k_results['novel_variants'][gene]
            row.update({
                'Novel_Variants_Count': len(variants),
                'Variant_Types': ';'.join(set([v['variant_type'] for v in variants]))
            })
        else:
            row.update({
                'Novel_Variants_Count': 0,
                'Variant_Types': None
            })
        
        # Add normalized 2-field allele names for easier comparison
        row['Allele_1_2field'] = self.normalize_allele_name(row.get('Allele_1'))
        row['Allele_2_2field'] = self.normalize_allele_name(row.get('Allele_2'))
        
        # Errors if any
        if t1k_results['errors']:
            row['Errors'] = ';'.join(t1k_results['errors'])
        else:
            row['Errors'] = None
        
        return row
    
    def normalize_allele_name(self, allele: str, target_resolution: str = "2-field") -> str:
        """Normalize allele names to 2-field resolution"""
        if not allele or allele == '.':
            return None
        
        # Remove HLA- prefix if present
        if allele.startswith('HLA-'):
            allele = allele[4:]
        
        # Split by asterisk and colon
        parts = allele.split('*')
        if len(parts) < 2:
            return allele
        
        gene = parts[0]
        allele_parts = parts[1].split(':')
        
        if target_resolution == "2-field" and len(allele_parts) >= 2:
            return f"{gene}*{allele_parts[0]}:{allele_parts[1]}"
        elif target_resolution == "1-field" and len(allele_parts) >= 1:
            return f"{gene}*{allele_parts[0]}"
        
        return allele
    
    def create_summary_statistics(self, results_df: pd.DataFrame) -> Dict:
        """Generate summary statistics"""
        stats = {
            'overview': {
                'total_genes_analyzed': len(results_df),
                'samples_analyzed': results_df['Sample_ID'].nunique(),
                'classical_genes': len(results_df[results_df['Gene_Type'] == 'Classical']),
                'nonclassical_genes': len(results_df[results_df['Gene_Type'] == 'Non-classical'])
            },
            'typing_success': {
                'successfully_typed': (results_df['Allele_1'].notna()).sum(),
                'not_typed': (results_df['Allele_1'].isna()).sum(),
                'high_quality': (results_df['Typing_Quality'] == 'High').sum(),
                'medium_quality': (results_df['Typing_Quality'] == 'Medium').sum(),
                'low_quality': (results_df['Typing_Quality'] == 'Low').sum()
            },
            'gene_type_analysis': {
                'classical_typed': (results_df[results_df['Gene_Type'] == 'Classical']['Allele_1'].notna()).sum(),
                'nonclassical_typed': (results_df[results_df['Gene_Type'] == 'Non-classical']['Allele_1'].notna()).sum()
            },
            'variant_analysis': {},
            'read_coverage': {}
        }
        
        # Variant analysis
        variant_data = results_df[results_df['Novel_Variants_Count'] > 0]
        if not variant_data.empty:
            stats['variant_analysis'] = {
                'genes_with_variants': len(variant_data),
                'total_variants': variant_data['Novel_Variants_Count'].sum(),
                'average_variants_per_gene': variant_data['Novel_Variants_Count'].mean()
            }
        
        # Read coverage analysis
        coverage_data = results_df[results_df['Total_Reads'].notna()]
        if not coverage_data.empty:
            stats['read_coverage'] = {
                'average_total_reads': coverage_data['Total_Reads'].mean(),
                'median_total_reads': coverage_data['Total_Reads'].median(),
                'min_reads': coverage_data['Total_Reads'].min(),
                'max_reads': coverage_data['Total_Reads'].max()
            }
        
        return stats
    
    def generate_report(self, sample_list: List[str], output_dir: str = ".") -> None:
        """Generate comprehensive T1K extraction report"""
        print("Starting T1K Results Extraction...")
        
        # Create results table
        results_df = self.create_results_table(sample_list)
        
        # Generate statistics
        stats = self.create_summary_statistics(results_df)
        
        # Create output directory
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        # Save main results table
        results_df.to_csv(output_path / "t1k_results_complete.csv", index=False)
        
        # Separate tables for different gene types
        classical_df = results_df[results_df['Gene_Type'] == 'Classical']
        nonclassical_df = results_df[results_df['Gene_Type'] == 'Non-classical']
        
        classical_df.to_csv(output_path / "t1k_results_classical.csv", index=False)
        nonclassical_df.to_csv(output_path / "t1k_results_nonclassical.csv", index=False)
        
        # Immune evasion genes (HLA-E, F, G)
        immune_genes = ['HLA-E', 'HLA-F', 'HLA-G']
        immune_df = results_df[results_df['Gene'].isin(immune_genes)]
        immune_df.to_csv(output_path / "t1k_results_immune_evasion.csv", index=False)
        
        # Variant analysis summary
        variant_df = results_df[results_df['Novel_Variants_Count'] > 0]
        if not variant_df.empty:
            variant_summary = variant_df[['Sample_ID', 'Gene', 'Novel_Variants_Count', 'Variant_Types']]
            variant_summary.to_csv(output_path / "t1k_variant_summary.csv", index=False)
        
        # Summary report
        with open(output_path / "t1k_extraction_summary.txt", 'w') as f:
            f.write("T1K HLA Results Extraction Summary\n")
            f.write("=" * 50 + "\n\n")
            
            # Overview
            f.write("OVERVIEW\n")
            f.write("-" * 20 + "\n")
            f.write(f"Total Genes Analyzed: {stats['overview']['total_genes_analyzed']}\n")
            f.write(f"Samples Analyzed: {stats['overview']['samples_analyzed']}\n")
            f.write(f"Classical Genes: {stats['overview']['classical_genes']}\n")
            f.write(f"Non-classical Genes: {stats['overview']['nonclassical_genes']}\n\n")
            
            # Typing Success
            f.write("TYPING SUCCESS\n")
            f.write("-" * 20 + "\n")
            f.write(f"Successfully Typed: {stats['typing_success']['successfully_typed']}\n")
            f.write(f"Not Typed: {stats['typing_success']['not_typed']}\n")
            f.write(f"High Quality: {stats['typing_success']['high_quality']}\n")
            f.write(f"Medium Quality: {stats['typing_success']['medium_quality']}\n")
            f.write(f"Low Quality: {stats['typing_success']['low_quality']}\n\n")
            
            # Gene Type Analysis
            f.write("GENE TYPE ANALYSIS\n")
            f.write("-" * 25 + "\n")
            f.write(f"Classical Genes Typed: {stats['gene_type_analysis']['classical_typed']}\n")
            f.write(f"Non-classical Genes Typed: {stats['gene_type_analysis']['nonclassical_typed']}\n\n")
            
            # Variant Analysis
            if stats['variant_analysis']:
                f.write("VARIANT ANALYSIS\n")
                f.write("-" * 20 + "\n")
                f.write(f"Genes with Variants: {stats['variant_analysis']['genes_with_variants']}\n")
                f.write(f"Total Variants: {stats['variant_analysis']['total_variants']}\n")
                f.write(f"Average Variants per Gene: {stats['variant_analysis']['average_variants_per_gene']:.1f}\n\n")
            
            # Read Coverage
            if stats['read_coverage']:
                f.write("READ COVERAGE STATISTICS\n")
                f.write("-" * 30 + "\n")
                f.write(f"Average Total Reads: {stats['read_coverage']['average_total_reads']:.1f}\n")
                f.write(f"Median Total Reads: {stats['read_coverage']['median_total_reads']:.1f}\n")
                f.write(f"Min Reads: {stats['read_coverage']['min_reads']}\n")
                f.write(f"Max Reads: {stats['read_coverage']['max_reads']}\n\n")
        
        print(f"\nT1K extraction complete. Results saved to {output_path}")
        print("\nGenerated files:")
        print("  - t1k_results_complete.csv (All results)")
        print("  - t1k_results_classical.csv (Classical genes only)")
        print("  - t1k_results_nonclassical.csv (Non-classical genes only)")
        print("  - t1k_results_immune_evasion.csv (HLA-E/F/G)")
        print("  - t1k_variant_summary.csv (Novel variants)")
        print("  - t1k_extraction_summary.txt (Summary statistics)")


def main():
    """Main execution function"""
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python t1k_results_extraction.py <t1k_results_directory> [sample1 sample2 ...]")
        print("Example: python t1k_results_extraction.py ~/hla_analysis/results_t1k 001_HypExomexBRCA1416_run19")
        sys.exit(1)
    
    base_dir = sys.argv[1]
    
    # Get sample list from command line or auto-detect
    if len(sys.argv) > 2:
        sample_list = sys.argv[2:]
    else:
        # Auto-detect samples from T1K genotype files
        t1k_dir = Path(base_dir)
        if not t1k_dir.exists():
            print(f"Error: Directory {base_dir} does not exist")
            sys.exit(1)
        
        t1k_files = list(t1k_dir.glob('*_t1k_result_genotype.tsv'))
        sample_list = [f.name.replace('_t1k_result_genotype.tsv', '') for f in t1k_files]
        sample_list = sorted(list(set(sample_list)))
    
    if not sample_list:
        print("No T1K samples found. Please specify sample names manually.")
        sys.exit(1)
    
    print(f"Found samples: {sample_list}")
    
    # Create extractor and run analysis
    extractor = T1KExtractor(base_dir)
    extractor.generate_report(sample_list, output_dir=f"{base_dir}/t1k_extraction_results")


if __name__ == "__main__":
    main()
