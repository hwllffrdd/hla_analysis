#!/usr/bin/env python3
"""
Enhanced HLA Tools Results Comparison Script for Oncogenetics
Extracts comprehensive data from SpecHLA, T1K, and HLA-HD including:
- Non-classical HLA genes (HLA-E, F, G, H, J, K, L) for immune evasion analysis
- Enhanced LOH information with allelic frequencies
- Novel variant detection from VCF files
- Extended gene coverage comparison
"""

import pandas as pd
import os
import re
from pathlib import Path
from typing import Dict, List, Tuple, Optional

class EnhancedHLAExtractor:
    def __init__(self, base_dir: str):
        self.base_dir = Path(base_dir)
        self.spechla_dir = self.base_dir / "results_spechla"
        self.t1k_dir = self.base_dir / "results_t1k" 
        self.hlahd_dir = self.base_dir / "results_hlahd"
        
        # Classical HLA genes for comparison
        self.classical_genes = [
            'HLA-A', 'HLA-B', 'HLA-C',           # Class I
            'HLA-DRB1', 'HLA-DQB1', 'HLA-DPB1', # Class II major
            'HLA-DQA1', 'HLA-DPA1'              # Class II alpha
        ]
        
        # Non-classical HLA genes of oncological interest
        self.nonclassical_genes = [
            'HLA-E', 'HLA-F', 'HLA-G',  # Immune evasion genes
            'HLA-H', 'HLA-J', 'HLA-K', 'HLA-L'  # Potential biomarkers
        ]
        
        # Dynamically discover all HLA-HD genes from any sample
        self.hlahd_additional_genes = self._discover_hlahd_genes()
        
        # All genes to analyze
        self.all_genes = self.classical_genes + self.nonclassical_genes + self.hlahd_additional_genes
    
    def _discover_hlahd_genes(self) -> List[str]:
        """Discover all HLA genes that HLA-HD can analyze by checking final.result.txt files"""
        all_hlahd_genes = set()
        
        # Check all HLA-HD result directories
        if self.hlahd_dir.exists():
            for sample_dir in self.hlahd_dir.iterdir():
                if sample_dir.is_dir():
                    result_file = sample_dir / sample_dir.name.replace('_hlahd', '') / "result" / f"{sample_dir.name.replace('_hlahd', '')}_final.result.txt"
                    if result_file.exists():
                        with open(result_file, 'r') as f:
                            for line in f:
                                parts = line.strip().split('\t')
                                if len(parts) >= 2:
                                    gene = parts[0]
                                    # Convert to HLA- format and add to set
                                    full_gene_name = f"HLA-{gene}" if not gene.startswith('HLA') else gene
                                    all_hlahd_genes.add(full_gene_name)
        
        # Remove genes already in classical or nonclassical lists
        existing_genes = set(self.classical_genes + self.nonclassical_genes)
        additional_genes = list(all_hlahd_genes - existing_genes)
        additional_genes.sort()  # Sort for consistent ordering
        
        return additional_genes
    
    def extract_enhanced_spechla_results(self, sample_id: str) -> Dict:
        """Extract comprehensive SpecHLA results including enhanced LOH information"""
        sample_dir = self.spechla_dir / sample_id
        
        results = {
            'sample_id': sample_id,
            'tool': 'SpecHLA',
            'resolution': '4-field',
            'hla_typing': {},
            'g_group': {},
            'quality_scores': {},
            'loh_status': {},
            'allelic_frequencies': {},
            'het_variants_count': {},
            'errors': []
        }
        
        try:
            # Main HLA typing results
            main_result_file = sample_dir / "hla.result.txt"
            if main_result_file.exists():
                with open(main_result_file, 'r') as f:
                    lines = f.readlines()
                    if len(lines) >= 3:
                        header = lines[1].strip().split('\t')
                        data = lines[2].strip().split('\t')
                        
                        for i, gene_field in enumerate(header[1:], 1):
                            if i < len(data):
                                gene_parts = gene_field.split('_')
                                if len(gene_parts) >= 3:
                                    gene = f"{gene_parts[0]}_{gene_parts[1]}"
                                    allele_num = gene_parts[2]
                                    
                                    if gene not in results['hla_typing']:
                                        results['hla_typing'][gene] = {}
                                    results['hla_typing'][gene][f'allele_{allele_num}'] = data[i]
            
            # G-group results
            ggroup_file = sample_dir / "hla.result.g.group.txt"
            if ggroup_file.exists():
                with open(ggroup_file, 'r') as f:
                    lines = f.readlines()
                    if len(lines) >= 3:
                        header = lines[1].strip().split('\t')
                        data = lines[2].strip().split('\t')
                        
                        for i, gene_field in enumerate(header[1:], 1):
                            if i < len(data):
                                gene_parts = gene_field.split('_')
                                if len(gene_parts) >= 3:
                                    gene = f"{gene_parts[0]}_{gene_parts[1]}"
                                    allele_num = gene_parts[2]
                                    
                                    if gene not in results['g_group']:
                                        results['g_group'][gene] = {}
                                    results['g_group'][gene][f'allele_{allele_num}'] = data[i]
            
            # Enhanced LOH detection with allelic frequencies
            for gene in ['HLA_A', 'HLA_B', 'HLA_C', 'HLA_DRB1', 'HLA_DQB1', 'HLA_DPB1', 'HLA_DQA1', 'HLA_DPA1']:
                freq_file = sample_dir / f"{gene}_freq.txt"
                if freq_file.exists():
                    with open(freq_file, 'r') as f:
                        lines = f.readlines()
                        frequencies = []
                        het_variants = 0
                        
                        for line in lines:
                            if line.startswith('#') and 'heterozygous variant' in line:
                                try:
                                    het_variants = int(re.search(r'(\d+)', line).group(1))
                                except:
                                    het_variants = 0
                            elif not line.startswith('#') and line.strip():
                                parts = line.strip().split()
                                if len(parts) >= 2:
                                    try:
                                        frequencies.append(float(parts[1]))
                                    except:
                                        continue
                        
                        if frequencies:
                            # Enhanced LOH analysis
                            freq1, freq2 = frequencies[0], frequencies[1] if len(frequencies) > 1 else 0.0
                            freq_diff = abs(freq1 - freq2)
                            
                            # More nuanced LOH classification
                            if freq_diff > 0.4:
                                loh_level = "Severe LOH"
                            elif freq_diff > 0.3:
                                loh_level = "Moderate LOH"
                            elif freq_diff > 0.2:
                                loh_level = "Mild LOH"
                            else:
                                loh_level = "No LOH"
                            
                            results['loh_status'][gene] = {
                                'loh_detected': freq_diff > 0.3,
                                'loh_level': loh_level,
                                'frequency_difference': freq_diff
                            }
                            results['allelic_frequencies'][gene] = {
                                'allele1_freq': freq1,
                                'allele2_freq': freq2,
                                'frequencies_list': frequencies
                            }
                            results['het_variants_count'][gene] = het_variants
        
        except Exception as e:
            results['errors'].append(f"Enhanced SpecHLA extraction error: {str(e)}")
        
        return results
    
    def extract_extended_t1k_results(self, sample_id: str) -> Dict:
        """Extract T1K results with extended gene coverage and variant detection"""
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
            # Main genotype results with extended gene classification
            genotype_file = self.t1k_dir / f"{sample_id}_t1k_result_genotype.tsv"
            if genotype_file.exists():
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
                            # VCF files from T1K are space-separated, not tab-separated
                            # Using split() without arguments handles multiple whitespace characters
                            parts = line.strip().split()
                            if len(parts) >= 5:
                                variant_count += 1
                                # Extract variant information
                                full_allele = parts[0]  # e.g., "HLA-J*01:01:01:04"
                    
                                # Extract just the gene name from the full allele
                                # Handle both "HLA-J*01:01:01:04" and "J*01:01:01:04" formats
                                if '*' in full_allele:
                                    gene_part = full_allele.split('*')[0]
                                    # Ensure it has HLA- prefix
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
                                    'reference_allele': full_allele  # Keep the full allele for reference
                                })
        
                    # Store summary statistics
                    results['variant_summary'] = {
                        'total_variants': variant_count,
                        'genes_with_variants': len(results['novel_variants'])
                    }
        
        except Exception as e:
            results['errors'].append(f"Extended T1K extraction error: {str(e)}")
        
        return results
    
    def extract_comprehensive_hlahd_results(self, sample_id: str) -> Dict:
        """Extract comprehensive HLA-HD results including non-classical genes"""
        results = {
            'sample_id': sample_id,
            'tool': 'HLA-HD',
            'resolution': '3-field',
            'classical_hla': {},
            'nonclassical_hla': {},
            'other_hla_genes': {},
            'detailed_calls': {},
            'read_counts': {},
            'not_typed_genes': [],
            'gene_coverage_stats': {},
            'errors': []
        }
        
        try:
            # Main result file
            result_dir = self.hlahd_dir / f"{sample_id}_hlahd" / sample_id / "result"
            final_result_file = result_dir / f"{sample_id}_final.result.txt"
            
            if final_result_file.exists():
                with open(final_result_file, 'r') as f:
                    for line in f:
                        parts = line.strip().split('\t')
                        if len(parts) >= 3:
                            gene = parts[0]
                            allele_1 = parts[1] if parts[1] != 'Not typed' else None
                            allele_2 = parts[2] if parts[2] != 'Not typed' and parts[2] != '-' else None
                            
                            # Handle additional alleles for some genes (like F)
                            additional_alleles = []
                            if len(parts) > 3:
                                for extra_part in parts[3:]:
                                    if extra_part and extra_part not in ['Not typed', '-']:
                                        additional_alleles.append(extra_part)
                            
                            gene_result = {
                                'allele_1': allele_1,
                                'allele_2': allele_2,
                                'additional_alleles': additional_alleles,
                                'homozygous': (allele_2 == '-' and allele_1 is not None),
                                'typing_status': 'typed' if allele_1 else 'not_typed'
                            }
                            
                            if allele_1 is None:
                                results['not_typed_genes'].append(gene)
                            
                            # Classify genes
                            full_gene_name = f"HLA-{gene}" if not gene.startswith('HLA') else gene
                            
                            if full_gene_name in self.classical_genes:
                                results['classical_hla'][full_gene_name] = gene_result
                            elif full_gene_name in self.nonclassical_genes:
                                results['nonclassical_hla'][full_gene_name] = gene_result
                            else:
                                results['other_hla_genes'][full_gene_name] = gene_result
            
            # Extract detailed results and read counts for all genes
            all_possible_genes = ['A', 'B', 'C', 'DRB1', 'DQA1', 'DQB1', 'DPA1', 'DPB1', 
                                'E', 'F', 'G', 'H', 'J', 'K', 'L', 'DRB2', 'DRB3', 'DRB4', 'DRB5']
            
            for gene in all_possible_genes:
                est_file = result_dir / f"{sample_id}_{gene}.est.txt"
                read_file = result_dir / f"{sample_id}_{gene}.read.txt"
                
                if est_file.exists():
                    with open(est_file, 'r') as f:
                        lines = f.readlines()
                        if len(lines) >= 3:
                            data_line = lines[2].strip().split('\t')
                            if len(data_line) >= 4:
                                results['detailed_calls'][f'HLA-{gene}'] = {
                                    'candidates_1': data_line[0].split(',') if data_line[0] else [],
                                    'candidates_2': data_line[1].split(',') if data_line[1] else [],
                                    'scores_1': data_line[2],
                                    'scores_2': data_line[3]
                                }
                
                if read_file.exists():
                    with open(read_file, 'r') as f:
                        lines = f.readlines()
                        gene_reads = {}
                        total_reads = 0
                        
                        for line in lines:
                            if line.startswith('HLA-'):
                                parts = line.strip().split('\t')
                                if len(parts) >= 2:
                                    allele = parts[0]
                                    read_count = int(parts[1])
                                    gene_reads[allele] = read_count
                                    total_reads += read_count
                        
                        if gene_reads:
                            results['read_counts'][f'HLA-{gene}'] = gene_reads
                            results['gene_coverage_stats'][f'HLA-{gene}'] = {
                                'total_reads': total_reads,
                                'unique_alleles': len(gene_reads),
                                'coverage_quality': 'high' if total_reads > 50 else 'medium' if total_reads > 20 else 'low'
                            }
        
        except Exception as e:
            results['errors'].append(f"Comprehensive HLA-HD extraction error: {str(e)}")
        
        return results
    
    def create_comprehensive_comparison_table(self, sample_list: List[str]) -> pd.DataFrame:
        """Create comprehensive comparison table including classical and non-classical genes"""
        comparison_data = []
        
        for sample_id in sample_list:
            print(f"Processing sample: {sample_id}")
            
            # Extract enhanced results from all three tools
            print("  Extracting enhanced SpecHLA results...")
            spechla_results = self.extract_enhanced_spechla_results(sample_id)
            
            print("  Extracting extended T1K results...")
            t1k_results = self.extract_extended_t1k_results(sample_id)
            
            print("  Extracting comprehensive HLA-HD results...")
            hlahd_results = self.extract_comprehensive_hlahd_results(sample_id)
            
            # Process all genes (classical + non-classical)
            for gene in self.all_genes:
                row = self.create_gene_comparison_row(sample_id, gene, spechla_results, t1k_results, hlahd_results)
                comparison_data.append(row)
        
        return pd.DataFrame(comparison_data)
    
    def create_gene_comparison_row(self, sample_id: str, gene: str, spechla_results: Dict, 
                                   t1k_results: Dict, hlahd_results: Dict) -> Dict:
        """Create a comprehensive comparison row for a single gene"""
        
        # Determine gene classification
        if gene in self.classical_genes:
            gene_type = 'Classical'
        elif gene in self.nonclassical_genes:
            gene_type = 'Non-classical'
        else:
            gene_type = 'HLA-HD Additional'
        
        row = {
            'Sample_ID': sample_id,
            'Gene': gene,
            'Gene_Type': gene_type,
        }
        
        # SpecHLA results (only for classical genes)
        spechla_gene = gene.replace('-', '_')
        if gene_type == 'Classical':
            if spechla_gene in spechla_results['hla_typing']:
                spechla_data = spechla_results['hla_typing'][spechla_gene]
                row.update({
                    'SpecHLA_Allele1': spechla_data.get('allele_1'),
                    'SpecHLA_Allele2': spechla_data.get('allele_2'),
                    'SpecHLA_Resolution': '4-field'
                })
                
                # Enhanced LOH information
                if spechla_gene in spechla_results['loh_status']:
                    loh_data = spechla_results['loh_status'][spechla_gene]
                    freq_data = spechla_results['allelic_frequencies'].get(spechla_gene, {})
                    row.update({
                        'SpecHLA_LOH_Level': loh_data.get('loh_level'),
                        'SpecHLA_Freq_Diff': round(loh_data.get('frequency_difference', 0), 3),
                        'SpecHLA_Allele1_Freq': round(freq_data.get('allele1_freq', 0), 3),
                        'SpecHLA_Allele2_Freq': round(freq_data.get('allele2_freq', 0), 3),
                        'SpecHLA_Het_Variants': spechla_results['het_variants_count'].get(spechla_gene, 0)
                    })
                
                # G-group information
                if spechla_gene in spechla_results['g_group']:
                    g_data = spechla_results['g_group'][spechla_gene]
                    row.update({
                        'SpecHLA_GGroup1': g_data.get('allele_1'),
                        'SpecHLA_GGroup2': g_data.get('allele_2')
                    })
            else:
                row.update({
                    'SpecHLA_Allele1': None, 'SpecHLA_Allele2': None,
                    'SpecHLA_Resolution': 'N/A', 'SpecHLA_LOH_Level': None,
                    'SpecHLA_Freq_Diff': None, 'SpecHLA_Het_Variants': None
                })
        else:
            # Non-classical genes not covered by SpecHLA
            row.update({
                'SpecHLA_Allele1': 'Not covered', 'SpecHLA_Allele2': 'Not covered',
                'SpecHLA_Resolution': 'N/A', 'SpecHLA_LOH_Level': 'N/A'
            })
        
        # T1K results (covers both classical and non-classical)
        t1k_source = t1k_results['classical_hla'] if gene_type == 'Classical' else t1k_results['nonclassical_hla']
        
        if gene in t1k_source:
            t1k_data = t1k_source[gene]
            row.update({
                'T1K_Allele1': t1k_data.get('allele_1'),
                'T1K_Allele2': t1k_data.get('allele_2'),
                'T1K_Score1': round(t1k_data.get('score_1', 0), 2),
                'T1K_Score2': round(t1k_data.get('score_2', 0), 2),
                'T1K_Reads1': t1k_data.get('reads_1'),
                'T1K_Reads2': t1k_data.get('reads_2'),
                'T1K_Num_Alleles': t1k_data.get('num_alleles'),
                'T1K_Resolution': '4-field'
            })
            
            # Alternative alleles information
            if 'alternatives' in t1k_data:
                alt_alleles = [alt['allele'] for alt in t1k_data['alternatives']]
                row['T1K_Alternative_Alleles'] = ';'.join(alt_alleles[:3])  # Top 3 alternatives
            else:
                row['T1K_Alternative_Alleles'] = None
        else:
            row.update({
                'T1K_Allele1': None, 'T1K_Allele2': None,
                'T1K_Score1': None, 'T1K_Score2': None,
                'T1K_Reads1': None, 'T1K_Reads2': None,
                'T1K_Num_Alleles': None, 'T1K_Resolution': 'N/A',
                'T1K_Alternative_Alleles': None
            })
        
        # HLA-HD results (covers classical, non-classical, and additional genes)
        if gene_type == 'Classical':
            hlahd_source = hlahd_results['classical_hla']
        elif gene_type == 'Non-classical':
            hlahd_source = hlahd_results['nonclassical_hla']
        else:
            # Additional HLA-HD genes
            hlahd_source = hlahd_results['other_hla_genes']
        
        if gene in hlahd_source:
            hlahd_data = hlahd_source[gene]
            row.update({
                'HLAHD_Allele1': hlahd_data.get('allele_1'),
                'HLAHD_Allele2': hlahd_data.get('allele_2'),
                'HLAHD_Additional_Alleles': ';'.join(hlahd_data.get('additional_alleles', [])),
                'HLAHD_Homozygous': hlahd_data.get('homozygous', False),
                'HLAHD_Typing_Status': hlahd_data.get('typing_status'),
                'HLAHD_Resolution': '3-field'
            })
            
            # Coverage statistics
            if gene in hlahd_results['gene_coverage_stats']:
                coverage_data = hlahd_results['gene_coverage_stats'][gene]
                row.update({
                    'HLAHD_Total_Reads': coverage_data.get('total_reads'),
                    'HLAHD_Coverage_Quality': coverage_data.get('coverage_quality'),
                    'HLAHD_Unique_Alleles': coverage_data.get('unique_alleles')
                })
        else:
            row.update({
                'HLAHD_Allele1': None, 'HLAHD_Allele2': None,
                'HLAHD_Homozygous': None, 'HLAHD_Typing_Status': 'not_analyzed',
                'HLAHD_Resolution': 'N/A', 'HLAHD_Total_Reads': None,
                'HLAHD_Coverage_Quality': None
            })
        
        # Variant information from T1K
        if gene in t1k_results.get('novel_variants', {}):
            variants = t1k_results['novel_variants'][gene]
            row.update({
                'T1K_Novel_Variants_Count': len(variants),
                'T1K_Variant_Types': ';'.join(set([v['variant_type'] for v in variants]))
            })
        else:
            row.update({
                'T1K_Novel_Variants_Count': 0,
                'T1K_Variant_Types': None
            })
        
        # Agreement analysis (adapted for different gene types)
        if gene_type == 'Classical':
            # All three tools
            row['Agreement_All_Tools'] = self.check_tool_agreement_enhanced(
                [row.get('SpecHLA_Allele1'), row.get('SpecHLA_Allele2')],
                [row.get('T1K_Allele1'), row.get('T1K_Allele2')],
                [row.get('HLAHD_Allele1'), row.get('HLAHD_Allele2')]
            )
            
            row['Agreement_SpecHLA_T1K'] = self.check_pairwise_agreement_enhanced(
                [row.get('SpecHLA_Allele1'), row.get('SpecHLA_Allele2')],
                [row.get('T1K_Allele1'), row.get('T1K_Allele2')]
            )
            
            row['Agreement_SpecHLA_HLAHD'] = self.check_pairwise_agreement_enhanced(
                [row.get('SpecHLA_Allele1'), row.get('SpecHLA_Allele2')],
                [row.get('HLAHD_Allele1'), row.get('HLAHD_Allele2')]
            )
        else:
            # Non-classical genes: only T1K vs HLA-HD
            row['Agreement_All_Tools'] = 'N/A (SpecHLA not applicable)'
            row['Agreement_SpecHLA_T1K'] = 'N/A (SpecHLA not applicable)'
            row['Agreement_SpecHLA_HLAHD'] = 'N/A (SpecHLA not applicable)'
        
        # T1K vs HLA-HD agreement (for both gene types)
        row['Agreement_T1K_HLAHD'] = self.check_pairwise_agreement_enhanced(
            [row.get('T1K_Allele1'), row.get('T1K_Allele2')],
            [row.get('HLAHD_Allele1'), row.get('HLAHD_Allele2')]
        )
        
        # Add normalized 2-field comparison for easier analysis
        row['SpecHLA_2field_1'] = self.normalize_allele_name(row.get('SpecHLA_Allele1'))
        row['SpecHLA_2field_2'] = self.normalize_allele_name(row.get('SpecHLA_Allele2'))
        row['T1K_2field_1'] = self.normalize_allele_name(row.get('T1K_Allele1'))
        row['T1K_2field_2'] = self.normalize_allele_name(row.get('T1K_Allele2'))
        row['HLAHD_2field_1'] = self.normalize_allele_name(row.get('HLAHD_Allele1'))
        row['HLAHD_2field_2'] = self.normalize_allele_name(row.get('HLAHD_Allele2'))
        
        return row
    
    def check_tool_agreement_enhanced(self, alleles1: List, alleles2: List, alleles3: List) -> str:
        """Enhanced agreement checking with better handling of None values"""
        # Remove None values and sort for comparison
        set1 = set([self.normalize_allele_name(a) for a in alleles1 if a is not None and a != 'Not covered'])
        set2 = set([self.normalize_allele_name(a) for a in alleles2 if a is not None])
        set3 = set([self.normalize_allele_name(a) for a in alleles3 if a is not None])
        
        # Remove None from normalized results
        set1 = set([a for a in set1 if a is not None])
        set2 = set([a for a in set2 if a is not None])
        set3 = set([a for a in set3 if a is not None])
        
        tools_with_data = sum([bool(s) for s in [set1, set2, set3]])
        
        if tools_with_data >= 2:
            if set1 and set2 and set3:
                # All three tools have data
                if set1 == set2 == set3:
                    return "Full Agreement (3 tools)"
                elif (set1 == set2) or (set1 == set3) or (set2 == set3):
                    return "Partial Agreement (2/3 tools agree)"
                elif len(set1.intersection(set2).intersection(set3)) > 0:
                    return "Partial Agreement (some overlap)"
                else:
                    return "No Agreement (3 tools)"
            else:
                # Only two tools have data
                available_sets = [s for s in [set1, set2, set3] if s]
                if len(available_sets) == 2:
                    if available_sets[0] == available_sets[1]:
                        return "Full Agreement (2 tools)"
                    elif len(available_sets[0].intersection(available_sets[1])) > 0:
                        return "Partial Agreement (2 tools)"
                    else:
                        return "No Agreement (2 tools)"
                else:
                    return "Single Tool Only"
        else:
            return "Insufficient Data"
    
    def check_pairwise_agreement_enhanced(self, alleles1: List, alleles2: List) -> str:
        """Enhanced pairwise agreement checking"""
        set1 = set([self.normalize_allele_name(a) for a in alleles1 if a is not None and a != 'Not covered'])
        set2 = set([self.normalize_allele_name(a) for a in alleles2 if a is not None])
        
        # Remove None from normalized results
        set1 = set([a for a in set1 if a is not None])
        set2 = set([a for a in set2 if a is not None])
        
        if set1 and set2:
            if set1 == set2:
                return "Full Agreement"
            elif len(set1.intersection(set2)) > 0:
                return "Partial Agreement"
            else:
                return "No Agreement"
        elif set1 or set2:
            return "One Tool Only"
        else:
            return "No Data"
    
    def normalize_allele_name(self, allele: str, target_resolution: str = "2-field") -> str:
        """Normalize allele names for comparison"""
        if not allele or allele in ['Not typed', 'Not covered', '.']:
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
    
    def create_enhanced_summary_statistics(self, comparison_df: pd.DataFrame) -> Dict:
        """Generate enhanced summary statistics including non-classical genes analysis"""
        stats = {
            'overview': {
                'total_comparisons': len(comparison_df),
                'samples_analyzed': comparison_df['Sample_ID'].nunique(),
                'classical_genes_analyzed': len(comparison_df[comparison_df['Gene_Type'] == 'Classical']),
                'nonclassical_genes_analyzed': len(comparison_df[comparison_df['Gene_Type'] == 'Non-classical'])
            },
            'tool_coverage': {},
            'agreement_stats': {},
            'loh_analysis': {},
            'variant_analysis': {},
            'gene_type_analysis': {}
        }
        
        # Tool coverage analysis
        classical_df = comparison_df[comparison_df['Gene_Type'] == 'Classical']
        nonclassical_df = comparison_df[comparison_df['Gene_Type'] == 'Non-classical']
        
        stats['tool_coverage'] = {
            'classical_genes': {
                'spechla_coverage': (classical_df['SpecHLA_Allele1'].notna() & 
                                   (classical_df['SpecHLA_Allele1'] != 'Not covered')).sum(),
                't1k_coverage': classical_df['T1K_Allele1'].notna().sum(),
                'hlahd_coverage': classical_df['HLAHD_Allele1'].notna().sum()
            },
            'nonclassical_genes': {
                'spechla_coverage': 0,  # SpecHLA doesn't cover non-classical
                't1k_coverage': nonclassical_df['T1K_Allele1'].notna().sum(),
                'hlahd_coverage': nonclassical_df['HLAHD_Allele1'].notna().sum()
            }
        }
        
        # Agreement analysis
        all_agreement = comparison_df['Agreement_All_Tools']
        stats['agreement_stats'] = {
            'classical_genes': {
                'full_agreement': (classical_df['Agreement_All_Tools'].str.contains('Full Agreement', na=False)).sum(),
                'partial_agreement': (classical_df['Agreement_All_Tools'].str.contains('Partial Agreement', na=False)).sum(),
                'no_agreement': (classical_df['Agreement_All_Tools'].str.contains('No Agreement', na=False)).sum()
            },
            'nonclassical_genes': {
                't1k_hlahd_agreement': (nonclassical_df['Agreement_T1K_HLAHD'] == 'Full Agreement').sum(),
                't1k_hlahd_partial': (nonclassical_df['Agreement_T1K_HLAHD'] == 'Partial Agreement').sum(),
                't1k_hlahd_no_agreement': (nonclassical_df['Agreement_T1K_HLAHD'] == 'No Agreement').sum()
            }
        }
        
        # LOH analysis (SpecHLA only, classical genes)
        loh_data = classical_df[classical_df['SpecHLA_LOH_Level'].notna()]
        if not loh_data.empty:
            stats['loh_analysis'] = {
                'total_genes_analyzed': len(loh_data),
                'severe_loh': (loh_data['SpecHLA_LOH_Level'] == 'Severe LOH').sum(),
                'moderate_loh': (loh_data['SpecHLA_LOH_Level'] == 'Moderate LOH').sum(),
                'mild_loh': (loh_data['SpecHLA_LOH_Level'] == 'Mild LOH').sum(),
                'no_loh': (loh_data['SpecHLA_LOH_Level'] == 'No LOH').sum(),
                'average_freq_difference': loh_data['SpecHLA_Freq_Diff'].mean(),
                'average_het_variants': loh_data['SpecHLA_Het_Variants'].mean()
            }
        
        # Variant analysis (T1K)
        variant_data = comparison_df[comparison_df['T1K_Novel_Variants_Count'].notna()]
        if not variant_data.empty:
            stats['variant_analysis'] = {
                'genes_with_variants': (variant_data['T1K_Novel_Variants_Count'] > 0).sum(),
                'total_variants_detected': variant_data['T1K_Novel_Variants_Count'].sum(),
                'average_variants_per_gene': variant_data['T1K_Novel_Variants_Count'].mean(),
                'variant_types_detected': set([vt for vt_string in variant_data['T1K_Variant_Types'].dropna() 
                                             for vt in vt_string.split(';')])
            }
        
        # Gene type analysis
        stats['gene_type_analysis'] = {
            'classical_genes_successfully_typed': {
                'spechla': (classical_df['SpecHLA_Allele1'].notna() & 
                           (classical_df['SpecHLA_Allele1'] != 'Not covered')).sum(),
                't1k': classical_df['T1K_Allele1'].notna().sum(),
                'hlahd': classical_df['HLAHD_Allele1'].notna().sum()
            },
            'nonclassical_genes_successfully_typed': {
                't1k': nonclassical_df['T1K_Allele1'].notna().sum(),
                'hlahd': nonclassical_df['HLAHD_Allele1'].notna().sum()
            },
            'immune_evasion_genes_detected': {
                'hla_e': len(comparison_df[(comparison_df['Gene'] == 'HLA-E') & 
                                         (comparison_df['T1K_Allele1'].notna() | comparison_df['HLAHD_Allele1'].notna())]),
                'hla_f': len(comparison_df[(comparison_df['Gene'] == 'HLA-F') & 
                                         (comparison_df['T1K_Allele1'].notna() | comparison_df['HLAHD_Allele1'].notna())]),
                'hla_g': len(comparison_df[(comparison_df['Gene'] == 'HLA-G') & 
                                         (comparison_df['T1K_Allele1'].notna() | comparison_df['HLAHD_Allele1'].notna())])
            }
        }
        
        return stats
    
    def create_oncogenetics_focused_report(self, comparison_df: pd.DataFrame, stats: Dict, output_dir: str) -> None:
        """Create oncogenetics-focused analysis report"""
        output_path = Path(output_dir)
        
        # Main enhanced comparison table
        comparison_df.to_csv(output_path / "enhanced_hla_comparison.csv", index=False)
        
        # Separate tables for different gene types
        classical_df = comparison_df[comparison_df['Gene_Type'] == 'Classical']
        nonclassical_df = comparison_df[comparison_df['Gene_Type'] == 'Non-classical']
        
        classical_df.to_csv(output_path / "classical_hla_comparison.csv", index=False)
        nonclassical_df.to_csv(output_path / "nonclassical_hla_comparison.csv", index=False)
        
        # LOH-focused analysis (SpecHLA classical genes only)
        loh_df = classical_df[classical_df['SpecHLA_LOH_Level'].notna()].copy()
        if not loh_df.empty:
            loh_summary = loh_df[['Sample_ID', 'Gene', 'SpecHLA_LOH_Level', 'SpecHLA_Freq_Diff', 
                                'SpecHLA_Allele1_Freq', 'SpecHLA_Allele2_Freq', 'SpecHLA_Het_Variants']]
            loh_summary.to_csv(output_path / "loh_analysis_summary.csv", index=False)
        
        # Immune evasion genes summary
        immune_evasion_genes = ['HLA-E', 'HLA-F', 'HLA-G']
        immune_df = comparison_df[comparison_df['Gene'].isin(immune_evasion_genes)]
        immune_df.to_csv(output_path / "immune_evasion_genes.csv", index=False)
        
        # Variant analysis summary
        variant_df = comparison_df[comparison_df['T1K_Novel_Variants_Count'] > 0]
        if not variant_df.empty:
            variant_summary = variant_df[['Sample_ID', 'Gene', 'T1K_Novel_Variants_Count', 'T1K_Variant_Types']]
            variant_summary.to_csv(output_path / "variant_analysis_summary.csv", index=False)
        
        # Enhanced summary report
        with open(output_path / "enhanced_comparison_summary.txt", 'w') as f:
            f.write("Enhanced HLA Tools Comparison Summary for Oncogenetics\n")
            f.write("=" * 60 + "\n\n")
            
            # Overview
            f.write("OVERVIEW\n")
            f.write("-" * 20 + "\n")
            f.write(f"Total Comparisons: {stats['overview']['total_comparisons']}\n")
            f.write(f"Samples Analyzed: {stats['overview']['samples_analyzed']}\n")
            f.write(f"Classical Genes Analyzed: {stats['overview']['classical_genes_analyzed']}\n")
            f.write(f"Non-classical Genes Analyzed: {stats['overview']['nonclassical_genes_analyzed']}\n\n")
            
            # Tool Coverage
            f.write("TOOL COVERAGE ANALYSIS\n")
            f.write("-" * 30 + "\n")
            f.write("Classical HLA Genes:\n")
            f.write(f"  SpecHLA Coverage: {stats['tool_coverage']['classical_genes']['spechla_coverage']}\n")
            f.write(f"  T1K Coverage: {stats['tool_coverage']['classical_genes']['t1k_coverage']}\n")
            f.write(f"  HLA-HD Coverage: {stats['tool_coverage']['classical_genes']['hlahd_coverage']}\n\n")
            
            f.write("Non-classical HLA Genes:\n")
            f.write(f"  SpecHLA Coverage: {stats['tool_coverage']['nonclassical_genes']['spechla_coverage']} (Not applicable)\n")
            f.write(f"  T1K Coverage: {stats['tool_coverage']['nonclassical_genes']['t1k_coverage']}\n")
            f.write(f"  HLA-HD Coverage: {stats['tool_coverage']['nonclassical_genes']['hlahd_coverage']}\n\n")
            
            # Agreement Analysis
            f.write("AGREEMENT ANALYSIS\n")
            f.write("-" * 25 + "\n")
            f.write("Classical Genes (3-tool comparison):\n")
            f.write(f"  Full Agreement: {stats['agreement_stats']['classical_genes']['full_agreement']}\n")
            f.write(f"  Partial Agreement: {stats['agreement_stats']['classical_genes']['partial_agreement']}\n")
            f.write(f"  No Agreement: {stats['agreement_stats']['classical_genes']['no_agreement']}\n\n")
            
            f.write("Non-classical Genes (T1K vs HLA-HD):\n")
            f.write(f"  Full Agreement: {stats['agreement_stats']['nonclassical_genes']['t1k_hlahd_agreement']}\n")
            f.write(f"  Partial Agreement: {stats['agreement_stats']['nonclassical_genes']['t1k_hlahd_partial']}\n")
            f.write(f"  No Agreement: {stats['agreement_stats']['nonclassical_genes']['t1k_hlahd_no_agreement']}\n\n")
            
            # LOH Analysis
            if 'loh_analysis' in stats:
                f.write("LOH ANALYSIS (SpecHLA - Classical Genes Only)\n")
                f.write("-" * 45 + "\n")
                f.write(f"Genes Analyzed for LOH: {stats['loh_analysis']['total_genes_analyzed']}\n")
                f.write(f"Severe LOH Events: {stats['loh_analysis']['severe_loh']}\n")
                f.write(f"Moderate LOH Events: {stats['loh_analysis']['moderate_loh']}\n")
                f.write(f"Mild LOH Events: {stats['loh_analysis']['mild_loh']}\n")
                f.write(f"No LOH Detected: {stats['loh_analysis']['no_loh']}\n")
                f.write(f"Average Frequency Difference: {stats['loh_analysis']['average_freq_difference']:.3f}\n")
                f.write(f"Average Heterozygous Variants: {stats['loh_analysis']['average_het_variants']:.1f}\n\n")
            
            # Variant Analysis
            if 'variant_analysis' in stats:
                f.write("VARIANT ANALYSIS (T1K)\n")
                f.write("-" * 25 + "\n")
                f.write(f"Genes with Novel Variants: {stats['variant_analysis']['genes_with_variants']}\n")
                f.write(f"Total Variants Detected: {stats['variant_analysis']['total_variants_detected']}\n")
                f.write(f"Average Variants per Gene: {stats['variant_analysis']['average_variants_per_gene']:.1f}\n")
                f.write(f"Variant Types Detected: {', '.join(stats['variant_analysis']['variant_types_detected'])}\n\n")
            
            # Immune Evasion Genes
            f.write("IMMUNE EVASION GENES ANALYSIS\n")
            f.write("-" * 35 + "\n")
            immune_stats = stats['gene_type_analysis']['immune_evasion_genes_detected']
            f.write(f"HLA-E Successfully Typed: {immune_stats['hla_e']} samples\n")
            f.write(f"HLA-F Successfully Typed: {immune_stats['hla_f']} samples\n")
            f.write(f"HLA-G Successfully Typed: {immune_stats['hla_g']} samples\n\n")
            
            # Recommendations for Oncogenetics
            f.write("RECOMMENDATIONS FOR ONCOGENETICS APPLICATIONS\n")
            f.write("-" * 50 + "\n")
            f.write("1. For comprehensive HLA typing: Use T1K or HLA-HD for broader gene coverage\n")
            f.write("2. For LOH detection: SpecHLA provides detailed allelic frequency analysis\n")
            f.write("3. For immune evasion analysis: Focus on HLA-E, -F, -G results from T1K/HLA-HD\n")
            f.write("4. For variant discovery: T1K provides novel variant detection capabilities\n")
            f.write("5. For clinical applications: Consider tool combination based on specific needs\n")
    
    def generate_enhanced_report(self, sample_list: List[str], output_dir: str = ".") -> None:
        """Generate comprehensive enhanced comparison report"""
        print("Starting Enhanced HLA Tools Comparison Analysis...")
        
        # Create comprehensive comparison table
        comparison_df = self.create_comprehensive_comparison_table(sample_list)
        
        # Generate enhanced statistics
        stats = self.create_enhanced_summary_statistics(comparison_df)
        
        # Create output directory
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)
        
        # Generate oncogenetics-focused report
        self.create_oncogenetics_focused_report(comparison_df, stats, output_dir)
        
        print(f"Enhanced analysis complete. Results saved to {output_path}")
        print("Generated files:")
        print("  - enhanced_hla_comparison.csv (Complete comparison)")
        print("  - classical_hla_comparison.csv (Classical genes only)")
        print("  - nonclassical_hla_comparison.csv (Non-classical genes only)")
        print("  - loh_analysis_summary.csv (LOH events analysis)")
        print("  - immune_evasion_genes.csv (HLA-E/F/G analysis)")
        print("  - variant_analysis_summary.csv (Novel variants)")
        print("  - enhanced_comparison_summary.txt (Detailed summary)")


def main():
    """Main execution function with enhanced functionality"""
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python enhanced_hla_comparison.py <base_directory> [sample1 sample2 ...]")
        print("Example: python enhanced_hla_comparison.py C:/hla_analysis 007_CZEPRS2xBRCA10386_run133")
        sys.exit(1)
    
    base_dir = sys.argv[1]
    
    # Get sample list from command line or auto-detect
    if len(sys.argv) > 2:
        sample_list = sys.argv[2:]
    else:
        # Auto-detect samples (same logic as original)
        extractor = EnhancedHLAExtractor(base_dir)
        sample_list = []
        
        # Check directories for samples
        if extractor.spechla_dir.exists():
            spechla_samples = [d.name for d in extractor.spechla_dir.iterdir() if d.is_dir()]
            sample_list.extend(spechla_samples)
        
        if extractor.t1k_dir.exists():
            t1k_files = [f.name for f in extractor.t1k_dir.iterdir() if f.name.endswith('_genotype.tsv')]
            t1k_samples = [f.replace('_t1k_result_genotype.tsv', '') for f in t1k_files]
            sample_list.extend(t1k_samples)
        
        if extractor.hlahd_dir.exists():
            hlahd_dirs = [d.name.replace('_hlahd', '') for d in extractor.hlahd_dir.iterdir() 
                         if d.is_dir() and d.name.endswith('_hlahd')]
            sample_list.extend(hlahd_dirs)
        
        sample_list = sorted(list(set(sample_list)))
    
    if not sample_list:
        print("No samples found. Please specify sample names manually.")
        sys.exit(1)
    
    print(f"Found samples: {sample_list}")
    
    # Create enhanced extractor and run analysis
    extractor = EnhancedHLAExtractor(base_dir)
    extractor.generate_enhanced_report(sample_list, output_dir=f"{base_dir}/enhanced_comparison_results")


if __name__ == "__main__":
    main()
