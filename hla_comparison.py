#!/usr/bin/env python3
"""
HLA Tools Results Comparison Script
Extracts and compares results from SpecHLA, T1K, and HLA-HD
"""

import pandas as pd
import os
import re
from pathlib import Path
from typing import Dict, List, Tuple, Optional

class HLAResultsExtractor:
    def __init__(self, base_dir: str):
        self.base_dir = Path(base_dir)
        self.spechla_dir = self.base_dir / "results_spechla"  # Fixed directory name
        self.t1k_dir = self.base_dir / "results_t1k" 
        self.hlahd_dir = self.base_dir / "results_hlahd"
        
        # Classical HLA genes for comparison
        self.classical_genes = [
            'HLA-A', 'HLA-B', 'HLA-C',           # Class I
            'HLA-DRB1', 'HLA-DQB1', 'HLA-DPB1', # Class II major
            'HLA-DQA1', 'HLA-DPA1'              # Class II alpha
        ]
    
    def extract_spechla_results(self, sample_id: str) -> Dict:
        """Extract SpecHLA results for a given sample"""
        sample_dir = self.spechla_dir / sample_id
        
        results = {
            'sample_id': sample_id,
            'tool': 'SpecHLA',
            'resolution': '4-field',
            'hla_typing': {},
            'g_group': {},
            'quality_scores': {},
            'loh_status': {},
            'errors': []
        }
        
        try:
            # Main HLA typing results
            main_result_file = sample_dir / "hla.result.txt"
            if main_result_file.exists():
                with open(main_result_file, 'r') as f:
                    lines = f.readlines()
                    if len(lines) >= 3:  # Header comment + column headers + data
                        header = lines[1].strip().split('\t')  # Skip version line
                        data = lines[2].strip().split('\t')    # Data is on third line
                        
                        print(f"SpecHLA header: {header[:5]}...")  # Debug
                        print(f"SpecHLA data length: {len(data)}")  # Debug
                        
                        # Parse HLA alleles - skip first column which is sample name
                        for i, gene_field in enumerate(header[1:], 1):
                            if i < len(data):
                                gene_parts = gene_field.split('_')
                                if len(gene_parts) >= 3:
                                    gene = f"{gene_parts[0]}_{gene_parts[1]}"  # HLA_A
                                    allele_num = gene_parts[2]  # 1 or 2
                                    
                                    if gene not in results['hla_typing']:
                                        results['hla_typing'][gene] = {}
                                    results['hla_typing'][gene][f'allele_{allele_num}'] = data[i]
            
            # G-group results
            ggroup_file = sample_dir / "hla.result.g.group.txt"
            if ggroup_file.exists():
                with open(ggroup_file, 'r') as f:
                    lines = f.readlines()
                    if len(lines) >= 3:  # Header comment + column headers + data
                        header = lines[1].strip().split('\t')  # Skip version line
                        data = lines[2].strip().split('\t')    # Data is on third line
                        
                        for i, gene_field in enumerate(header[1:], 1):
                            if i < len(data):
                                gene_parts = gene_field.split('_')
                                if len(gene_parts) >= 3:
                                    gene = f"{gene_parts[0]}_{gene_parts[1]}"
                                    allele_num = gene_parts[2]
                                    
                                    if gene not in results['g_group']:
                                        results['g_group'][gene] = {}
                                    results['g_group'][gene][f'allele_{allele_num}'] = data[i]
            
            # Quality scores and details
            details_file = sample_dir / "hla.result.details.txt"
            if details_file.exists():
                with open(details_file, 'r') as f:
                    lines = f.readlines()
                    for line in lines[1:]:  # Skip header
                        if line.startswith('#'):
                            continue
                        parts = line.strip().split('\t')
                        if len(parts) >= 4:
                            gene = parts[0]
                            score_info = parts[3].split(';')
                            if len(score_info) >= 2:
                                results['quality_scores'][gene] = {
                                    'score': score_info[1],
                                    'details': parts[3]
                                }
            
            # LOH detection from frequency files
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
                            # Check for LOH (significant deviation from 0.5)
                            max_freq = max(frequencies)
                            min_freq = min(frequencies)
                            loh_detected = abs(max_freq - min_freq) > 0.3  # Threshold
                            
                            results['loh_status'][gene] = {
                                'frequencies': frequencies,
                                'het_variants': het_variants,
                                'loh_detected': loh_detected
                            }
                else:
                    # File doesn't exist - mark as no LOH data available
                    results['loh_status'][gene] = {
                        'frequencies': [],
                        'het_variants': 0,
                        'loh_detected': None  # None indicates no data, not no LOH
                    }
        
        except Exception as e:
            results['errors'].append(f"SpecHLA extraction error: {str(e)}")
        
        return results
    
    def extract_t1k_results(self, sample_id: str) -> Dict:
        """Extract T1K results for a given sample"""
        results = {
            'sample_id': sample_id,
            'tool': 'T1K',
            'resolution': '4-field',
            'hla_typing': {},
            'quality_scores': {},
            'read_counts': {},
            'alternatives': {},
            'extended_genes': {},
            'errors': []
        }
        
        try:
            # Main genotype results
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
                            
                            # Separate classical vs extended genes
                            gene_normalized = gene.replace('-', '_')  # Convert HLA-A to HLA_A for comparison
                            if gene_normalized in [g.replace('-', '_') for g in self.classical_genes]:
                                results['hla_typing'][gene] = gene_result
                            else:
                                results['extended_genes'][gene] = gene_result
            
            # Detailed allele list
            allele_file = self.t1k_dir / f"{sample_id}_t1k_result_allele.tsv"
            if allele_file.exists():
                with open(allele_file, 'r') as f:
                    for line in f:
                        parts = line.strip().split('\t')
                        if len(parts) >= 2:
                            allele = parts[0]
                            read_count = int(parts[1])
                            results['read_counts'][allele] = read_count
        
        except Exception as e:
            results['errors'].append(f"T1K extraction error: {str(e)}")
        
        return results
    
    def extract_hlahd_results(self, sample_id: str) -> Dict:
        """Extract HLA-HD results for a given sample"""
        results = {
            'sample_id': sample_id,
            'tool': 'HLA-HD',
            'resolution': '3-field',
            'hla_typing': {},
            'detailed_calls': {},
            'read_counts': {},
            'not_typed_genes': [],
            'all_genes': {},
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
                            
                            gene_result = {
                                'allele_1': allele_1,
                                'allele_2': allele_2,
                                'homozygous': (allele_2 == '-' and allele_1 is not None)
                            }
                            
                            if allele_1 is None:
                                results['not_typed_genes'].append(gene)
                            
                            # Store all genes, but separate classical for comparison
                            full_gene_name = f"HLA-{gene}" if not gene.startswith('HLA') else gene
                            results['all_genes'][full_gene_name] = gene_result
                            
                            if full_gene_name in self.classical_genes:
                                results['hla_typing'][full_gene_name] = gene_result
            
            # Extract detailed results for classical genes
            for gene in ['A', 'B', 'C', 'DRB1', 'DQA1', 'DQB1', 'DPA1', 'DPB1']:
                est_file = result_dir / f"{sample_id}_{gene}.est.txt"
                read_file = result_dir / f"{sample_id}_{gene}.read.txt"
                
                if est_file.exists():
                    with open(est_file, 'r') as f:
                        lines = f.readlines()
                        if len(lines) >= 3:  # Header lines + data
                            data_line = lines[2].strip().split('\t')
                            if len(data_line) >= 4:
                                results['detailed_calls'][f'HLA-{gene}'] = {
                                    'candidates_1': data_line[0].split(','),
                                    'candidates_2': data_line[1].split(','),
                                    'scores_1': data_line[2],
                                    'scores_2': data_line[3]
                                }
                
                if read_file.exists():
                    with open(read_file, 'r') as f:
                        lines = f.readlines()
                        gene_reads = {}
                        current_allele = None
                        
                        for line in lines:
                            if line.startswith('HLA-'):
                                parts = line.strip().split('\t')
                                if len(parts) >= 2:
                                    current_allele = parts[0]
                                    gene_reads[current_allele] = int(parts[1])
                        
                        if gene_reads:
                            results['read_counts'][f'HLA-{gene}'] = gene_reads
        
        except Exception as e:
            results['errors'].append(f"HLA-HD extraction error: {str(e)}")
        
        return results
    
    def normalize_allele_name(self, allele: str, target_resolution: str = "2-field") -> str:
        """Normalize allele names for comparison"""
        if not allele or allele == 'Not typed':
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
    
    def create_comparison_table(self, sample_list: List[str]) -> pd.DataFrame:
        """Create comprehensive comparison table for all samples and tools"""
        comparison_data = []
        
        for sample_id in sample_list:
            print(f"Processing sample: {sample_id}")
            
            # Extract results from all three tools
            print("  Extracting SpecHLA results...")
            spechla_results = self.extract_spechla_results(sample_id)
            print(f"    SpecHLA genes found: {list(spechla_results['hla_typing'].keys())}")
            
            print("  Extracting T1K results...")
            t1k_results = self.extract_t1k_results(sample_id)
            print(f"    T1K genes found: {list(t1k_results['hla_typing'].keys())}")
            
            print("  Extracting HLA-HD results...")
            hlahd_results = self.extract_hlahd_results(sample_id)
            print(f"    HLA-HD genes found: {list(hlahd_results['hla_typing'].keys())}")
            
            # Create comparison for each classical gene
            for gene in self.classical_genes:
                # Normalize gene names for each tool
                spechla_gene = gene.replace('-', '_')  # HLA-A -> HLA_A  
                t1k_gene = gene                        # HLA-A (keep as-is)
                hlahd_gene = gene                      # HLA-A (keep as-is)
                
                row = {
                    'Sample_ID': sample_id,
                    'Gene': gene,
                }
                
                print(f"    Processing gene {gene}:")
                print(f"      SpecHLA lookup: {spechla_gene}")
                print(f"      T1K lookup: {t1k_gene}")
                print(f"      HLA-HD lookup: {hlahd_gene}")
                
                # SpecHLA results
                if spechla_gene in spechla_results['hla_typing']:
                    spechla_data = spechla_results['hla_typing'][spechla_gene]
                    row['SpecHLA_Allele1'] = spechla_data.get('allele_1')
                    row['SpecHLA_Allele2'] = spechla_data.get('allele_2')
                    row['SpecHLA_Resolution'] = '4-field'
                    print(f"        SpecHLA found: {row['SpecHLA_Allele1']}, {row['SpecHLA_Allele2']}")
                    
                    # Add G-group if available
                    if spechla_gene in spechla_results['g_group']:
                        g_data = spechla_results['g_group'][spechla_gene]
                        row['SpecHLA_GGroup1'] = g_data.get('allele_1')
                        row['SpecHLA_GGroup2'] = g_data.get('allele_2')
                    
                    # LOH status
                    if spechla_gene in spechla_results['loh_status']:
                        loh_data = spechla_results['loh_status'][spechla_gene]
                        row['SpecHLA_LOH_Detected'] = loh_data['loh_detected']
                        row['SpecHLA_Het_Variants'] = loh_data['het_variants']
                else:
                    row.update({
                        'SpecHLA_Allele1': None, 'SpecHLA_Allele2': None,
                        'SpecHLA_Resolution': '4-field', 'SpecHLA_LOH_Detected': None
                    })
                    print(f"        SpecHLA not found")
                
                # T1K results  
                if t1k_gene in t1k_results['hla_typing']:
                    t1k_data = t1k_results['hla_typing'][t1k_gene]
                    row['T1K_Allele1'] = t1k_data.get('allele_1')
                    row['T1K_Allele2'] = t1k_data.get('allele_2')
                    row['T1K_Score1'] = t1k_data.get('score_1')
                    row['T1K_Score2'] = t1k_data.get('score_2')
                    row['T1K_Reads1'] = t1k_data.get('reads_1')
                    row['T1K_Reads2'] = t1k_data.get('reads_2')
                    row['T1K_Resolution'] = '4-field'
                    print(f"        T1K found: {row['T1K_Allele1']}, {row['T1K_Allele2']}")
                else:
                    row.update({
                        'T1K_Allele1': None, 'T1K_Allele2': None,
                        'T1K_Score1': None, 'T1K_Score2': None,
                        'T1K_Reads1': None, 'T1K_Reads2': None,
                        'T1K_Resolution': '4-field'
                    })
                    print(f"        T1K not found")
                
                # HLA-HD results  
                if hlahd_gene in hlahd_results['hla_typing']:
                    hlahd_data = hlahd_results['hla_typing'][hlahd_gene]
                    row['HLAHD_Allele1'] = hlahd_data.get('allele_1')
                    row['HLAHD_Allele2'] = hlahd_data.get('allele_2')
                    row['HLAHD_Homozygous'] = hlahd_data.get('homozygous', False)
                    row['HLAHD_Resolution'] = '3-field'
                    print(f"        HLA-HD found: {row['HLAHD_Allele1']}, {row['HLAHD_Allele2']}")
                else:
                    row.update({
                        'HLAHD_Allele1': None, 'HLAHD_Allele2': None,
                        'HLAHD_Homozygous': None, 'HLAHD_Resolution': '3-field'
                    })
                    print(f"        HLA-HD not found")
                
                # Add normalized 2-field comparison
                row['SpecHLA_2field_1'] = self.normalize_allele_name(row.get('SpecHLA_Allele1'))
                row['SpecHLA_2field_2'] = self.normalize_allele_name(row.get('SpecHLA_Allele2'))
                row['T1K_2field_1'] = self.normalize_allele_name(row.get('T1K_Allele1'))
                row['T1K_2field_2'] = self.normalize_allele_name(row.get('T1K_Allele2'))
                row['HLAHD_2field_1'] = self.normalize_allele_name(row.get('HLAHD_Allele1'))
                row['HLAHD_2field_2'] = self.normalize_allele_name(row.get('HLAHD_Allele2'))
                
                # Agreement analysis
                row['Agreement_All_Tools'] = self.check_tool_agreement(
                    [row['SpecHLA_2field_1'], row['SpecHLA_2field_2']],
                    [row['T1K_2field_1'], row['T1K_2field_2']],
                    [row['HLAHD_2field_1'], row['HLAHD_2field_2']]
                )
                
                row['Agreement_SpecHLA_T1K'] = self.check_pairwise_agreement(
                    [row['SpecHLA_2field_1'], row['SpecHLA_2field_2']],
                    [row['T1K_2field_1'], row['T1K_2field_2']]
                )
                
                row['Agreement_SpecHLA_HLAHD'] = self.check_pairwise_agreement(
                    [row['SpecHLA_2field_1'], row['SpecHLA_2field_2']],
                    [row['HLAHD_2field_1'], row['HLAHD_2field_2']]
                )
                
                row['Agreement_T1K_HLAHD'] = self.check_pairwise_agreement(
                    [row['T1K_2field_1'], row['T1K_2field_2']],
                    [row['HLAHD_2field_1'], row['HLAHD_2field_2']]
                )
                
                comparison_data.append(row)
        
        return pd.DataFrame(comparison_data)
    
    def check_tool_agreement(self, alleles1: List, alleles2: List, alleles3: List) -> str:
        """Check agreement between all three tools"""
        # Remove None values and sort for comparison
        set1 = set([a for a in alleles1 if a is not None])
        set2 = set([a for a in alleles2 if a is not None])
        set3 = set([a for a in alleles3 if a is not None])
        
        # Count how many tools have data
        tools_with_data = sum([bool(s) for s in [set1, set2, set3]])
        
        if tools_with_data >= 2:
            if set1 and set2 and set3:
                # All three tools have data
                if set1 == set2 == set3:
                    return "Full Agreement"
                elif (set1 == set2) or (set1 == set3) or (set2 == set3):
                    return "Partial Agreement (2/3 tools)"
                elif len(set1.intersection(set2).intersection(set3)) > 0:
                    return "Partial Agreement (some overlap)"
                else:
                    return "No Agreement"
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
    
    def check_pairwise_agreement(self, alleles1: List, alleles2: List) -> str:
        """Check agreement between two tools"""
        set1 = set([a for a in alleles1 if a is not None])
        set2 = set([a for a in alleles2 if a is not None])
        
        if set1 and set2:
            if set1 == set2:
                return "Full Agreement"
            elif len(set1.intersection(set2)) > 0:
                return "Partial Agreement"
            else:
                return "No Agreement"
        else:
            return "Insufficient Data"
    
    def create_summary_statistics(self, comparison_df: pd.DataFrame) -> Dict:
        """Generate summary statistics for the comparison"""
        stats = {
            'total_comparisons': len(comparison_df),
            'samples_analyzed': comparison_df['Sample_ID'].nunique(),
            'genes_analyzed': comparison_df['Gene'].nunique(),
        }
        
        # Agreement statistics - updated to match new categories
        all_agreement = comparison_df['Agreement_All_Tools']
        
        # Count full agreement (both 3 tools and 2 tools versions)
        stats['full_agreement_all_tools'] = (
            (all_agreement == 'Full Agreement') |
            (all_agreement == 'Full Agreement (2 tools)') |
            (all_agreement == 'Full Agreement (3 tools)')
        ).sum()
        
        # Count partial agreement (all variants)
        stats['partial_agreement_all_tools'] = (
            (all_agreement == 'Partial Agreement') |
            (all_agreement == 'Partial Agreement (2/3 tools)') |
            (all_agreement == 'Partial Agreement (some overlap)') |
            (all_agreement == 'Partial Agreement (2 tools)')
        ).sum()
        
        # Count no agreement
        stats['no_agreement_all_tools'] = (
            (all_agreement == 'No Agreement') |
            (all_agreement == 'No Agreement (2 tools)')
        ).sum()
        
        # Count insufficient data cases
        stats['insufficient_data'] = (
            (all_agreement == 'Insufficient Data') |
            (all_agreement == 'Single Tool Only')
        ).sum()
        
        # Pairwise agreement - count exact matches for "Full Agreement"
        stats['spechla_t1k_agreement'] = (comparison_df['Agreement_SpecHLA_T1K'] == 'Full Agreement').sum()
        stats['spechla_hlahd_agreement'] = (comparison_df['Agreement_SpecHLA_HLAHD'] == 'Full Agreement').sum()
        stats['t1k_hlahd_agreement'] = (comparison_df['Agreement_T1K_HLAHD'] == 'Full Agreement').sum()
        
        # Coverage statistics
        stats['spechla_coverage'] = (comparison_df['SpecHLA_Allele1'].notna()).sum()
        stats['t1k_coverage'] = (comparison_df['T1K_Allele1'].notna()).sum()
        stats['hlahd_coverage'] = (comparison_df['HLAHD_Allele1'].notna()).sum()
        
        # LOH detection (SpecHLA only)
        loh_detected = (comparison_df['SpecHLA_LOH_Detected'] == True).sum()
        stats['loh_events_detected'] = loh_detected
        
        # Detailed breakdown of agreement categories
        stats['agreement_breakdown'] = {}
        for category in all_agreement.unique():
            stats['agreement_breakdown'][category] = (all_agreement == category).sum()
        
        return stats
    
    def create_detailed_analysis(self, comparison_df: pd.DataFrame) -> Dict:
        """Create detailed analysis including discrepancies and tool-specific features"""
        analysis = {
            'discrepancies': [],
            'tool_specific_features': {},
            'resolution_differences': [],
            'quality_metrics': {}
        }
        
        # Find discrepancies
        disagreement_mask = comparison_df['Agreement_All_Tools'].isin(['No Agreement', 'Partial Agreement'])
        disagreements = comparison_df[disagreement_mask]
        
        for _, row in disagreements.iterrows():
            discrepancy = {
                'sample': row['Sample_ID'],
                'gene': row['Gene'],
                'spechla_calls': [row['SpecHLA_2field_1'], row['SpecHLA_2field_2']],
                't1k_calls': [row['T1K_2field_1'], row['T1K_2field_2']],
                'hlahd_calls': [row['HLAHD_2field_1'], row['HLAHD_2field_2']],
                'agreement_type': row['Agreement_All_Tools']
            }
            analysis['discrepancies'].append(discrepancy)
        
        # Tool-specific features
        analysis['tool_specific_features'] = {
            'SpecHLA': {
                'loh_detection': True,
                'g_group_assignment': True,
                'variant_phasing': True,
                'genes_covered': 8,
                'unique_features': ['LOH detection', 'G-group classification', 'Population frequencies']
            },
            'T1K': {
                'extended_gene_coverage': True,
                'read_count_metrics': True,
                'alternative_alleles': True,
                'genes_covered': '40+',
                'unique_features': ['Non-classical HLA genes', 'KIR typing capability', 'Read count support']
            },
            'HLA-HD': {
                'comprehensive_hla_coverage': True,
                'robust_quality_control': True,
                'homozygous_detection': True,
                'genes_covered': 29,
                'unique_features': ['Most comprehensive HLA coverage', 'Detailed quality metrics', 'Not typed reporting']
            }
        }
        
        return analysis
    
    def generate_report(self, sample_list: List[str], output_dir: str = ".") -> None:
        """Generate comprehensive comparison report"""
        print("Starting HLA Tools Comparison Analysis...")
        
        # Create comparison table
        comparison_df = self.create_comparison_table(sample_list)
        
        # Generate statistics
        stats = self.create_summary_statistics(comparison_df)
        detailed_analysis = self.create_detailed_analysis(comparison_df)
        
        # Save results
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)
        
        # Main comparison table
        comparison_df.to_csv(output_path / "hla_tools_comparison.csv", index=False)
        
        # Summary statistics
        with open(output_path / "comparison_summary.txt", 'w') as f:
            f.write("HLA Tools Comparison Summary\n")
            f.write("=" * 40 + "\n\n")
            
            f.write(f"Total Comparisons: {stats['total_comparisons']}\n")
            f.write(f"Samples Analyzed: {stats['samples_analyzed']}\n")
            f.write(f"Genes Analyzed: {stats['genes_analyzed']}\n\n")
            
            f.write("Agreement Statistics:\n")
            f.write(f"  Full Agreement (All combinations): {stats['full_agreement_all_tools']}\n")
            f.write(f"  Partial Agreement (All combinations): {stats['partial_agreement_all_tools']}\n")
            f.write(f"  No Agreement: {stats['no_agreement_all_tools']}\n")
            f.write(f"  Insufficient Data: {stats['insufficient_data']}\n\n")
            
            f.write("Detailed Agreement Breakdown:\n")
            for category, count in stats['agreement_breakdown'].items():
                f.write(f"  {category}: {count}\n")
            f.write("\n")
            
            f.write("Pairwise Agreement (Full Agreement Only):\n")
            f.write(f"  SpecHLA vs T1K: {stats['spechla_t1k_agreement']}\n")
            f.write(f"  SpecHLA vs HLA-HD: {stats['spechla_hlahd_agreement']}\n")
            f.write(f"  T1K vs HLA-HD: {stats['t1k_hlahd_agreement']}\n\n")
            
            f.write("Coverage Statistics:\n")
            f.write(f"  SpecHLA Coverage: {stats['spechla_coverage']}\n")
            f.write(f"  T1K Coverage: {stats['t1k_coverage']}\n")
            f.write(f"  HLA-HD Coverage: {stats['hlahd_coverage']}\n\n")
            
            f.write(f"LOH Events Detected (SpecHLA): {stats['loh_events_detected']}\n")
        
        # Discrepancies report
        if detailed_analysis['discrepancies']:
            discrepancy_df = pd.DataFrame(detailed_analysis['discrepancies'])
            discrepancy_df.to_csv(output_path / "discrepancies.csv", index=False)
        
        print(f"Analysis complete. Results saved to {output_path}")
        print(f"Main results: hla_tools_comparison.csv")
        print(f"Summary: comparison_summary.txt")
        if detailed_analysis['discrepancies']:
            print(f"Discrepancies: discrepancies.csv")


def main():
    """Main execution function"""
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python hla_comparison.py <base_directory> [sample1 sample2 ...]")
        print("Example: python hla_comparison.py C:/hla_analysis 007_CZEPRS2xBRCA10386_run133")
        sys.exit(1)
    
    base_dir = sys.argv[1]
    
    # Get sample list from command line or auto-detect
    if len(sys.argv) > 2:
        sample_list = sys.argv[2:]
    else:
        # Auto-detect samples from directory structure
        extractor = HLAResultsExtractor(base_dir)
        sample_list = []
        
        # Check SpecHLA directory for sample names
        if extractor.spechla_dir.exists():
            spechla_samples = [d.name for d in extractor.spechla_dir.iterdir() if d.is_dir()]
            print(f"SpecHLA samples found: {spechla_samples}")
            sample_list.extend(spechla_samples)
        
        # Check T1K directory for sample patterns
        if extractor.t1k_dir.exists():
            t1k_files = [f.name for f in extractor.t1k_dir.iterdir() if f.name.endswith('_genotype.tsv')]
            print(f"T1K files found: {t1k_files}")
            t1k_samples = [f.replace('_t1k_result_genotype.tsv', '') for f in t1k_files]
            print(f"T1K samples extracted: {t1k_samples}")
            sample_list.extend(t1k_samples)
        
        # Check HLA-HD directory
        if extractor.hlahd_dir.exists():
            hlahd_dirs = [d.name.replace('_hlahd', '') for d in extractor.hlahd_dir.iterdir() 
                         if d.is_dir() and d.name.endswith('_hlahd')]
            sample_list.extend(hlahd_dirs)
        
        # Remove duplicates and sort
        sample_list = sorted(list(set(sample_list)))
    
    if not sample_list:
        print("No samples found. Please specify sample names manually.")
        sys.exit(1)
    
    print(f"Found samples: {sample_list}")
    
    # Create extractor and run analysis
    extractor = HLAResultsExtractor(base_dir)
    extractor.generate_report(sample_list, output_dir=f"{base_dir}/comparison_results")


if __name__ == "__main__":
    main()
