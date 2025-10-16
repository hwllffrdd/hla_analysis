#!/usr/bin/env python3
"""
Working MHCflurry Promiscuity Calculator
Fixed the array dimension issue - models are working perfectly!
"""

import pandas as pd
import numpy as np
import argparse
import sys
from typing import List, Dict, Optional

try:
    from mhcflurry import Class1AffinityPredictor
except ImportError:
    print("❌ MHCflurry not found. Make sure you're in the virtual environment:")
    print("source mhcflurry_env/bin/activate")
    sys.exit(1)

class MHCflurryPromiscuityCalculator:
    def __init__(self):
        """Initialize MHCflurry predictor"""
        print("🧬 Loading MHCflurry models...")
        try:
            self.predictor = Class1AffinityPredictor.load()
            print("✅ MHCflurry pan-allele models loaded successfully")
        except Exception as e:
            print(f"❌ Error loading MHCflurry: {e}")
            sys.exit(1)
    
    def load_hla_data(self, hla_file: str) -> pd.DataFrame:
        """Load HLA typing data"""
        try:
            df = pd.read_csv(hla_file, sep='\t')
            print(f"✅ Loaded HLA data: {len(df)} samples")
            return df
        except Exception as e:
            print(f"❌ Error loading HLA file: {e}")
            sys.exit(1)
    
    def format_alleles(self, hla_df: pd.DataFrame) -> List[str]:
        """Extract and format unique HLA alleles for MHCflurry"""
        alleles = set()
        hla_columns = ['A1', 'A2', 'B1', 'B2', 'C1', 'C2']
        
        for col in hla_columns:
            if col in hla_df.columns:
                col_alleles = hla_df[col].dropna().astype(str)
                for allele in col_alleles:
                    if allele and allele != 'nan':
                        # Convert A*24:02 -> HLA-A24:02
                        if allele.startswith(('A*', 'B*', 'C*')):
                            formatted = 'HLA-' + allele.replace('*', '')
                            alleles.add(formatted)
        
        result = sorted(list(alleles))
        print(f"✅ Found {len(result)} unique HLA alleles")
        return result
    
    def get_sample_alleles(self, hla_df: pd.DataFrame, sample_id: str) -> List[str]:
        """Get HLA alleles for a specific sample"""
        sample_row = hla_df[hla_df['Sample'] == sample_id]
        if sample_row.empty:
            return []
        
        alleles = []
        hla_columns = ['A1', 'A2', 'B1', 'B2', 'C1', 'C2']
        
        for col in hla_columns:
            if col in sample_row.columns:
                allele = sample_row[col].iloc[0]
                if pd.notna(allele) and str(allele) != '':
                    if str(allele).startswith(('A*', 'B*', 'C*')):
                        formatted = 'HLA-' + str(allele).replace('*', '')
                        alleles.append(formatted)
        
        return list(set(alleles))
    
    def predict_binding(self, peptides: List[str], alleles: List[str]) -> pd.DataFrame:
        """Run MHCflurry predictions"""
        if not peptides or not alleles:
            return pd.DataFrame()
        
        print(f"🔬 Running predictions: {len(peptides)} peptides × {len(alleles)} alleles = {len(peptides)*len(alleles)} total")
        
        try:
            # MHCflurry batch prediction - need to create all combinations manually
            all_peptides = []
            all_alleles = []
            
            # Create all peptide-allele combinations
            for peptide in peptides:
                for allele in alleles:
                    all_peptides.append(peptide)
                    all_alleles.append(allele)
            
            print(f"🔍 Created {len(all_peptides)} peptide-allele pairs")
            
            # Now both arrays have the same length
            predictions = self.predictor.predict(all_peptides, all_alleles)
            
            print(f"✅ Got {len(predictions)} predictions")
            
            # Convert to DataFrame
            results = []
            for i in range(len(all_peptides)):
                results.append({
                    'peptide': all_peptides[i],
                    'allele': all_alleles[i],
                    'binding_affinity_nM': predictions[i]
                })
            
            df = pd.DataFrame(results)
            
            # Calculate percentile ranks per allele (like NetMHCpan %Rank)
            for allele in df['allele'].unique():
                mask = df['allele'] == allele
                affinities = df.loc[mask, 'binding_affinity_nM']
                # Lower affinity = better binding = lower percentile
                percentiles = affinities.rank(pct=True) * 100
                df.loc[mask, 'percentile_rank'] = percentiles
            
            print(f"✅ Predictions complete: {len(df)} results")
            return df
            
        except Exception as e:
            print(f"❌ Batch prediction failed: {e}")
            
            # Fallback to individual predictions
            print("🔄 Trying individual predictions...")
            results = []
            
            for peptide in peptides:
                for allele in alleles:
                    try:
                        single_pred = self.predictor.predict([peptide], [allele])
                        results.append({
                            'peptide': peptide,
                            'allele': allele,
                            'binding_affinity_nM': single_pred[0]
                        })
                    except Exception as single_error:
                        print(f"⚠️  Failed prediction for {peptide} + {allele}: {single_error}")
                        continue
            
            if not results:
                print("❌ No predictions generated")
                return pd.DataFrame()
            
            df = pd.DataFrame(results)
            
            # Calculate percentile ranks
            for allele in df['allele'].unique():
                mask = df['allele'] == allele
                affinities = df.loc[mask, 'binding_affinity_nM']
                percentiles = affinities.rank(pct=True) * 100
                df.loc[mask, 'percentile_rank'] = percentiles
            
            print(f"✅ Individual predictions complete: {len(df)} results")
            return df
    
    def calculate_promiscuity(self, predictions_df: pd.DataFrame, 
                            threshold_percentile: float = 2.0,
                            threshold_affinity: Optional[float] = None) -> pd.DataFrame:
        """Calculate promiscuity scores"""
        if predictions_df.empty:
            return pd.DataFrame()
        
        print(f"📊 Calculating promiscuity scores...")
        
        promiscuity_results = []
        
        for peptide in predictions_df['peptide'].unique():
            peptide_data = predictions_df[predictions_df['peptide'] == peptide]
            total_alleles = len(peptide_data)
            
            # Determine binders
            if threshold_affinity is not None:
                binders = peptide_data[peptide_data['binding_affinity_nM'] <= threshold_affinity]
                threshold_desc = f"{threshold_affinity} nM"
            else:
                binders = peptide_data[peptide_data['percentile_rank'] <= threshold_percentile]
                threshold_desc = f"{threshold_percentile}% rank"
            
            num_binders = len(binders)
            promiscuity_score = num_binders / total_alleles if total_alleles > 0 else 0
            
            # Statistics
            median_affinity = peptide_data['binding_affinity_nM'].median()
            median_percentile = peptide_data['percentile_rank'].median()
            binding_alleles = binders['allele'].tolist()
            
            promiscuity_results.append({
                'Peptide': peptide,
                'Promiscuity_Score': promiscuity_score,
                'Binding_Percentage': promiscuity_score * 100,
                'Num_Binders': num_binders,
                'Total_Alleles': total_alleles,
                'Binding_Alleles': ','.join(binding_alleles),
                'Median_Affinity_nM': median_affinity,
                'Median_Percentile': median_percentile,
                'Threshold_Used': threshold_desc
            })
        
        df = pd.DataFrame(promiscuity_results)
        df = df.sort_values('Promiscuity_Score', ascending=False).reset_index(drop=True)
        
        print(f"✅ Promiscuity calculated for {len(df)} peptides")
        return df

def load_peptides(peptides_file: str) -> List[str]:
    """Load peptides from file"""
    try:
        with open(peptides_file, 'r') as f:
            peptides = [line.strip().upper() for line in f if line.strip()]
        
        # Filter valid peptides (8-15 amino acids, letters only)
        valid_peptides = []
        for peptide in peptides:
            if 8 <= len(peptide) <= 15 and peptide.isalpha():
                valid_peptides.append(peptide)
            else:
                print(f"⚠️  Skipping invalid peptide: {peptide}")
        
        print(f"✅ Loaded {len(valid_peptides)} valid peptides")
        return valid_peptides
        
    except FileNotFoundError:
        print(f"❌ Peptides file not found: {peptides_file}")
        sys.exit(1)

def print_summary(results_df: pd.DataFrame):
    """Print analysis summary"""
    if results_df.empty:
        return
    
    print(f"\n📈 Analysis Summary:")
    print(f"   Total peptides: {len(results_df)}")
    print(f"   Average promiscuity: {results_df['Promiscuity_Score'].mean():.3f}")
    print(f"   High promiscuity (>50%): {len(results_df[results_df['Promiscuity_Score'] > 0.5])}")
    print(f"   Very high promiscuity (>75%): {len(results_df[results_df['Promiscuity_Score'] > 0.75])}")
    
    print(f"\n🏆 Top 5 most promiscuous peptides:")
    top5 = results_df.head(5)
    for _, row in top5.iterrows():
        print(f"   {row['Peptide']:12} {row['Promiscuity_Score']:6.3f} "
              f"({row['Num_Binders']:2}/{row['Total_Alleles']:2} alleles) "
              f"{row['Median_Affinity_nM']:6.0f} nM")

def test_installation():
    """Quick test to verify everything works"""
    print("🧪 Testing MHCflurry installation...")
    
    calc = MHCflurryPromiscuityCalculator()
    
    # Test data
    test_peptides = ['SIINFEKL', 'SLYNTVATL']
    test_alleles = ['HLA-A02:01', 'HLA-B07:02']
    
    predictions_df = calc.predict_binding(test_peptides, test_alleles)
    
    if not predictions_df.empty:
        promiscuity_df = calc.calculate_promiscuity(predictions_df)
        print("\n🎉 Test successful! Sample results:")
        print(promiscuity_df[['Peptide', 'Promiscuity_Score', 'Num_Binders', 'Total_Alleles']])
        return True
    else:
        print("❌ Test failed")
        return False

def main():
    parser = argparse.ArgumentParser(
        description='MHCflurry Promiscuity Calculator',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Test installation
  python mhcflurry_promiscuity.py --test
  
  # Basic analysis  
  python mhcflurry_promiscuity.py hla_results.tsv peptides.txt
  
  # Population analysis with custom threshold
  python mhcflurry_promiscuity.py hla_results.tsv peptides.txt --threshold-percentile 1.0
  
  # Per-sample analysis
  python mhcflurry_promiscuity.py hla_results.tsv peptides.txt --per-sample
        """
    )
    
    parser.add_argument('hla_file', nargs='?', help='HLA typing results file')
    parser.add_argument('peptides_file', nargs='?', help='Peptides file (one per line)')
    parser.add_argument('-o', '--output', help='Output file')
    parser.add_argument('--threshold-percentile', type=float, default=2.0,
                       help='Percentile threshold for binding (default: 2.0)')
    parser.add_argument('--threshold-affinity', type=float,
                       help='Affinity threshold in nM (overrides percentile)')
    parser.add_argument('--per-sample', action='store_true',
                       help='Analyze per sample instead of population')
    parser.add_argument('--test', action='store_true',
                       help='Run installation test')
    
    args = parser.parse_args()
    
    # Test mode
    if args.test:
        if test_installation():
            print("\n✅ MHCflurry is working perfectly!")
            print("You can now run promiscuity analysis with your data.")
        sys.exit(0)
    
    # Regular analysis mode
    if not args.hla_file or not args.peptides_file:
        parser.error("HLA file and peptides file are required (or use --test)")
    
    # Load data
    peptides = load_peptides(args.peptides_file)
    if not peptides:
        print("❌ No valid peptides found")
        sys.exit(1)
    
    # Initialize calculator
    calc = MHCflurryPromiscuityCalculator()
    hla_df = calc.load_hla_data(args.hla_file)
    
    # Generate output filename if not provided
    if not args.output:
        base_name = args.hla_file.replace('.tsv', '').replace('.txt', '')
        threshold = args.threshold_affinity if args.threshold_affinity else args.threshold_percentile
        mode = "per_sample" if args.per_sample else "population"
        args.output = f"mhcflurry_promiscuity_{base_name}_{mode}_{threshold}.tsv"
    
    print(f"📁 Output file: {args.output}")
    
    if args.per_sample:
        # Per-sample analysis
        print("\n👥 Running per-sample analysis...")
        all_results = []
        
        for _, row in hla_df.iterrows():
            sample_id = row['Sample']
            print(f"\n🔬 Analyzing sample: {sample_id}")
            
            sample_alleles = calc.get_sample_alleles(hla_df, sample_id)
            if not sample_alleles:
                print(f"⚠️  No alleles found for {sample_id}")
                continue
            
            print(f"   HLA alleles: {', '.join(sample_alleles)}")
            
            predictions_df = calc.predict_binding(peptides, sample_alleles)
            if not predictions_df.empty:
                promiscuity_df = calc.calculate_promiscuity(
                    predictions_df, args.threshold_percentile, args.threshold_affinity
                )
                promiscuity_df['Sample'] = sample_id
                all_results.append(promiscuity_df)
        
        if all_results:
            combined_df = pd.concat(all_results, ignore_index=True)
            combined_df.to_csv(args.output, sep='\t', index=False)
            print(f"\n✅ Results saved to {args.output}")
            
            # Sample-wise summary
            print(f"\n👥 Per-sample summary:")
            for sample_id in combined_df['Sample'].unique():
                sample_data = combined_df[combined_df['Sample'] == sample_id]
                avg_prom = sample_data['Promiscuity_Score'].mean()
                high_prom = len(sample_data[sample_data['Promiscuity_Score'] > 0.5])
                print(f"   {sample_id}: avg={avg_prom:.3f}, high_prom={high_prom}")
        else:
            print("❌ No results generated")
    
    else:
        # Population analysis
        print("\n🌍 Running population analysis...")
        unique_alleles = calc.format_alleles(hla_df)
        
        predictions_df = calc.predict_binding(peptides, unique_alleles)
        if not predictions_df.empty:
            promiscuity_df = calc.calculate_promiscuity(
                predictions_df, args.threshold_percentile, args.threshold_affinity
            )
            
            promiscuity_df.to_csv(args.output, sep='\t', index=False)
            print(f"\n✅ Results saved to {args.output}")
            print_summary(promiscuity_df)
        else:
            print("❌ No results generated")

if __name__ == "__main__":
    main()
