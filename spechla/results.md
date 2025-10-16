# SpecHLA Results Analysis Guide

## Key Result Files (Start Here!)

### 1. **Primary HLA Typing Results**

#### `hla.result.txt` - **MAIN RESULT FILE** ⭐
This is your primary output with final HLA typing calls:
```bash
# View the main results
cat hla.result.txt
```
**Contains:** Final HLA allele calls for each locus (HLA-A, -B, -C, -DPA1, -DPB1, -DQA1, -DQB1, -DRB1)
**Format:** Each line shows locus and two alleles (maternal/paternal)

#### `hla.result.g.group.txt` - **G-GROUP RESOLUTION**
```bash
# View G-group results (standard clinical resolution)
cat hla.result.g.group.txt
```
**Contains:** HLA typing at G-group resolution (functionally equivalent alleles grouped)
**Use for:** Clinical interpretation and compatibility matching

#### `hla.result.details.txt` - **DETAILED TYPING INFO**
```bash
# View detailed typing information
cat hla.result.details.txt
```
**Contains:** Confidence scores, coverage depth, and quality metrics for each call

## 2. **Loss of Heterozygosity (LOH) Analysis**

SpecHLA's unique capability - check for immune evasion events:

### HLA Allele Sequence Files
#### `hla.allele.1.[LOCUS].fasta` and `hla.allele.2.[LOCUS].fasta`
```bash
# Check if both alleles are present (no LOH) or one is missing (LOH detected)
ls -la hla.allele.*.HLA_A.fasta
ls -la hla.allele.*.HLA_B.fasta
ls -la hla.allele.*.HLA_C.fasta
```

**Interpretation:**
- **Two files present** = Heterozygous (normal)
- **One file missing or significantly shorter** = Potential LOH
- **Check file sizes** - dramatic size differences indicate LOH

### Coverage Analysis
#### `007_CZEPRS2xBRCA10386_run133.realign.sort.bam.depth`
```bash
# Examine coverage depth across HLA regions
head -20 007_CZEPRS2xBRCA10386_run133.realign.sort.bam.depth
```
**Contains:** Per-position coverage depth - look for regions with ~50% expected coverage (indicates LOH)

## 3. **Variant and Phasing Information**

### VCF Files (Variants)
#### `HLA_[LOCUS].vcf.gz` - Raw variants
```bash
# View variants for specific locus
zcat HLA_A.vcf.gz | head -20
zcat HLA_B.vcf.gz | head -20
```

#### `HLA_[LOCUS].rephase.vcf.gz` - Phased variants
```bash
# View phased variants (shows which allele each variant belongs to)
zcat HLA_A.rephase.vcf.gz | grep -v "^#" | head -10
```

### Break Point Analysis
#### `HLA_[LOCUS]_break_points_phased.txt`
```bash
# View recombination breakpoints between alleles
cat HLA_A_break_points_phased.txt
cat HLA_B_break_points_phased.txt
```
**Contains:** Genomic coordinates where allele sequences switch

## 4. **Quality Control and Metrics**

### Frequency Information
#### `HLA_[LOCUS]_freq.txt`
```bash
# Check allele frequencies in your population
cat HLA_A_freq.txt
cat HLA_B_freq.txt
```
**Contains:** Population frequency data for identified alleles

### Low Coverage Regions
#### `low_depth.bed`
```bash
# Check regions with insufficient coverage
cat low_depth.bed
```
**Contains:** Genomic intervals with low sequencing coverage that might affect typing accuracy

## Practical Analysis Commands

### Quick Summary Check
```bash
# Get a quick overview of your results
echo "=== MAIN HLA TYPING RESULTS ==="
cat hla.result.txt

echo -e "\n=== G-GROUP RESULTS (Clinical) ==="
cat hla.result.g.group.txt

echo -e "\n=== DETAILED RESULTS ==="
cat hla.result.details.txt
```

### LOH Detection Analysis
```bash
# Check for potential LOH events
echo "=== LOH ANALYSIS ==="
for locus in A B C DPA1 DPB1 DQA1 DQB1 DRB1; do
    echo "HLA-$locus:"
    allele1_size=$(wc -c < hla.allele.1.HLA_${locus}.fasta 2>/dev/null || echo "0")
    allele2_size=$(wc -c < hla.allele.2.HLA_${locus}.fasta 2>/dev/null || echo "0")
    echo "  Allele 1 size: $allele1_size"
    echo "  Allele 2 size: $allele2_size"
    if [ $allele1_size -eq 0 ] || [ $allele2_size -eq 0 ]; then
        echo "  ⚠️  POTENTIAL LOH DETECTED"
    elif [ $((allele1_size * 2)) -lt $allele2_size ] || [ $((allele2_size * 2)) -lt $allele1_size ]; then
        echo "  ⚠️  SIGNIFICANT SIZE DIFFERENCE - CHECK MANUALLY"
    else
        echo "  ✅  Normal heterozygosity"
    fi
    echo
done
```

### Coverage Assessment
```bash
# Check overall coverage quality
echo "=== COVERAGE SUMMARY ==="
if [ -f "007_CZEPRS2xBRCA10386_run133.realign.sort.bam.depth" ]; then
    total_positions=$(wc -l < 007_CZEPRS2xBRCA10386_run133.realign.sort.bam.depth)
    avg_depth=$(awk '{sum+=$3} END {printf "%.1f", sum/NR}' 007_CZEPRS2xBRCA10386_run133.realign.sort.bam.depth)
    echo "Total HLA positions covered: $total_positions"
    echo "Average coverage depth: ${avg_depth}x"
fi

if [ -f "low_depth.bed" ]; then
    low_depth_regions=$(wc -l < low_depth.bed)
    echo "Low coverage regions: $low_depth_regions"
    if [ $low_depth_regions -gt 0 ]; then
        echo "⚠️  Some regions have low coverage - check low_depth.bed"
    fi
fi
```

## Key Interpretation Points

### 1. **Normal Results Should Show:**
- ✅ Two different alleles per locus in `hla.result.txt`
- ✅ Similar file sizes for `hla.allele.1.*` and `hla.allele.2.*`
- ✅ Even coverage across HLA regions
- ✅ Minimal entries in `low_depth.bed`

### 2. **LOH Warning Signs:**
- ⚠️ Missing allele files for any locus
- ⚠️ Dramatically different file sizes between allele 1 and 2
- ⚠️ Coverage drops to ~50% of expected in depth files
- ⚠️ Unusual break point patterns

### 3. **Quality Indicators:**
- **High confidence:** Deep coverage (>30x), clear allele calls, minimal low-depth regions
- **Moderate confidence:** Some low-depth regions but clear typing calls
- **Low confidence:** Extensive low-depth regions, ambiguous calls

## Files You Can Ignore (Unless Debugging)
- `*.log` files - processing logs
- `*_break_points_score.txt` - internal scoring metrics
- `*_break_points_spechap.txt` - internal algorithm data
- `select.DRB1.seq.txt` - internal sequence selection data

## Next Steps for Your Comparison Study

1. **Extract key results:** Focus on `hla.result.txt` and `hla.result.g.group.txt`
2. **Document LOH status:** Run the LOH detection script above
3. **Save summary:** Create a standardized output format for comparison with T1K and HLA-HD
4. **Note confidence levels:** Record any quality issues from `hla.result.details.txt`

This comprehensive output from SpecHLA gives you both standard HLA typing AND unique LOH detection capabilities that the other tools don't provide!
