# CheckV vs Chimeric Detective Comparison Guide

This guide explains how to compare CheckV and Chimeric Detective results to get comprehensive assembly quality insights.

## Quick Start

```bash
# Run the comparison
python compare_checkv_ccd.py \
    --checkv-dir /path/to/checkv_output/ \
    --ccd-dir /path/to/chimeric_detective_output/ \
    --output-dir comparison_results/
```

## What the Tools Detect

| Tool | Primary Focus | Detected Issues |
|------|---------------|-----------------|
| **CheckV** | Viral genome completeness & contamination | Host gene contamination, incomplete genomes |
| **Chimeric Detective** | Assembly structural integrity | Viral-viral chimeras, PCR artifacts, misassemblies |

## Input Requirements

### CheckV Output Files
- `quality_summary.tsv` - Overall quality assessment
- `contamination.tsv` - Host contamination details  
- `complete_genomes.tsv` - Complete viral genomes (optional)
- `proviruses.tsv` - Integrated proviruses (optional)

### Chimeric Detective Output Files
- `splitting_decisions.tsv` - Chimera detection results
- `chimeric_detective_results.json` - Detailed analysis (optional)

## Output Files

The comparison generates:

- **`checkv_ccd_comparison.tsv`** - Detailed per-contig comparison
- **`comparison_summary.tsv`** - Summary statistics
- **`manual_inspection_targets.tsv`** - High-priority contigs for review
- **`checkv_ccd_comparison.png`** - Visualization plots

## Interpreting Results

### Agreement Categories

#### ✅ **High Confidence Issues** (Both tools agree)
- CheckV: "Low-quality" or "Contaminated"  
- CCD: "Split"
- **Action**: Safe to split these contigs

#### ✅ **Conservative Preservation** (Both tools agree)
- CheckV: "High-quality" or "Medium-quality"
- CCD: "Preserved"  
- **Action**: Keep these contigs intact

#### ⚠️ **Disagreements** (Manual review needed)

**Case 1: CheckV likes it, CCD splits it**
- CheckV: "High-quality"
- CCD: "Split"
- **Possible reasons**:
  - Biological recombination (real, should preserve)
  - Technical chimera with good viral content (should split)
  - Assembly spanning multiple viral strains (context-dependent)

**Case 2: CheckV flags it, CCD preserves it**
- CheckV: "Contaminated" or "Low-quality"
- CCD: "Preserved"
- **Possible reasons**:
  - Host contamination but structurally sound
  - Low completeness but not chimeric
  - Different quality thresholds

### Manual Inspection Priorities

The script identifies contigs needing manual review:

1. **High CheckV quality + CCD split** - Potential biological recombination
2. **CheckV contamination + CCD preserved** - Host genes but intact structure  
3. **High CCD confidence + CheckV issues** - Strong chimera evidence

## Example Analysis Workflow

```bash
# 1. Run both tools on your assembly
checkv end_to_end assembly.fasta checkv_output/
chimeric_detective --assembly assembly.fasta --reads-dir reads/ --out ccd_output/

# 2. Compare results  
python compare_checkv_ccd.py \
    --checkv-dir checkv_output/ \
    --ccd-dir ccd_output/ \
    --output-dir comparison/

# 3. Review summary
cat comparison/comparison_summary.tsv

# 4. Check high-priority targets
head -20 comparison/manual_inspection_targets.tsv

# 5. Make informed decisions
# - Split high-confidence problems
# - Preserve conservative agreements  
# - Manually review disagreements
```

## Decision Framework

### For Each Contig Category:

**Both tools flag issues** → **Split** (high confidence)

**Both tools approve** → **Preserve** (high confidence)  

**CheckV good + CCD split** → **Manual review**:
- Check for biological plausibility of recombination
- Examine breakpoint context
- Consider downstream analysis needs

**CheckV bad + CCD preserve** → **Manual review**:
- Assess impact of contamination vs structural integrity
- Consider filtering vs splitting options

## Common Patterns

### Viral Metagenomes
- **Many short contigs**: CheckV "Not-determined", CCD often skips (< min-length)
- **Host contamination**: CheckV detects, CCD may not (different focus)
- **Viral-viral chimeras**: CCD detects, CheckV may miss
- **Assembly artifacts**: Both tools often agree on problems

### Large Assemblies (>10K contigs)
- Focus on longer contigs (>2kb) where both tools provide most value
- Batch process disagreements by contig length/quality
- Prioritize contigs important for downstream analysis

## Tips for Your Analysis

1. **Start with agreement** - Handle high-confidence cases first
2. **Length matters** - Longer contigs deserve more attention  
3. **Context is key** - Consider your downstream analysis needs
4. **Validate manually** - Spot-check a few disagreements
5. **Document decisions** - Keep track of why you chose split vs preserve

## Common Questions

**Q: Why do the tools disagree?**
A: They detect different types of problems with different methods. Disagreement often indicates complex cases needing expert judgment.

**Q: Which tool should I trust more?**  
A: Neither - use both! CheckV excels at contamination detection, CCD at structural chimeras. Combined they give comprehensive quality assessment.

**Q: How many disagreements are normal?**
A: 10-30% disagreement is typical, depending on assembly quality and viral diversity.

**Q: Should I always split contigs both tools flag?**
A: Generally yes for high-confidence cases, but always consider biological context.