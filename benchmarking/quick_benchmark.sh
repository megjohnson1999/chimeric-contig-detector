#!/bin/bash
"""
Quick benchmark script for comparing chimera detection tools.
This script provides a simplified way to run common comparisons.
"""

set -e

# Default parameters
ASSEMBLY=""
READS1=""
READS2=""
BAM=""
OUTPUT_DIR="quick_benchmark_$(date +%Y%m%d_%H%M%S)"
THREADS=4

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -a|--assembly)
            ASSEMBLY="$2"
            shift 2
            ;;
        -1|--reads1)
            READS1="$2"
            shift 2
            ;;
        -2|--reads2)
            READS2="$2"
            shift 2
            ;;
        -b|--bam)
            BAM="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 -a assembly.fasta -1 reads_R1.fastq.gz [-2 reads_R2.fastq.gz] [-b alignment.bam] [-o output_dir] [-t threads]"
            echo ""
            echo "This script compares multiple chimera detection tools:"
            echo "  - Chimeric Detective (our tool)"
            echo "  - VSEARCH (UCHIME algorithm)"
            echo "  - CheckV (viral contamination detection)"
            echo "  - QUAST (assembly quality assessment)"
            echo ""
            echo "Required tools (install with conda/mamba):"
            echo "  conda install -c bioconda vsearch checkv quast bwa samtools"
            exit 0
            ;;
        *)
            echo "Unknown parameter: $1"
            exit 1
            ;;
    esac
done

# Check required inputs
if [[ -z "$ASSEMBLY" ]]; then
    echo "Error: Assembly file (-a) is required"
    exit 1
fi

if [[ -z "$READS1" && -z "$BAM" ]]; then
    echo "Error: Either reads files (-1/-2) or BAM file (-b) is required"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

echo "==================================================================="
echo "QUICK CHIMERA DETECTION BENCHMARK"
echo "==================================================================="
echo "Assembly: $ASSEMBLY"
echo "Reads1: $READS1"
echo "Reads2: $READS2"
echo "BAM: $BAM"
echo "Output: $OUTPUT_DIR"
echo "Threads: $THREADS"
echo "Date: $(date)"
echo ""

# Function to check if command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Function to run a tool and capture timing
run_tool() {
    local tool_name="$1"
    local command="$2"
    local output_dir="$3"
    
    echo "-------------------------------------------------------------------"
    echo "Running $tool_name..."
    echo "Command: $command"
    echo "-------------------------------------------------------------------"
    
    mkdir -p "$output_dir"
    
    local start_time=$(date +%s)
    
    if eval "$command" > "$output_dir/stdout.log" 2> "$output_dir/stderr.log"; then
        local end_time=$(date +%s)
        local runtime=$((end_time - start_time))
        echo "✅ $tool_name completed successfully in ${runtime}s"
        echo "$runtime" > "$output_dir/runtime.txt"
        echo "SUCCESS" > "$output_dir/status.txt"
    else
        local end_time=$(date +%s)
        local runtime=$((end_time - start_time))
        echo "❌ $tool_name failed after ${runtime}s"
        echo "$runtime" > "$output_dir/runtime.txt"
        echo "FAILED" > "$output_dir/status.txt"
        echo "Error details in $output_dir/stderr.log"
    fi
    echo ""
}

# Prepare BAM file if needed
if [[ -z "$BAM" && -n "$READS1" ]]; then
    echo "Creating BAM file from reads..."
    
    # Index assembly
    if ! command_exists bwa; then
        echo "❌ BWA not found. Install with: conda install -c bioconda bwa"
        exit 1
    fi
    
    if ! command_exists samtools; then
        echo "❌ samtools not found. Install with: conda install -c bioconda samtools"
        exit 1
    fi
    
    bwa index "$ASSEMBLY"
    
    if [[ -n "$READS2" ]]; then
        bwa mem -t "$THREADS" "$ASSEMBLY" "$READS1" "$READS2" | samtools sort -@ "$THREADS" -o aligned_reads.bam -
    else
        bwa mem -t "$THREADS" "$ASSEMBLY" "$READS1" | samtools sort -@ "$THREADS" -o aligned_reads.bam -
    fi
    
    samtools index aligned_reads.bam
    BAM="aligned_reads.bam"
    echo "✅ BAM file created: $BAM"
    echo ""
fi

# 1. Run Chimeric Detective
if command_exists chimeric_detective; then
    if [[ -n "$READS2" ]]; then
        run_tool "Chimeric Detective" \
                "chimeric_detective -a $ASSEMBLY -1 $READS1 -2 $READS2 -o chimeric_detective_results -t $THREADS" \
                "chimeric_detective_results"
    else
        run_tool "Chimeric Detective" \
                "chimeric_detective -a $ASSEMBLY -r $READS1 -o chimeric_detective_results -t $THREADS" \
                "chimeric_detective_results"
    fi
else
    echo "❌ Chimeric Detective not found. Install from: https://github.com/megjohnson1999/chimeric-contig-detector"
    echo ""
fi

# 2. Run VSEARCH (UCHIME)
if command_exists vsearch; then
    run_tool "VSEARCH (UCHIME)" \
            "vsearch --uchime_denovo $ASSEMBLY --chimeras vsearch_results/chimeras.fasta --nonchimeras vsearch_results/nonchimeras.fasta --uchimeout vsearch_results/uchime.out --threads $THREADS" \
            "vsearch_results"
else
    echo "❌ VSEARCH not found. Install with: conda install -c bioconda vsearch"
    echo ""
fi

# 3. Run CheckV
if command_exists checkv; then
    run_tool "CheckV" \
            "checkv end_to_end $ASSEMBLY checkv_results -t $THREADS" \
            "checkv_results"
else
    echo "❌ CheckV not found. Install with: conda install -c bioconda checkv"
    echo ""
fi

# 4. Run QUAST
if command_exists quast.py; then
    run_tool "QUAST" \
            "quast.py $ASSEMBLY -o quast_results --threads $THREADS" \
            "quast_results"
else
    echo "❌ QUAST not found. Install with: conda install -c bioconda quast"
    echo ""
fi

# Generate summary report
echo "==================================================================="
echo "GENERATING SUMMARY REPORT"
echo "==================================================================="

# Create summary table
{
    echo -e "Tool\tStatus\tRuntime(s)\tChimeras_Detected\tNotes"
    
    # Chimeric Detective
    if [[ -f "chimeric_detective_results/status.txt" ]]; then
        status=$(cat chimeric_detective_results/status.txt)
        runtime=$(cat chimeric_detective_results/runtime.txt)
        if [[ -f "chimeric_detective_results/results.json" ]]; then
            chimeras=$(python3 -c "import json; data=json.load(open('chimeric_detective_results/results.json')); print(len(data.get('chimeric_contigs', [])))" 2>/dev/null || echo "0")
        else
            chimeras="0"
        fi
        echo -e "Chimeric_Detective\t$status\t$runtime\t$chimeras\tComprehensive analysis with resolution"
    else
        echo -e "Chimeric_Detective\tNOT_RUN\t-\t-\tTool not available"
    fi
    
    # VSEARCH
    if [[ -f "vsearch_results/status.txt" ]]; then
        status=$(cat vsearch_results/status.txt)
        runtime=$(cat vsearch_results/runtime.txt)
        if [[ -f "vsearch_results/chimeras.fasta" ]]; then
            chimeras=$(grep -c "^>" vsearch_results/chimeras.fasta 2>/dev/null || echo "0")
        else
            chimeras="0"
        fi
        echo -e "VSEARCH_UCHIME\t$status\t$runtime\t$chimeras\tDe novo chimera detection"
    else
        echo -e "VSEARCH_UCHIME\tNOT_RUN\t-\t-\tTool not available"
    fi
    
    # CheckV
    if [[ -f "checkv_results/status.txt" ]]; then
        status=$(cat checkv_results/status.txt)
        runtime=$(cat checkv_results/runtime.txt)
        if [[ -f "checkv_results/contamination.tsv" ]]; then
            chimeras=$(python3 -c "import pandas as pd; df=pd.read_csv('checkv_results/contamination.tsv', sep='\t'); print(len(df[df['host_genes'] > 0]))" 2>/dev/null || echo "0")
        else
            chimeras="0"
        fi
        echo -e "CheckV\t$status\t$runtime\t$chimeras\tViral contamination detection"
    else
        echo -e "CheckV\tNOT_RUN\t-\t-\tTool not available"
    fi
    
    # QUAST
    if [[ -f "quast_results/status.txt" ]]; then
        status=$(cat quast_results/status.txt)
        runtime=$(cat quast_results/runtime.txt)
        if [[ -f "quast_results/report.tsv" ]]; then
            chimeras=$(python3 -c "import pandas as pd; df=pd.read_csv('quast_results/report.tsv', sep='\t'); print(df.get('# misassemblies', [0])[0] if not df.empty else 0)" 2>/dev/null || echo "0")
        else
            chimeras="0"
        fi
        echo -e "QUAST\t$status\t$runtime\t$chimeras\tAssembly quality metrics"
    else
        echo -e "QUAST\tNOT_RUN\t-\t-\tTool not available"
    fi
    
} > benchmark_summary.tsv

# Display summary
echo ""
echo "BENCHMARK SUMMARY:"
echo "=================="
column -t -s $'\t' benchmark_summary.tsv

# Create markdown report
{
    echo "# Quick Chimera Detection Benchmark Report"
    echo ""
    echo "**Date:** $(date)"
    echo "**Assembly:** $ASSEMBLY"
    echo "**Reads:** $READS1"
    if [[ -n "$READS2" ]]; then
        echo "**Reads2:** $READS2"
    fi
    echo "**Threads:** $THREADS"
    echo ""
    echo "## Results Summary"
    echo ""
    echo "| Tool | Status | Runtime(s) | Chimeras Detected | Notes |"
    echo "|------|--------|------------|-------------------|-------|"
    
    tail -n +2 benchmark_summary.tsv | while IFS=$'\t' read -r tool status runtime chimeras notes; do
        echo "| $tool | $status | $runtime | $chimeras | $notes |"
    done
    
    echo ""
    echo "## Tool Details"
    echo ""
    
    for tool_dir in */; do
        if [[ -d "$tool_dir" && -f "$tool_dir/status.txt" ]]; then
            tool_name=$(basename "$tool_dir")
            echo "### $tool_name"
            echo ""
            echo "**Status:** $(cat "$tool_dir/status.txt")"
            echo "**Runtime:** $(cat "$tool_dir/runtime.txt")s"
            echo ""
            if [[ -f "$tool_dir/stderr.log" && -s "$tool_dir/stderr.log" ]]; then
                echo "**Errors/Warnings:**"
                echo "```"
                head -n 20 "$tool_dir/stderr.log"
                echo "```"
                echo ""
            fi
        fi
    done
    
    echo "## Installation Instructions"
    echo ""
    echo "To install missing tools:"
    echo ""
    echo "```bash"
    echo "# Install conda/mamba if not available"
    echo "conda install -c conda-forge mamba"
    echo ""
    echo "# Install bioinformatics tools"
    echo "mamba install -c bioconda vsearch checkv quast bwa samtools"
    echo ""
    echo "# Install Chimeric Detective"
    echo "git clone https://github.com/megjohnson1999/chimeric-contig-detector.git"
    echo "cd chimeric-contig-detector"
    echo "conda env create -f environment.yml"
    echo "conda activate chimeric-detective"
    echo "pip install -e ."
    echo "```"
    
} > benchmark_report.md

echo ""
echo "==================================================================="
echo "BENCHMARK COMPLETE!"
echo "==================================================================="
echo "Results saved to: $(pwd)"
echo ""
echo "📊 Summary table: benchmark_summary.tsv"
echo "📋 Full report: benchmark_report.md"
echo ""
echo "💡 Interpretation tips:"
echo "  - Higher chimera counts don't always mean better detection"
echo "  - Compare with known contamination levels if available"
echo "  - Check individual tool outputs for detailed analysis"
echo "  - Consider false positive rates vs sensitivity trade-offs"
echo ""