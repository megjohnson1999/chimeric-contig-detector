#!/usr/bin/env python3
"""
Comprehensive benchmarking script for Chimeric Detective.

This script compares Chimeric Detective against other chimera detection tools
and assembly quality assessment tools to evaluate performance.
"""

import os
import sys
import json
import time
import shutil
import logging
import subprocess
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass
import tempfile
import argparse

@dataclass
class BenchmarkResult:
    """Results from a single tool run."""
    tool_name: str
    runtime_seconds: float
    chimeras_detected: int
    output_files: List[str]
    error_message: Optional[str] = None
    success: bool = True

@dataclass
class ToolConfig:
    """Configuration for a benchmarking tool."""
    name: str
    command_template: str
    check_command: Optional[str] = None
    required_files: List[str] = None
    parse_function: Optional[callable] = None

class ChimeraBenchmarker:
    """Main benchmarking class for chimera detection tools."""
    
    def __init__(self, assembly_file: str, reads_file1: str, reads_file2: str = None, 
                 bam_file: str = None, output_dir: str = "benchmark_results", 
                 threads: int = 4):
        """
        Initialize benchmarker.
        
        Args:
            assembly_file: Path to assembly FASTA file
            reads_file1: Path to forward reads (or single-end reads)
            reads_file2: Path to reverse reads (optional)
            bam_file: Path to BAM file (optional, will be created if not provided)
            output_dir: Output directory for results
            threads: Number of threads to use
        """
        self.assembly_file = Path(assembly_file)
        self.reads_file1 = Path(reads_file1)
        self.reads_file2 = Path(reads_file2) if reads_file2 else None
        self.bam_file = Path(bam_file) if bam_file else None
        self.output_dir = Path(output_dir)
        self.threads = threads
        
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Setup logging
        self.setup_logging()
        
        # Tool configurations
        self.tools = self._setup_tools()
        
        # Results storage
        self.results = []
        
    def setup_logging(self):
        """Setup logging configuration."""
        log_file = self.output_dir / "benchmark.log"
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler()
            ]
        )
        self.logger = logging.getLogger(__name__)
        
    def _setup_tools(self) -> Dict[str, ToolConfig]:
        """Setup tool configurations."""
        return {
            'chimeric_detective': ToolConfig(
                name='Chimeric Detective',
                command_template='chimeric_detective -a {assembly} -1 {reads1} -2 {reads2} -o {output} -t {threads}',
                check_command='chimeric_detective --version',
                required_files=['cleaned_assembly.fasta', 'results.json'],
                parse_function=self._parse_chimeric_detective_results
            ),
            'checkv': ToolConfig(
                name='CheckV',
                command_template='checkv end_to_end {assembly} {output} -t {threads}',
                check_command='checkv -h',
                required_files=['quality_summary.tsv', 'contamination.tsv'],
                parse_function=self._parse_checkv_results
            ),
            'virsorter2': ToolConfig(
                name='VirSorter2',
                command_template='virsorter run -w {output} -i {assembly} --min-length 1000 -j {threads}',
                check_command='virsorter -h',
                required_files=['final-viral-score.tsv', 'final-viral-combined.fa'],
                parse_function=self._parse_virsorter2_results
            ),
            'metaquast': ToolConfig(
                name='metaQUAST',
                command_template='metaquast.py {assembly} -o {output} --threads {threads} --max-ref-number 0 --min-contig 500',
                check_command='metaquast.py --version',
                required_files=['report.tsv', 'transposed_report.tsv'],
                parse_function=self._parse_metaquast_results
            ),
            'quast': ToolConfig(
                name='QUAST',
                command_template='quast.py {assembly} -o {output} --threads {threads} --min-contig 500',
                check_command='quast.py --version',
                required_files=['report.tsv'],
                parse_function=self._parse_quast_results
            )
        }
    
    def check_tool_availability(self) -> Dict[str, bool]:
        """Check which tools are available on the system."""
        available = {}
        
        for tool_key, tool_config in self.tools.items():
            if tool_config.check_command:
                try:
                    result = subprocess.run(
                        tool_config.check_command.split(),
                        capture_output=True,
                        text=True,
                        timeout=30
                    )
                    available[tool_key] = result.returncode == 0
                    if available[tool_key]:
                        self.logger.info(f"✅ {tool_config.name} is available")
                    else:
                        self.logger.warning(f"❌ {tool_config.name} is not available")
                except (subprocess.TimeoutExpired, FileNotFoundError):
                    available[tool_key] = False
                    self.logger.warning(f"❌ {tool_config.name} is not available")
            else:
                available[tool_key] = True  # Assume available if no check command
        
        return available
    
    def prepare_bam_file(self) -> Path:
        """Create BAM file if not provided."""
        if self.bam_file and self.bam_file.exists():
            return self.bam_file
        
        self.logger.info("Creating BAM file from reads...")
        bam_output = self.output_dir / "aligned_reads.bam"
        
        # Index assembly
        index_cmd = f"bwa index {self.assembly_file}"
        subprocess.run(index_cmd.split(), check=True)
        
        # Align reads
        if self.reads_file2:
            align_cmd = f"bwa mem -t {self.threads} {self.assembly_file} {self.reads_file1} {self.reads_file2}"
        else:
            align_cmd = f"bwa mem -t {self.threads} {self.assembly_file} {self.reads_file1}"
        
        # Convert to BAM and sort
        with open(bam_output, 'w') as bam_file:
            align_proc = subprocess.Popen(align_cmd.split(), stdout=subprocess.PIPE)
            sort_cmd = f"samtools sort -@ {self.threads} -o {bam_output} -"
            subprocess.run(sort_cmd.split(), stdin=align_proc.stdout, check=True)
            align_proc.wait()
        
        # Index BAM
        subprocess.run(f"samtools index {bam_output}".split(), check=True)
        
        self.bam_file = bam_output
        return bam_output
    
    def run_tool(self, tool_key: str) -> BenchmarkResult:
        """Run a single tool and return results."""
        tool_config = self.tools[tool_key]
        tool_output_dir = self.output_dir / tool_key
        tool_output_dir.mkdir(exist_ok=True)
        
        self.logger.info(f"Running {tool_config.name}...")
        
        # Prepare command
        if tool_key == 'chimeric_detective':
            if self.reads_file2:
                command = tool_config.command_template.format(
                    assembly=self.assembly_file,
                    reads1=self.reads_file1,
                    reads2=f"-2 {self.reads_file2}",
                    output=tool_output_dir,
                    threads=self.threads
                )
            else:
                command = tool_config.command_template.format(
                    assembly=self.assembly_file,
                    reads1=self.reads_file1,
                    reads2="",
                    output=tool_output_dir,
                    threads=self.threads
                )
        else:
            command = tool_config.command_template.format(
                assembly=self.assembly_file,
                output=tool_output_dir,
                threads=self.threads
            )
        
        # Run tool
        start_time = time.time()
        try:
            result = subprocess.run(
                command.split(),
                capture_output=True,
                text=True,
                timeout=3600  # 1 hour timeout
            )
            runtime = time.time() - start_time
            
            if result.returncode != 0:
                error_msg = f"Command failed with return code {result.returncode}: {result.stderr}"
                self.logger.error(f"{tool_config.name} failed: {error_msg}")
                return BenchmarkResult(
                    tool_name=tool_config.name,
                    runtime_seconds=runtime,
                    chimeras_detected=0,
                    output_files=[],
                    error_message=error_msg,
                    success=False
                )
            
            # Parse results
            chimeras_detected = 0
            output_files = []
            
            if tool_config.parse_function:
                chimeras_detected, output_files = tool_config.parse_function(tool_output_dir)
            
            self.logger.info(f"{tool_config.name} completed in {runtime:.2f}s, detected {chimeras_detected} chimeras")
            
            return BenchmarkResult(
                tool_name=tool_config.name,
                runtime_seconds=runtime,
                chimeras_detected=chimeras_detected,
                output_files=output_files,
                success=True
            )
            
        except subprocess.TimeoutExpired:
            runtime = time.time() - start_time
            error_msg = "Tool timed out after 1 hour"
            self.logger.error(f"{tool_config.name} timed out")
            return BenchmarkResult(
                tool_name=tool_config.name,
                runtime_seconds=runtime,
                chimeras_detected=0,
                output_files=[],
                error_message=error_msg,
                success=False
            )
        except Exception as e:
            runtime = time.time() - start_time
            error_msg = str(e)
            self.logger.error(f"{tool_config.name} failed with exception: {error_msg}")
            return BenchmarkResult(
                tool_name=tool_config.name,
                runtime_seconds=runtime,
                chimeras_detected=0,
                output_files=[],
                error_message=error_msg,
                success=False
            )
    
    def _parse_chimeric_detective_results(self, output_dir: Path) -> Tuple[int, List[str]]:
        """Parse Chimeric Detective results."""
        results_file = output_dir / "results.json"
        output_files = []
        
        if results_file.exists():
            with open(results_file) as f:
                data = json.load(f)
            chimeras_detected = len(data.get('chimeric_contigs', []))
            output_files = [str(f) for f in output_dir.iterdir() if f.is_file()]
        else:
            chimeras_detected = 0
        
        return chimeras_detected, output_files
    
    def _parse_virsorter2_results(self, output_dir: Path) -> Tuple[int, List[str]]:
        """Parse VirSorter2 results."""
        score_file = output_dir / "final-viral-score.tsv"
        output_files = []
        
        if score_file.exists():
            df = pd.read_csv(score_file, sep='\t')
            # Count sequences identified as viral (potential chimeras are those with mixed scores)
            # VirSorter2 doesn't directly detect chimeras, but low-confidence viral sequences
            # might indicate chimeric content
            chimeras_detected = len(df[df['max_score'] < 0.7])  # Low confidence threshold
            output_files = [str(f) for f in output_dir.iterdir() if f.is_file()]
        else:
            chimeras_detected = 0
        
        return chimeras_detected, output_files
    
    def _parse_checkv_results(self, output_dir: Path) -> Tuple[int, List[str]]:
        """Parse CheckV results."""
        contamination_file = output_dir / "contamination.tsv"
        output_files = []
        
        if contamination_file.exists():
            df = pd.read_csv(contamination_file, sep='\t')
            # Count contigs with contamination
            chimeras_detected = len(df[df['host_genes'] > 0])
            output_files = [str(f) for f in output_dir.iterdir() if f.is_file()]
        else:
            chimeras_detected = 0
        
        return chimeras_detected, output_files
    
    def _parse_metaquast_results(self, output_dir: Path) -> Tuple[int, List[str]]:
        """Parse metaQUAST results."""
        report_file = output_dir / "report.tsv"
        output_files = []
        
        if report_file.exists():
            df = pd.read_csv(report_file, sep='\t')
            # Look for misassemblies
            misassemblies = df.get('# misassemblies', [0])[0] if not df.empty else 0
            chimeras_detected = misassemblies
            output_files = [str(f) for f in output_dir.iterdir() if f.is_file()]
        else:
            chimeras_detected = 0
        
        return chimeras_detected, output_files
    
    def _parse_quast_results(self, output_dir: Path) -> Tuple[int, List[str]]:
        """Parse QUAST results."""
        report_file = output_dir / "report.tsv"
        output_files = []
        
        if report_file.exists():
            df = pd.read_csv(report_file, sep='\t')
            # Look for misassemblies
            misassemblies = df.get('# misassemblies', [0])[0] if not df.empty else 0
            chimeras_detected = misassemblies
            output_files = [str(f) for f in output_dir.iterdir() if f.is_file()]
        else:
            chimeras_detected = 0
        
        return chimeras_detected, output_files
    
    def run_benchmark(self, tools_to_run: List[str] = None) -> List[BenchmarkResult]:
        """Run benchmark comparison."""
        if tools_to_run is None:
            tools_to_run = list(self.tools.keys())
        
        # Check tool availability
        available_tools = self.check_tool_availability()
        tools_to_run = [t for t in tools_to_run if available_tools.get(t, False)]
        
        if not tools_to_run:
            self.logger.error("No tools available to run!")
            return []
        
        # Prepare BAM file if needed
        if 'chimeric_detective' in tools_to_run:
            self.prepare_bam_file()
        
        # Run each tool
        results = []
        for tool_key in tools_to_run:
            result = self.run_tool(tool_key)
            results.append(result)
            self.results.append(result)
        
        # Generate comparison report
        self.generate_report(results)
        
        return results
    
    def generate_report(self, results: List[BenchmarkResult]):
        """Generate benchmark comparison report."""
        # Create summary table
        summary_data = []
        for result in results:
            summary_data.append({
                'Tool': result.tool_name,
                'Success': result.success,
                'Runtime (seconds)': round(result.runtime_seconds, 2),
                'Chimeras Detected': result.chimeras_detected,
                'Error': result.error_message or 'None'
            })
        
        df = pd.DataFrame(summary_data)
        
        # Save summary
        summary_file = self.output_dir / "benchmark_summary.tsv"
        df.to_csv(summary_file, sep='\t', index=False)
        
        # Create detailed report
        report_file = self.output_dir / "benchmark_report.md"
        with open(report_file, 'w') as f:
            f.write("# Chimeric Detective Benchmark Report\n\n")
            f.write(f"**Assembly:** {self.assembly_file}\n")
            f.write(f"**Reads:** {self.reads_file1}")
            if self.reads_file2:
                f.write(f", {self.reads_file2}")
            f.write(f"\n**Threads:** {self.threads}\n")
            f.write(f"**Date:** {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            f.write("## Summary\n\n")
            f.write(df.to_markdown(index=False))
            f.write("\n\n")
            
            f.write("## Tool Details\n\n")
            for result in results:
                f.write(f"### {result.tool_name}\n\n")
                f.write(f"- **Success:** {result.success}\n")
                f.write(f"- **Runtime:** {result.runtime_seconds:.2f} seconds\n")
                f.write(f"- **Chimeras Detected:** {result.chimeras_detected}\n")
                if result.error_message:
                    f.write(f"- **Error:** {result.error_message}\n")
                f.write(f"- **Output Files:** {len(result.output_files)} files\n\n")
        
        self.logger.info(f"Benchmark report saved to {report_file}")
        self.logger.info(f"Summary table saved to {summary_file}")
        
        # Print summary to console
        print("\n" + "="*60)
        print("BENCHMARK RESULTS SUMMARY")
        print("="*60)
        print(df.to_string(index=False))
        print("="*60)

def create_synthetic_test_data(output_dir: Path, num_contigs: int = 10, 
                              chimera_fraction: float = 0.3) -> Tuple[Path, Path, Path]:
    """Create synthetic test data with known chimeras."""
    output_dir.mkdir(parents=True, exist_ok=True)
    
    assembly_file = output_dir / "synthetic_assembly.fasta"
    reads_file = output_dir / "synthetic_reads.fastq"
    truth_file = output_dir / "chimera_truth.json"
    
    import random
    import string
    
    # Generate random sequences
    def random_dna(length: int) -> str:
        return ''.join(random.choices('ATCG', k=length))
    
    chimeric_contigs = []
    contigs = []
    
    # Create normal contigs
    num_normal = int(num_contigs * (1 - chimera_fraction))
    for i in range(num_normal):
        seq = random_dna(random.randint(1000, 5000))
        contigs.append(f">contig_{i}\n{seq}")
    
    # Create chimeric contigs
    num_chimeric = num_contigs - num_normal
    for i in range(num_chimeric):
        # Create chimera by joining two sequences
        seq1 = random_dna(random.randint(500, 2000))
        seq2 = random_dna(random.randint(500, 2000))
        chimeric_seq = seq1 + seq2
        contig_id = f"contig_{num_normal + i}"
        contigs.append(f">{contig_id}\n{chimeric_seq}")
        
        chimeric_contigs.append({
            'contig_id': contig_id,
            'breakpoint': len(seq1),
            'total_length': len(chimeric_seq),
            'type': 'synthetic_chimera'
        })
    
    # Write assembly
    with open(assembly_file, 'w') as f:
        f.write('\n'.join(contigs))
    
    # Create simple reads (just substrings of contigs)
    with open(reads_file, 'w') as f:
        read_id = 0
        for contig in contigs:
            if contig.startswith('>'):
                continue
            seq = contig.strip()
            # Generate reads from this sequence
            for pos in range(0, len(seq) - 100, 200):
                read_seq = seq[pos:pos+100]
                f.write(f"@read_{read_id}\n{read_seq}\n+\n{'I'*len(read_seq)}\n")
                read_id += 1
    
    # Write truth data
    with open(truth_file, 'w') as f:
        json.dump({
            'chimeric_contigs': chimeric_contigs,
            'total_contigs': num_contigs,
            'chimera_fraction': chimera_fraction
        }, f, indent=2)
    
    return assembly_file, reads_file, truth_file

def main():
    """Main function."""
    parser = argparse.ArgumentParser(description="Benchmark Chimeric Detective against other tools")
    parser.add_argument('-a', '--assembly', required=True, help="Assembly FASTA file")
    parser.add_argument('-1', '--reads1', help="Forward reads file")
    parser.add_argument('-2', '--reads2', help="Reverse reads file")
    parser.add_argument('-b', '--bam', help="BAM file (optional)")
    parser.add_argument('-o', '--output', default="benchmark_results", help="Output directory")
    parser.add_argument('-t', '--threads', type=int, default=4, help="Number of threads")
    parser.add_argument('--tools', nargs='+', 
                       choices=['chimeric_detective', 'checkv', 'virsorter2', 'metaquast', 'quast'],
                       help="Tools to run (default: all available)")
    parser.add_argument('--create-test-data', action='store_true', 
                       help="Create synthetic test data")
    parser.add_argument('--test-contigs', type=int, default=20, 
                       help="Number of contigs in test data")
    parser.add_argument('--chimera-fraction', type=float, default=0.3,
                       help="Fraction of chimeric contigs in test data")
    
    args = parser.parse_args()
    
    if args.create_test_data:
        print("Creating synthetic test data...")
        test_dir = Path(args.output) / "test_data"
        assembly_file, reads_file, truth_file = create_synthetic_test_data(
            test_dir, args.test_contigs, args.chimera_fraction
        )
        print(f"Created test data in {test_dir}")
        print(f"Assembly: {assembly_file}")
        print(f"Reads: {reads_file}")
        print(f"Truth: {truth_file}")
        
        # Use created test data for benchmarking
        args.assembly = str(assembly_file)
        args.reads1 = str(reads_file)
    
    if not args.reads1 and not args.bam:
        print("Error: Must provide either reads files or BAM file")
        sys.exit(1)
    
    # Run benchmark
    benchmarker = ChimeraBenchmarker(
        assembly_file=args.assembly,
        reads_file1=args.reads1,
        reads_file2=args.reads2,
        bam_file=args.bam,
        output_dir=args.output,
        threads=args.threads
    )
    
    results = benchmarker.run_benchmark(args.tools)
    
    if not results:
        print("No tools could be run successfully!")
        sys.exit(1)
    
    print(f"\nBenchmark completed! Results saved to {args.output}")

if __name__ == "__main__":
    main()