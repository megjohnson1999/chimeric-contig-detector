#!/usr/bin/env python3
"""
Example benchmarking script demonstrating how to compare chimera detection tools.
This script shows how to use the benchmarking suite with real or synthetic data.
"""

import sys
import subprocess
from pathlib import Path

def run_example_benchmark():
    """Run an example benchmark with synthetic data."""
    
    print("🧬 Chimeric Detective Benchmarking Example")
    print("=" * 50)
    
    # Check if main benchmark script exists
    benchmark_script = Path(__file__).parent / "benchmark_chimeric_detective.py"
    if not benchmark_script.exists():
        print("❌ Benchmark script not found!")
        return False
    
    print("📊 Creating synthetic test data and running benchmark...")
    print("")
    
    # Run benchmark with synthetic data
    cmd = [
        sys.executable, str(benchmark_script),
        "--create-test-data",
        "--test-contigs", "20",
        "--chimera-fraction", "0.25",
        "-o", "example_benchmark_output",
        "-t", "2"
    ]
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print("✅ Benchmark completed successfully!")
        print("")
        print("📁 Results saved to: example_benchmark_output/")
        print("")
        print("📋 Quick summary:")
        
        # Try to read summary file
        summary_file = Path("example_benchmark_output/benchmark_summary.tsv")
        if summary_file.exists():
            with open(summary_file) as f:
                lines = f.readlines()
                for line in lines[:6]:  # Show first few lines
                    print(f"   {line.strip()}")
        
        print("")
        print("💡 Next steps:")
        print("   1. Check the detailed report: example_benchmark_output/benchmark_report.md")
        print("   2. Review individual tool outputs in their respective directories")
        print("   3. Run with your own data using: python benchmark_chimeric_detective.py -a your_assembly.fasta -1 your_reads.fastq.gz")
        
        return True
        
    except subprocess.CalledProcessError as e:
        print(f"❌ Benchmark failed: {e}")
        print(f"Error output: {e.stderr}")
        return False
    except Exception as e:
        print(f"❌ Unexpected error: {e}")
        return False

def show_quick_benchmark_example():
    """Show example of using the quick benchmark script."""
    
    print("\n" + "=" * 50)
    print("🚀 Quick Benchmark Script Example")
    print("=" * 50)
    
    quick_script = Path(__file__).parent / "quick_benchmark.sh"
    if not quick_script.exists():
        print("❌ Quick benchmark script not found!")
        return
    
    print("The quick benchmark script can be used like this:")
    print("")
    print("# With your own data:")
    print(f"   {quick_script} -a your_assembly.fasta -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz")
    print("")
    print("# With single-end reads:")
    print(f"   {quick_script} -a your_assembly.fasta -1 reads.fastq.gz")
    print("")
    print("# With existing BAM file:")
    print(f"   {quick_script} -a your_assembly.fasta -b aligned_reads.bam")
    print("")
    print("# Custom settings:")
    print(f"   {quick_script} -a assembly.fasta -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -o my_results -t 8")
    print("")

def check_dependencies():
    """Check if required tools are available."""
    
    print("🔍 Checking tool availability...")
    print("-" * 30)
    
    tools = {
        'chimeric_detective': 'chimeric_detective --version',
        'checkv': 'checkv -h',
        'virsorter': 'virsorter -h',
        'quast': 'quast.py --version',
        'bwa': 'bwa',
        'samtools': 'samtools --version'
    }
    
    available = []
    missing = []
    
    for tool, cmd in tools.items():
        try:
            result = subprocess.run(cmd.split(), capture_output=True, timeout=10)
            if result.returncode == 0:
                available.append(tool)
                print(f"✅ {tool}")
            else:
                missing.append(tool)
                print(f"❌ {tool}")
        except (subprocess.TimeoutExpired, FileNotFoundError, subprocess.CalledProcessError):
            missing.append(tool)
            print(f"❌ {tool}")
    
    print("")
    print(f"Available tools: {len(available)}/{len(tools)}")
    
    if missing:
        print("")
        print("🔧 To install missing tools:")
        conda_tools = [t for t in missing if t != 'chimeric_detective']
        if conda_tools:
            # Special handling for virsorter2
            if 'virsorter' in conda_tools:
                conda_tools.remove('virsorter')
                conda_tools.append('virsorter=2')
            print("   conda install -c bioconda " + " ".join(conda_tools))
        print("")
        if 'chimeric_detective' in missing:
            print("   For Chimeric Detective:")
            print("   git clone https://github.com/megjohnson1999/chimeric-contig-detector.git")
            print("   cd chimeric-contig-detector")
            print("   conda env create -f environment.yml")
            print("   conda activate chimeric-detective")
            print("   pip install -e .")
    
    return len(missing) == 0

def main():
    """Main function."""
    
    print("🧬 Chimeric Detective Benchmarking Suite")
    print("=========================================")
    print("")
    print("This example demonstrates how to benchmark chimera detection tools.")
    print("")
    
    # Check dependencies
    all_available = check_dependencies()
    
    if not all_available:
        print("⚠️  Some tools are missing. The benchmark will run with available tools only.")
        print("")
    
    # Show quick benchmark example
    show_quick_benchmark_example()
    
    # Ask user if they want to run the example
    try:
        response = input("Do you want to run an example benchmark with synthetic data? (y/n): ")
        if response.lower().startswith('y'):
            success = run_example_benchmark()
            if success:
                print("\n🎉 Example completed successfully!")
            else:
                print("\n💥 Example failed. Check error messages above.")
        else:
            print("\n💡 You can run the benchmark manually using the commands shown above.")
    except KeyboardInterrupt:
        print("\n\n👋 Benchmark cancelled by user.")
    except Exception as e:
        print(f"\n❌ Error: {e}")

if __name__ == "__main__":
    main()