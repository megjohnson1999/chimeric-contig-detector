#!/usr/bin/env python3
"""
Create a comprehensive test dataset for chimeric_detective.
Generates assembly with normal and chimeric contigs, plus synthetic reads.
"""

import os
import random
import gzip
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def generate_random_sequence(length, gc_content=0.5):
    """Generate a random DNA sequence with specified GC content."""
    gc_count = int(length * gc_content)
    at_count = length - gc_count
    
    bases = ['G', 'C'] * (gc_count // 2) + ['A', 'T'] * (at_count // 2)
    # Add remaining bases if odd numbers
    if gc_count % 2:
        bases.append('G')
    if at_count % 2:
        bases.append('A')
    
    random.shuffle(bases)
    return ''.join(bases)

def create_chimeric_contig(seq1, seq2, breakpoint_pos=None):
    """Create a chimeric contig by joining two sequences."""
    if breakpoint_pos is None:
        breakpoint_pos = len(seq1) // 2
    
    # Take first part of seq1 and second part of seq2
    chimeric_seq = seq1[:breakpoint_pos] + seq2[breakpoint_pos:]
    return chimeric_seq, breakpoint_pos

def generate_reads_from_sequence(sequence, read_length=150, coverage=20, error_rate=0.01):
    """Generate paired-end reads from a sequence."""
    seq_len = len(sequence)
    insert_size = 300
    read_pairs = []
    
    # Calculate number of reads needed for coverage
    num_reads = int((seq_len * coverage) / (read_length * 2))
    
    for _ in range(num_reads):
        # Random start position
        start_pos = random.randint(0, max(0, seq_len - insert_size))
        end_pos = min(start_pos + insert_size, seq_len)
        
        # Forward read
        r1_start = start_pos
        r1_end = min(r1_start + read_length, seq_len)
        read1 = sequence[r1_start:r1_end]
        
        # Reverse read (reverse complement of end region)
        r2_end = end_pos
        r2_start = max(r2_end - read_length, 0)
        read2_seq = sequence[r2_start:r2_end]
        read2 = str(Seq(read2_seq).reverse_complement())
        
        # Add sequencing errors
        if random.random() < error_rate:
            pos = random.randint(0, len(read1) - 1)
            bases = ['A', 'T', 'G', 'C']
            read1 = read1[:pos] + random.choice(bases) + read1[pos+1:]
        
        if random.random() < error_rate:
            pos = random.randint(0, len(read2) - 1)
            bases = ['A', 'T', 'G', 'C']
            read2 = read2[:pos] + random.choice(bases) + read2[pos+1:]
        
        read_pairs.append((read1, read2))
    
    return read_pairs

def write_fastq_gz(filename, reads, read_num):
    """Write reads to gzipped FASTQ file."""
    with gzip.open(filename, 'wt') as f:
        for i, read in enumerate(reads):
            f.write(f"@read_{i}/_{read_num}\n")
            f.write(f"{read}\n")
            f.write("+\n")
            f.write("I" * len(read) + "\n")  # High quality scores

def main():
    # Create test directory structure
    test_dir = Path("test_data_large")
    test_dir.mkdir(exist_ok=True)
    reads_dir = test_dir / "reads"
    reads_dir.mkdir(exist_ok=True)
    
    print("Creating large test dataset...")
    
    # Generate assembly with multiple contigs
    contigs = []
    chimeric_info = []
    
    # Normal contigs with varying GC content
    for i in range(20):
        gc_content = 0.3 + (i * 0.02)  # GC content from 30% to 68%
        length = random.randint(2000, 8000)
        seq = generate_random_sequence(length, gc_content)
        contigs.append(SeqRecord(Seq(seq), id=f"contig_{i+1:03d}", description=f"Normal contig {i+1}"))
    
    # Create some chimeric contigs
    chimeric_contigs = []
    for i in range(5):
        # Create two different sequences with different GC content
        seq1 = generate_random_sequence(3000, gc_content=0.3)
        seq2 = generate_random_sequence(3000, gc_content=0.7)
        
        chimeric_seq, breakpoint = create_chimeric_contig(seq1, seq2, 1500)
        contig_id = f"chimeric_{i+1:03d}"
        chimeric_contigs.append(SeqRecord(Seq(chimeric_seq), id=contig_id, 
                                        description=f"Chimeric contig {i+1}"))
        chimeric_info.append((contig_id, breakpoint, "GC_shift"))
    
    # Add coverage breakpoint chimeras
    for i in range(3):
        seq1 = generate_random_sequence(4000, gc_content=0.5)
        seq2 = generate_random_sequence(2000, gc_content=0.5)
        chimeric_seq, breakpoint = create_chimeric_contig(seq1, seq2, 2000)
        contig_id = f"cov_chimeric_{i+1:03d}"
        chimeric_contigs.append(SeqRecord(Seq(chimeric_seq), id=contig_id,
                                        description=f"Coverage chimeric contig {i+1}"))
        chimeric_info.append((contig_id, breakpoint, "coverage_shift"))
    
    all_contigs = contigs + chimeric_contigs
    
    # Write assembly file
    assembly_file = test_dir / "large_test_assembly.fasta"
    with open(assembly_file, 'w') as f:
        SeqIO.write(all_contigs, f, "fasta")
    
    print(f"Created assembly with {len(all_contigs)} contigs:")
    print(f"  - {len(contigs)} normal contigs")
    print(f"  - {len(chimeric_contigs)} chimeric contigs")
    
    # Generate reads for multiple samples
    sample_names = [f"sample_{i+1:03d}" for i in range(8)]
    
    for sample_name in sample_names:
        print(f"Generating reads for {sample_name}...")
        
        all_r1_reads = []
        all_r2_reads = []
        
        for contig in all_contigs:
            # Vary coverage for chimeric contigs to create coverage discontinuities
            if "chimeric" in contig.id:
                # Simulate different coverage for different parts
                coverage = random.randint(15, 40)
            else:
                coverage = random.randint(20, 50)
            
            read_pairs = generate_reads_from_sequence(str(contig.seq), coverage=coverage)
            
            for r1, r2 in read_pairs:
                all_r1_reads.append(r1)
                all_r2_reads.append(r2)
        
        # Write paired-end reads
        r1_file = reads_dir / f"{sample_name}_R1.fastq.gz"
        r2_file = reads_dir / f"{sample_name}_R2.fastq.gz"
        
        write_fastq_gz(r1_file, all_r1_reads, 1)
        write_fastq_gz(r2_file, all_r2_reads, 2)
    
    # Write chimeric contig info for validation
    info_file = test_dir / "chimeric_info.txt"
    with open(info_file, 'w') as f:
        f.write("contig_id\tbreakpoint\ttype\n")
        for contig_id, breakpoint, chimera_type in chimeric_info:
            f.write(f"{contig_id}\t{breakpoint}\t{chimera_type}\n")
    
    print(f"\nTest dataset created in {test_dir}/")
    print(f"Assembly: {assembly_file}")
    print(f"Reads: {reads_dir}/")
    print(f"Samples: {len(sample_names)}")
    print(f"Expected chimeric contigs: {len(chimeric_info)}")
    print(f"Validation info: {info_file}")

if __name__ == "__main__":
    main()