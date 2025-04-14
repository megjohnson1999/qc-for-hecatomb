#!/usr/bin/env python3
"""
Script to verify and report read counts in paired-end FASTQ files.
This script checks if the number of reads in R1 and R2 files are equal.
"""

import os
import sys
import gzip
import glob
import argparse
from collections import defaultdict


def count_reads_in_fastq(fastq_file):
    """
    Count the number of reads in a gzipped FASTQ file.
    
    Args:
        fastq_file: Path to the FASTQ file (gzipped)
        
    Returns:
        Number of reads in the file
    """
    # Count lines and divide by 4 (each read has 4 lines in FASTQ format)
    with gzip.open(fastq_file, 'rt') as f:
        line_count = sum(1 for _ in f)
    
    return line_count // 4


def verify_paired_files(input_dir, file_pattern_r1, file_pattern_r2):
    """
    Verify that paired FASTQ files have the same number of reads.
    
    Args:
        input_dir: Directory containing FASTQ files
        file_pattern_r1: Glob pattern for R1 files
        file_pattern_r2: Glob pattern for R2 files
        
    Returns:
        Dictionary with sample names as keys and dictionaries containing read counts as values
    """
    results = defaultdict(dict)
    
    # Generate file paths for R1 files
    r1_files = glob.glob(os.path.join(input_dir, file_pattern_r1))
    
    for r1_file in r1_files:
        # Extract sample name from R1 file
        basename = os.path.basename(r1_file)
        sample_name = basename.split('_R1')[0]
        
        # Construct R2 file path
        r2_file = r1_file.replace('_R1', '_R2')
        if not os.path.exists(r2_file):
            print(f"WARNING: No matching R2 file found for {r1_file}")
            continue
        
        # Count reads in both files
        r1_count = count_reads_in_fastq(r1_file)
        r2_count = count_reads_in_fastq(r2_file)
        
        # Store results
        results[sample_name] = {
            'r1_file': r1_file,
            'r2_file': r2_file,
            'r1_reads': r1_count,
            'r2_reads': r2_count,
            'equal': r1_count == r2_count,
            'difference': abs(r1_count - r2_count)
        }
    
    return results


def main():
    parser = argparse.ArgumentParser(description='Verify paired-end FASTQ read counts.')
    parser.add_argument('--input-dir', required=True, help='Directory containing FASTQ files')
    parser.add_argument('--r1-pattern', default='*_R1*.fastq.gz', help='Glob pattern for R1 files')
    parser.add_argument('--r2-pattern', default='*_R2*.fastq.gz', help='Glob pattern for R2 files')
    parser.add_argument('--output', help='Output file (TSV format)')
    parser.add_argument('--stage', default='unknown', help='Processing stage name')
    
    args = parser.parse_args()
    
    print(f"Checking read pairs in: {args.input_dir}")
    print(f"Stage: {args.stage}")
    
    results = verify_paired_files(args.input_dir, args.r1_pattern, args.r2_pattern)
    
    # Print summary table
    print("\nResults:")
    print(f"{'Sample':<20} {'R1 Reads':<12} {'R2 Reads':<12} {'Equal':<8} {'Difference':<10}")
    print("-" * 65)
    
    all_equal = True
    for sample, data in sorted(results.items()):
        if not data['equal']:
            all_equal = False
        
        print(f"{sample:<20} {data['r1_reads']:<12} {data['r2_reads']:<12} {data['equal']:<8} {data['difference']:<10}")
    
    # Overall result
    print("\nSummary:")
    if all_equal:
        print("✅ All paired files have equal read counts")
    else:
        print("❌ Unequal read counts detected in paired files")
    
    # Write to output file if specified
    if args.output:
        with open(args.output, 'w') as f:
            f.write(f"sample\tr1_file\tr2_file\tr1_reads\tr2_reads\tequal\tdifference\tstage\n")
            for sample, data in sorted(results.items()):
                f.write(f"{sample}\t{data['r1_file']}\t{data['r2_file']}\t{data['r1_reads']}\t{data['r2_reads']}\t{data['equal']}\t{data['difference']}\t{args.stage}\n")
        print(f"Results written to {args.output}")


if __name__ == "__main__":
    main()