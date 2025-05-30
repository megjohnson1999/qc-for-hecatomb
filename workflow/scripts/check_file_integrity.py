#!/usr/bin/env python3
"""
Simple script to check gzip file integrity and paired-end read counts.
"""

import os
import sys
import gzip
import glob
import argparse
import subprocess

def check_gzip_integrity(file_path):
    """Check if a gzipped file is intact"""
    try:
        result = subprocess.run(['gzip', '-t', file_path], 
                               capture_output=True, 
                               text=True, 
                               check=False)
        
        return result.returncode == 0, result.stderr.strip() if result.returncode != 0 else "OK"
    except Exception as e:
        return False, str(e)

def count_fastq_reads(file_path):
    """Count reads in a FASTQ file"""
    try:
        # First check integrity
        is_intact, error = check_gzip_integrity(file_path)
        if not is_intact:
            return -1, f"File corrupted: {error}"
        
        # Count lines and divide by 4
        result = subprocess.run(f"zcat {file_path} | wc -l", 
                              shell=True, 
                              capture_output=True, 
                              text=True, 
                              check=True)
        
        line_count = int(result.stdout.strip())
        read_count = line_count // 4
        
        if line_count % 4 != 0:
            return read_count, f"Warning: Line count not divisible by 4 ({line_count})"
        
        return read_count, "OK"
    
    except Exception as e:
        return -1, str(e)

def main():
    parser = argparse.ArgumentParser(description='Check FASTQ file integrity and read counts')
    parser.add_argument('-d', '--directory', required=True, help='Directory containing FASTQ files')
    parser.add_argument('-p', '--pattern', default='*.fastq.gz', help='File pattern to match')
    parser.add_argument('-o', '--output', help='Output file path')
    
    args = parser.parse_args()
    
    files = sorted(glob.glob(os.path.join(args.directory, args.pattern)))
    
    if not files:
        print(f"No files matching '{args.pattern}' found in {args.directory}")
        sys.exit(1)
    
    print(f"Checking {len(files)} files in {args.directory}")
    
    results = []
    for file_path in files:
        file_size = os.path.getsize(file_path) / (1024 * 1024)  # Size in MB
        is_intact, error = check_gzip_integrity(file_path)
        
        if is_intact:
            read_count, msg = count_fastq_reads(file_path)
            status = "✅ OK" if "OK" in msg else f"⚠️ {msg}"
        else:
            read_count = -1
            status = f"❌ Corrupted: {error}"
        
        results.append({
            'file': os.path.basename(file_path),
            'size_mb': file_size,
            'intact': is_intact,
            'reads': read_count,
            'status': status
        })
        
        print(f"{os.path.basename(file_path):40} {file_size:6.1f} MB  {read_count:10}  {status}")
    
    # Check paired files
    print("\nChecking paired files:")
    r1_files = [f for f in files if "_R1" in f]
    
    for r1_file in r1_files:
        r2_file = r1_file.replace("_R1", "_R2")
        if r2_file in files:
            r1_count, _ = count_fastq_reads(r1_file)
            r2_count, _ = count_fastq_reads(r2_file)
            
            if r1_count > 0 and r2_count > 0:
                if r1_count == r2_count:
                    print(f"{os.path.basename(r1_file)} and {os.path.basename(r2_file)}: ✅ Paired counts match ({r1_count})")
                else:
                    print(f"{os.path.basename(r1_file)} and {os.path.basename(r2_file)}: ❌ Count mismatch! R1: {r1_count}, R2: {r2_count}")
            else:
                print(f"{os.path.basename(r1_file)} and {os.path.basename(r2_file)}: ❌ Cannot verify (corrupted files)")
    
    if args.output:
        with open(args.output, 'w') as f:
            f.write("file\tsize_mb\tintact\treads\tstatus\n")
            for result in results:
                f.write(f"{result['file']}\t{result['size_mb']:.1f}\t{result['intact']}\t{result['reads']}\t{result['status']}\n")
        print(f"\nResults written to {args.output}")

if __name__ == "__main__":
    main()