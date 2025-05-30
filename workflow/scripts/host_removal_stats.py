import subprocess
import os
import re
import pandas as pd
import os.path

def extract_host_removal_stats(r1_fastq, r2_fastq):
    sample = os.path.basename(r1_fastq).replace('_hr_R1.fastq', '')
    data = {'sample': sample}

    # Get stats for reads AFTER host removal using seqkit
    r1_stats_cmd = f"seqkit stats {r1_fastq}"
    r1_stats_output = subprocess.check_output(r1_stats_cmd, shell=True).decode('utf-8').splitlines()
    r1_header = r1_stats_output[0].split()
    r1_values = r1_stats_output[1].split()
    unmapped_r1 = int(r1_values[r1_header.index('num_seqs')].replace(',', ''))  # Remove commas before converting

    r2_stats_cmd = f"seqkit stats {r2_fastq}"
    r2_stats_output = subprocess.check_output(r2_stats_cmd, shell=True).decode('utf-8').splitlines()
    r2_header = r2_stats_output[0].split()
    r2_values = r2_stats_output[1].split()
    unmapped_r2 = int(r2_values[r2_header.index('num_seqs')].replace(',', ''))

    # Get paths to the original files (before host removal)
    # The input files are uncompressed fastq, but the original files are compressed
    original_r1 = r1_fastq.replace('host_removed', 'qc/rm_vector_contamination').replace('_hr_R1.fastq', '_R1_rm_vc.fastq.gz')
    original_r2 = r2_fastq.replace('host_removed', 'qc/rm_vector_contamination').replace('_hr_R2.fastq', '_R2_rm_vc.fastq.gz')

    # Get stats for reads BEFORE host removal
    if os.path.exists(original_r1) and os.path.exists(original_r2):
        original_r1_stats_cmd = f"seqkit stats {original_r1}"
        original_r1_stats_output = subprocess.check_output(original_r1_stats_cmd, shell=True).decode('utf-8').splitlines()
        original_r1_header = original_r1_stats_output[0].split()
        original_r1_values = original_r1_stats_output[1].split()
        total_r1 = int(original_r1_values[original_r1_header.index('num_seqs')].replace(',', ''))

        original_r2_stats_cmd = f"seqkit stats {original_r2}"
        original_r2_stats_output = subprocess.check_output(original_r2_stats_cmd, shell=True).decode('utf-8').splitlines()
        original_r2_header = original_r2_stats_output[0].split()
        original_r2_values = original_r2_stats_output[1].split()
        total_r2 = int(original_r2_values[original_r2_header.index('num_seqs')].replace(',', ''))
    else:
        # If we can't find the original files, estimate
        total_r1 = unmapped_r1  # Assume no host reads were found
        total_r2 = unmapped_r2  # Assume no host reads were found

    # Calculate host reads (mapped = total - unmapped)
    mapped_r1 = total_r1 - unmapped_r1
    mapped_r2 = total_r2 - unmapped_r2
    
    # Make sure mapped reads is not negative (can happen due to estimation errors)
    mapped_r1 = max(0, mapped_r1)
    mapped_r2 = max(0, mapped_r2)
    
    # Calculate percentage of host reads
    percent_host_r1 = 100 * (mapped_r1 / total_r1) if total_r1 > 0 else 0
    percent_host_r2 = 100 * (mapped_r2 / total_r2) if total_r2 > 0 else 0

    # Ensure pair counts match
    total_pairs = min(total_r1, total_r2)  # Use the min to ensure pairs match
    unmapped_pairs = min(unmapped_r1, unmapped_r2)  # Use the min to ensure pairs match
    mapped_pairs = total_pairs - unmapped_pairs
    
    # Make sure mapped pairs is not negative
    mapped_pairs = max(0, mapped_pairs)
    
    # Calculate percentage of host reads
    percent_host = 100 * (mapped_pairs / total_pairs) if total_pairs > 0 else 0

    # Store values in the data dictionary
    data.update({
        'total_read_pairs': total_pairs,
        'unmapped_read_pairs': unmapped_pairs,
        'mapped_read_pairs': mapped_pairs,
        'percent_host': percent_host,
        
        'total_reads': total_pairs * 2,
        'unmapped_reads': unmapped_pairs * 2,
        'mapped_reads': mapped_pairs * 2
    })

    return data

# Snakemake variables
r1_fastq_files = snakemake.input.r1
r2_fastq_files = snakemake.input.r2
output_file = snakemake.output[0]

# Extract data from all FASTQ files
data_list = []
for r1_fastq, r2_fastq in zip(r1_fastq_files, r2_fastq_files):
    data = extract_host_removal_stats(r1_fastq, r2_fastq)
    data_list.append(data)

# Convert to DataFrame
df = pd.DataFrame(data_list)

# Save the DataFrame to a tabular file (TSV)
df.to_csv(output_file, sep='\t', index=False)
