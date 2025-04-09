import subprocess
import os
import re
import pandas as pd
import os.path

def extract_host_removal_stats(merged_fastq, unmerged_fastq_R1, unmerged_fastq_R2):
    sample = os.path.basename(merged_fastq).replace('_merged_hr.fastq.gz', '')
    data = {'sample': sample}

    # Get stats for merged reads AFTER host removal using seqkit
    merged_stats_cmd = f"seqkit stats {merged_fastq}"
    merged_stats_output = subprocess.check_output(merged_stats_cmd, shell=True).decode('utf-8').splitlines()
    merged_header = merged_stats_output[0].split()
    merged_values = merged_stats_output[1].split()
    unmapped_merged = int(merged_values[merged_header.index('num_seqs')].replace(',', ''))  # Remove commas before converting

    # Get stats for unmerged reads AFTER host removal using seqkit
    unmerged_stats_cmd_R1 = f"seqkit stats {unmerged_fastq_R1}"
    unmerged_stats_output_R1 = subprocess.check_output(unmerged_stats_cmd_R1, shell=True).decode('utf-8').splitlines()
    unmerged_header_R1 = unmerged_stats_output_R1[0].split()
    unmerged_values_R1 = unmerged_stats_output_R1[1].split()
    unmapped_unmerged_R1 = int(unmerged_values_R1[unmerged_header_R1.index('num_seqs')].replace(',', ''))

    unmerged_stats_cmd_R2 = f"seqkit stats {unmerged_fastq_R2}"
    unmerged_stats_output_R2 = subprocess.check_output(unmerged_stats_cmd_R2, shell=True).decode('utf-8').splitlines()
    unmerged_header_R2 = unmerged_stats_output_R2[0].split()
    unmerged_values_R2 = unmerged_stats_output_R2[1].split()
    unmapped_unmerged_R2 = int(unmerged_values_R2[unmerged_header_R2.index('num_seqs')].replace(',', ''))

    # Get paths to the original files (before host removal)
    merged_input = merged_fastq.replace('host_removed', 'bbmerge').replace('_hr', '')
    unmerged_input_R1 = unmerged_fastq_R1.replace('host_removed', 'bbmerge').replace('_hr_R1', '_R1_unmerged')
    unmerged_input_R2 = unmerged_fastq_R2.replace('host_removed', 'bbmerge').replace('_hr_R2', '_R2_unmerged')

    # Get stats for merged reads BEFORE host removal
    if os.path.exists(merged_input):
        merged_in_stats_cmd = f"seqkit stats {merged_input}"
        merged_in_stats_output = subprocess.check_output(merged_in_stats_cmd, shell=True).decode('utf-8').splitlines()
        merged_in_header = merged_in_stats_output[0].split()
        merged_in_values = merged_in_stats_output[1].split()
        total_merged = int(merged_in_values[merged_in_header.index('num_seqs')].replace(',', ''))
    else:
        # If we can't find the original file, estimate
        total_merged = unmapped_merged  # Assume no host reads were found (will adjust if needed)

    # Get stats for unmerged reads BEFORE host removal
    if os.path.exists(unmerged_input_R1) and os.path.exists(unmerged_input_R2):
        unmerged_in_stats_cmd_R1 = f"seqkit stats {unmerged_input_R1}"
        unmerged_in_stats_output_R1 = subprocess.check_output(unmerged_in_stats_cmd_R1, shell=True).decode('utf-8').splitlines()
        unmerged_in_header_R1 = unmerged_in_stats_output_R1[0].split()
        unmerged_in_values_R1 = unmerged_in_stats_output_R1[1].split()
        total_unmerged_R1 = int(unmerged_in_values_R1[unmerged_in_header_R1.index('num_seqs')].replace(',', ''))

        unmerged_in_stats_cmd_R2 = f"seqkit stats {unmerged_input_R2}"
        unmerged_in_stats_output_R2 = subprocess.check_output(unmerged_in_stats_cmd_R2, shell=True).decode('utf-8').splitlines()
        unmerged_in_header_R2 = unmerged_in_stats_output_R2[0].split()
        unmerged_in_values_R2 = unmerged_in_stats_output_R2[1].split()
        total_unmerged_R2 = int(unmerged_in_values_R2[unmerged_in_header_R2.index('num_seqs')].replace(',', ''))
    else:
        # If we can't find the original files, estimate
        total_unmerged_R1 = unmapped_unmerged_R1  # Assume no host reads were found
        total_unmerged_R2 = unmapped_unmerged_R2  # Assume no host reads were found

    # Calculate host reads (mapped = total - unmapped)
    mapped_merged = total_merged - unmapped_merged
    
    # Make sure mapped reads is not negative (can happen due to estimation errors)
    mapped_merged = max(0, mapped_merged)
    
    # Calculate percentage of host reads in merged
    percent_host_in_merged = 100 * (mapped_merged / total_merged) if total_merged > 0 else 0

    # Calculate host reads for unmerged (ensuring pair counts match)
    unmapped_unmerged_pairs = min(unmapped_unmerged_R1, unmapped_unmerged_R2)  # Use the min to ensure pairs match
    total_unmerged_pairs = min(total_unmerged_R1, total_unmerged_R2)  # Use the min to ensure pairs match
    mapped_unmerged_pairs = total_unmerged_pairs - unmapped_unmerged_pairs
    
    # Make sure mapped pairs is not negative
    mapped_unmerged_pairs = max(0, mapped_unmerged_pairs)
    
    # Calculate percentage of host reads in unmerged
    percent_host_in_unmerged = 100 * (mapped_unmerged_pairs / total_unmerged_pairs) if total_unmerged_pairs > 0 else 0

    # Calculate combined statistics - count each read in a pair
    total_reads = total_merged + (total_unmerged_pairs * 2)
    unmapped_reads = unmapped_merged + (unmapped_unmerged_pairs * 2)
    mapped_reads = mapped_merged + (mapped_unmerged_pairs * 2)
    
    # Calculate overall host percentage
    percent_host_overall = 100 * (mapped_reads / total_reads) if total_reads > 0 else 0

    # Store values in the data dictionary
    data.update({
        'total_merged_reads': total_merged,
        'unmapped_merged_reads': unmapped_merged,
        'mapped_merged_reads': mapped_merged,
        'percent_host_in_merged': percent_host_in_merged,
        
        'total_unmerged_pairs': total_unmerged_pairs,
        'unmapped_unmerged_pairs': unmapped_unmerged_pairs,
        'mapped_unmerged_pairs': mapped_unmerged_pairs,
        'percent_host_in_unmerged': percent_host_in_unmerged,
        
        'total_reads': total_reads,
        'unmapped_reads': unmapped_reads,
        'mapped_reads': mapped_reads,
        'percent_host_overall': percent_host_overall
    })

    return data

# Snakemake variables
merged_fastq_files = snakemake.input.merged
unmerged_fastq_files_R1 = snakemake.input.unmerged_R1
unmerged_fastq_files_R2 = snakemake.input.unmerged_R2
output_file = snakemake.output[0]

# Extract data from all FASTQ files
data_list = []
for merged_fastq, unmerged_fastq_R1, unmerged_fastq_R2 in zip(merged_fastq_files, unmerged_fastq_files_R1, unmerged_fastq_files_R2):
    data = extract_host_removal_stats(merged_fastq, unmerged_fastq_R1, unmerged_fastq_R2)
    data_list.append(data)

# Convert to DataFrame
df = pd.DataFrame(data_list)

# Save the DataFrame to a tabular file (TSV)
df.to_csv(output_file, sep='\t', index=False)
