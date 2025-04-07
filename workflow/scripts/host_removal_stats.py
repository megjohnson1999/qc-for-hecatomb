import subprocess
import os
import re
import pandas as pd

def extract_host_removal_stats(merged_fastq, unmerged_fastq_R1, unmerged_fastq_R2):
    sample = os.path.basename(merged_fastq).replace('_merged_hr.fastq.gz', '')
    data = {'sample': sample}

    # Get stats for merged reads using seqkit
    merged_stats_cmd = f"seqkit stats {merged_fastq}"
    merged_stats = subprocess.check_output(merged_stats_cmd, shell=True).decode('utf-8').splitlines()
    total_merged = int(merged_stats[1].split()[3])  # Total reads

    # Calculate unmapped reads assuming all reads are not removed by host filtering in this example.
    mapped_merged = 0  # Placeholder, as we don't have direct mapping info
    unmapped_merged = total_merged - mapped_merged
    percent_host_in_merged = 100 * (mapped_merged / total_merged) if total_merged > 0 else 0

    # Get stats for unmerged reads using seqkit
    unmerged_stats_cmd_R1 = f"seqkit stats {unmerged_fastq_R1}"
    unmerged_stats_cmd_R2 = f"seqkit stats {unmerged_fastq_R2}"
    unmerged_stats_R1 = subprocess.check_output(unmerged_stats_cmd_R1, shell=True).decode('utf-8').splitlines()
    unmerged_stats_R2 = subprocess.check_output(unmerged_stats_cmd_R2, shell=True).decode('utf-8').splitlines()
    
    total_unmerged_R1 = int(unmerged_stats_R1[1].split()[3])  # Total reads R1
    total_unmerged_R2 = int(unmerged_stats_R2[1].split()[3])  # Total reads R2
    total_unmerged = total_unmerged_R1 + total_unmerged_R2

    # Calculate unmapped pairs assuming all reads are not removed by host filtering in this example.
    mapped_unmerged_pairs = 0  # Placeholder, as we don't have direct mapping info
    unmapped_unmerged_pairs = total_unmerged - mapped_unmerged_pairs
    percent_host_in_unmerged = 100 * (mapped_unmerged_pairs / total_unmerged) if total_unmerged > 0 else 0

    # Calculate combined statistics
    total_reads = total_merged + total_unmerged
    mapped_reads = mapped_merged + mapped_unmerged_pairs
    percent_host_overall = 100 * (mapped_reads / total_reads) if total_reads > 0 else 0

    # Store values in the data dictionary
    data.update({
        'total_merged_reads': total_merged,
        'unmapped_merged_reads': unmapped_merged,
        'percent_host_in_merged': percent_host_in_merged,
        'total_unmerged_reads': total_unmerged,
        'unmapped_unmerged_pairs': unmapped_unmerged_pairs,
        'percent_host_in_unmerged': percent_host_in_unmerged,
        'total_reads': total_reads,
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
