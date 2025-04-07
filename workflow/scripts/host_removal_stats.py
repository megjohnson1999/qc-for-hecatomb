import subprocess
import re
import os
import pandas as pd

def extract_host_removal_stats(merged_bam, unmerged_bam):
    sample = os.path.basename(merged_bam).replace('_merged.bam', '')
    data = {'sample': sample}

    # Get stats for merged reads
    merged_stats_cmd = f"samtools flagstat {merged_bam}"
    merged_stats = subprocess.check_output(merged_stats_cmd, shell=True).decode('utf-8')

    # Parse merged stats
    total_merged = int(re.search(r'(\d+) \+ \d+ in total', merged_stats).group(1))
    mapped_merged = int(re.search(r'(\d+) \+ \d+ mapped', merged_stats).group(1))
    unmapped_merged = total_merged - mapped_merged
    percent_host_in_merged = 100 * (mapped_merged / total_merged) if total_merged > 0 else 0

    # Get stats for unmerged paired reads
    unmerged_stats_cmd = f"samtools flagstat {unmerged_bam}"
    unmerged_stats = subprocess.check_output(unmerged_stats_cmd, shell=True).decode('utf-8')

    # Parse unmerged stats
    total_unmerged = int(re.search(r'(\d+) \+ \d+ in total', unmerged_stats).group(1))
    
    # For paired-end data, calculate unmapped pairs correctly
    unpaired_unmap_cmd = f"samtools view -f 4 -F 8 -c {unmerged_bam}"
    paired_unmap_cmd = f"samtools view -f 8 -F 4 -c {unmerged_bam}"
    unpaired_unmapped = int(subprocess.check_output(unpaired_unmap_cmd, shell=True).decode('utf-8').strip())
    paired_unmapped = int(subprocess.check_output(paired_unmap_cmd, shell=True).decode('utf-8').strip())
    unmapped_unmerged_pairs = (unpaired_unmapped + paired_unmapped // 2)
    
    mapped_unmerged_pairs = (total_unmerged // 2) - unmapped_unmerged_pairs
    percent_host_in_unmerged = 100 * (mapped_unmerged_pairs / (total_unmerged // 2)) if total_unmerged > 0 else 0

    # Calculate combined statistics
    total_reads = total_merged + total_unmerged
    mapped_reads = mapped_merged + (mapped_unmerged_pairs * 2)  # Each pair contributes two reads
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

# Example usage: this part mimics what Snakemake will provide
input_files_merged = ["sample1_merged.bam", "sample2_merged.bam"]  # Replace with actual merged BAMs from Snakemake
input_files_unmerged = ["sample1_unmerged.bam", "sample2_unmerged.bam"]  # Replace with actual unmerged BAMs from Snakemake
output_file = "host_removal_summary.tsv"

# Extract data from all BAM files
data_list = []
for merged_bam, unmerged_bam in zip(input_files_merged, input_files_unmerged):
    data = extract_host_removal_stats(merged_bam, unmerged_bam)
    data_list.append(data)

# Convert to DataFrame
df = pd.DataFrame(data_list)

# Save the DataFrame to a tabular file (TSV)
df.to_csv(output_file, sep='\t', index=False)
