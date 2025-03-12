import pandas as pd
import os
import subprocess
import re
from glob import glob

input_files_merged = snakemake.input.merged_bam
input_files_unmerged = snakemake.input.unmerged_bam
output_file = snakemake.output[0]

# Function to extract stats from BAM files
def extract_host_removal_stats(merged_bam, unmerged_bam):
    sample = os.path.basename(merged_bam).replace('_merged.bam', '')
    data = {'sample': sample}
    
    # Get stats for merged reads
    merged_stats_cmd = f"samtools flagstat {merged_bam}"
    merged_stats = subprocess.check_output(merged_stats_cmd, shell=True).decode('utf-8')
    
    # Parse merged stats
    total_merged = int(re.search(r'(\d+) \+ \d+ in total', merged_stats).group(1))
    unmapped_merged = int(re.search(r'(\d+) \+ \d+ mapped', merged_stats).group(1))
    
    # Get stats for unmerged paired reads
    unmerged_stats_cmd = f"samtools flagstat {unmerged_bam}"
    unmerged_stats = subprocess.check_output(unmerged_stats_cmd, shell=True).decode('utf-8')
    
    # Parse unmerged stats
    total_unmerged = int(re.search(r'(\d+) \+ \d+ in total', unmerged_stats).group(1))
    # For paired data, we need to look at read pairs that have both reads unmapped (flag 12)
    unmerged_unmapped_cmd = f"samtools view -c -f 12 {unmerged_bam}"
    unmapped_unmerged_pairs = int(subprocess.check_output(unmerged_unmapped_cmd, shell=True).decode('utf-8').strip())
    
    # Calculate statistics
    data['total_merged_reads'] = total_merged
    data['unmapped_merged_reads'] = unmapped_merged
    data['percent_host_in_merged'] = 100 * (1 - unmapped_merged / total_merged) if total_merged > 0 else 0
    
    data['total_unmerged_reads'] = total_unmerged
    data['unmapped_unmerged_pairs'] = unmapped_unmerged_pairs
    data['percent_host_in_unmerged'] = 100 * (1 - (unmapped_unmerged_pairs * 2) / total_unmerged) if total_unmerged > 0 else 0
    
    # Calculate combined statistics
    total_reads = total_merged + total_unmerged
    unmapped_reads = unmapped_merged + (unmapped_unmerged_pairs * 2)  # Count both reads in the pair
    data['total_reads'] = total_reads
    data['unmapped_reads'] = unmapped_reads
    data['percent_host_overall'] = 100 * (1 - unmapped_reads / total_reads) if total_reads > 0 else 0
    
    return data

# Extract data from all BAM files
data_list = []
for i in range(len(input_files_merged)):
    merged_bam = input_files_merged[i]
    unmerged_bam = input_files_unmerged[i]
    data = extract_host_removal_stats(merged_bam, unmerged_bam)
    data_list.append(data)

# Convert to DataFrame
df = pd.DataFrame(data_list)

# Save the DataFrame to a tabular file (TSV)
df.to_csv(output_file, sep='\t', index=False)