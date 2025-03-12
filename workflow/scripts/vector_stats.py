import pandas as pd
import re
from glob import glob

input_files = snakemake.input
output_file = snakemake.output[0]

# Function to extract read counts from BBDuk vector contamination stats files
def extract_vector_stats(stats_file):
    data = {'sample': stats_file.split('/')[-1].replace('_rm_vc.stats', '')}
    total_reads = 0
    reads_with_vector = 0
    vector_hits = {}
    
    with open(stats_file) as f:
        reading_vectors = False
        for line in f:
            if "Total:" in line and total_reads == 0:
                total_reads = int(line.split()[1])
            elif "Matched:" in line:
                reads_with_vector = int(line.split()[1])
            elif "#Ref" in line:
                reading_vectors = True
                continue
            elif reading_vectors and line.strip() and not line.startswith('#'):
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    vector_name = parts[0]
                    count = int(parts[1])
                    vector_hits[vector_name] = count
    
    data['total_reads'] = total_reads
    data['reads_with_vector'] = reads_with_vector
    data['percent_with_vector'] = (reads_with_vector / total_reads * 100) if total_reads > 0 else 0
    
    # Add the top 5 vectors by hit count
    top_vectors = sorted(vector_hits.items(), key=lambda x: x[1], reverse=True)[:5]
    for i, (vector, count) in enumerate(top_vectors, 1):
        data[f'top{i}_vector'] = vector
        data[f'top{i}_count'] = count
        data[f'top{i}_percent'] = (count / total_reads * 100) if total_reads > 0 else 0
    
    return data

# Extract data from all vector stats files
data_list = [extract_vector_stats(f) for f in input_files]

# Convert to DataFrame
df = pd.DataFrame(data_list)

# Save the DataFrame to a tabular file (TSV)
df.to_csv(output_file, sep='\t', index=False)