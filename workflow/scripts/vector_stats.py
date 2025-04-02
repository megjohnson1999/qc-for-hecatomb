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
            if line.startswith("#Total"):
                total_reads = int(line.split()[1])
            elif line.startswith("#Matched"):
                reads_with_vector = int(line.split()[1])
            elif line.startswith("#Name"):
                reading_vectors = True
                continue
            elif reading_vectors and line.strip() and not line.startswith('#'):
                parts = re.split(r'\s{2,}', line.strip())
                if len(parts) >= 2:
                    vector_name = parts[0]
                    count = int(parts[1])
                    vector_hits[vector_name] = count

    data['total_reads'] = total_reads
    data['reads_with_vector'] = reads_with_vector
    data['percent_with_vector'] = (reads_with_vector / total_reads * 100) if total_reads > 0 else 0

    # Determine the top vector
    if vector_hits:
        top_vector = max(vector_hits, key=vector_hits.get)
        data['top_vector'] = top_vector
    else:
        data['top_vector'] = 'None'

    return data

# Extract data from all vector stats files
data_list = [extract_vector_stats(f) for f in input_files]

# Convert to DataFrame
df = pd.DataFrame(data_list)

# Save the DataFrame to a tabular file (TSV)
df.to_csv(output_file, sep='\t', index=False)
