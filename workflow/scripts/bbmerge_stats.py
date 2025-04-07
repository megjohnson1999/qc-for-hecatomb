import pandas as pd
import re
from glob import glob
import numpy as np

input_files = snakemake.input
output_file = snakemake.output[0]
output_hist = snakemake.output[1]

# Function to extract merge stats from BBMerge output files
def extract_merge_stats(hist_file):
    data = {'sample': hist_file.split('/')[-1].replace('_bbmerge.out', '')}
    
    total_pairs = 0
    joined_pairs = 0
    avg_insert_size = 0
    insert_hist = {}
    median_insert_size = 0
    percent_merged = 0
    
    with open(hist_file) as f:
        for line in f:
            if line.startswith('#'):
                if "Mean" in line:
                    avg_insert_size = float(line.split()[1])
                elif "Median" in line:
                    median_insert_size = float(line.split()[1])
                elif "PercentOfPairs" in line:
                    percent_merged = float(line.split()[1])
            elif re.match(r'^\d+', line):  # This is a histogram line
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    insert_size = int(parts[0])
                    count = int(parts[1])
                    insert_hist[insert_size] = count
    
    data['total_pairs'] = total_pairs
    data['joined_pairs'] = joined_pairs
    data['percent_merged'] = percent_merged
    data['avg_insert_size'] = avg_insert_size
    data['median_insert_size'] = median_insert_size
    data['insert_histogram'] = insert_hist
    
    return data
    
# Extract data from all BBMerge output files
data_list = []
all_hist_data = {}

for f in input_files:
    sample_data = extract_merge_stats(f)
    data_list.append(sample_data)
    all_hist_data[sample_data['sample']] = sample_data['insert_histogram']
    
# Convert the main stats to a pandas DataFrame
df = pd.DataFrame(data_list)

# Drop the insert_histogram column since we'll save that separately
df = df.drop(columns=['insert_histogram'])

# Save the DataFrame to a tabular file (TSV)
df.to_csv(output_file, sep='\t', index=False)

# Create a DataFrame for the insert size histograms
# First, find the maximum insert size across all samples
max_insert = max([max(hist.keys()) for hist in all_hist_data.values() if hist])
insert_sizes = range(0, max_insert + 1)

# Create a DataFrame with insert sizes as rows and samples as columns
hist_df = pd.DataFrame(index=insert_sizes)

for sample, hist in all_hist_data.items():
    hist_series = pd.Series(hist)
    hist_series = hist_series.reindex(insert_sizes, fill_value=0)
    hist_df[sample] = hist_series

# Save the histogram DataFrame
hist_df.to_csv(output_hist, sep='\t')
