import pandas as pd
import re
from glob import glob

input_files_step1 = snakemake.input.step1
input_files_step2 = snakemake.input.step2
output_file = snakemake.output[0]

# Function to extract read counts from BBDuk stats files
def extract_stats(stats_file):
    data = {'sample': stats_file.split('/')[-1].replace('_s1.stats', '').replace('_s2.stats', '')}
    total_reads = 0
    reads_with_hits = 0
    
    with open(stats_file) as f:
        for line in f:
            if "Total:" in line and total_reads == 0:
                total_reads = int(line.split()[1])
            if "Matched:" in line:
                reads_with_hits = int(line.split()[1])
    
    data['total_reads'] = total_reads
    data['reads_with_hits'] = reads_with_hits
    data['percent_with_hits'] = (reads_with_hits / total_reads * 100) if total_reads > 0 else 0
    
    return data

# Extract data from all stats files for step 1
data_step1 = [extract_stats(f) for f in input_files_step1]
df_step1 = pd.DataFrame(data_step1)
df_step1 = df_step1.rename(columns={
    'total_reads': 'total_reads_step1',
    'reads_with_hits': 'reads_with_primerB_hits',
    'percent_with_hits': 'percent_with_primerB'
})

# Extract data from all stats files for step 2
data_step2 = [extract_stats(f) for f in input_files_step2]
df_step2 = pd.DataFrame(data_step2)
df_step2 = df_step2.rename(columns={
    'total_reads': 'total_reads_step2',
    'reads_with_hits': 'reads_with_primerB_rc_hits',
    'percent_with_hits': 'percent_with_primerB_rc'
})

# Merge the two dataframes on sample
df = pd.merge(df_step1, df_step2, on='sample')

# Calculate combined metrics
df['total_primerB_hits'] = df['reads_with_primerB_hits'] + df['reads_with_primerB_rc_hits']
df['percent_total_primerB'] = (df['total_primerB_hits'] / df['total_reads_step1'] * 100)

# Save the DataFrame to a tabular file (TSV)
df.to_csv(output_file, sep='\t', index=False)