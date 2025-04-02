import pandas as pd
import re
import logging

# Setup logging
logging.basicConfig(filename="vector_stats_debug.log", level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

input_files = snakemake.input
output_file = snakemake.output[0]

def extract_vector_stats(stats_file):
    data = {'sample': stats_file.split('/')[-1].replace('_rm_vc.stats', '')}
    total_reads = 0
    reads_with_vector = 0
    vector_hits = {}

    logging.debug(f"Processing file: {stats_file}")

    with open(stats_file) as f:
        reading_vectors = False
        for line in f:
            line = line.strip()
            if line.startswith("#Total"):
                total_reads = int(line.split()[1])
            elif line.startswith("#Matched"):
                reads_with_vector = int(line.split()[1])
            elif line.startswith("#Name"):
                reading_vectors = True
                continue
            elif reading_vectors and line:
                logging.debug(f"Processing line: {line}")
                parts = re.split(r'\s{2,}|\t+', line)
                if len(parts) >= 4:
                    vector_name = " ".join(parts[:-2]).strip()
                    count = int(parts[-2])
                    logging.debug(f"Vector name: {vector_name}, Count: {count}")
                    vector_hits[vector_name] = count
                else:
                    logging.warning(f"Line skipped: {line}, Parts: {parts}")

    data['total_reads'] = total_reads
    data['reads_with_vector'] = reads_with_vector
    data['percent_with_vector'] = (reads_with_vector / total_reads * 100) if total_reads > 0 else 0

    if vector_hits:
        top_vector = max(vector_hits, key=vector_hits.get)
        data['top_vector'] = top_vector
    else:
        data['top_vector'] = 'None'

    logging.info(f"Sample: {data['sample']}, Total Reads: {total_reads}, Reads with Vector: {reads_with_vector}, Percent with Vector: {data['percent_with_vector']}, Top Vector: {data['top_vector']}")

    return data

# Extract data from all vector stats files
data_list = [extract_vector_stats(f) for f in input_files]

# Convert to DataFrame
df = pd.DataFrame(data_list)

# Save the DataFrame to a tabular file (TSV)
df.to_csv(output_file, sep='\t', index=False)

logging.info(f"Output written to {output_file}")
