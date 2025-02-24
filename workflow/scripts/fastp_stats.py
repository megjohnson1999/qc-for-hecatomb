import json
import pandas as pd
from glob import glob

input_files = snakemake.input
output_file = snakemake.output[0]

# Specify the fields you want to extract
fields = [
    "summary.before_filtering.total_reads",
    "adapter_cutting.adapter_trimmed_reads"
]

# Function to extract the specific fields from a JSON file
def extract_fields(json_file):
    data = {}
    with open(json_file) as f:
        json_data = json.load(f)
        data['sample'] = json_file.split('/')[-1].replace('.json', '')
        for field in fields:
            keys = field.split('.')
            value = json_data
            for key in keys:
                value = value.get(key, None)
                if value is None:
                    break
            data[field] = value
    return data

# Extract data from all JSON files
data = [extract_fields(f) for f in input_files]

# Convert the data to a pandas DataFrame
df = pd.DataFrame(data)

# Save the DataFrame to a tabular file (TSV)
df.to_csv(output_file, sep='\t', index=False)
