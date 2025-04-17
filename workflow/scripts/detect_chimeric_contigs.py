#!/usr/bin/env python3
"""
Detect chimeric contigs by analyzing per-sample coverage patterns.
This script identifies contigs where different regions are predominantly covered by reads
from different samples, which suggests the contig may be artificially chimeric.
"""

import os
import sys
import pandas as pd
import numpy as np
from collections import defaultdict

# Get parameters from Snakemake
min_contig_length = snakemake.params.min_contig_length
min_coverage = snakemake.params.min_coverage
min_windows = snakemake.params.min_windows
sample_names = snakemake.params.sample_names
coverage_files = snakemake.input.coverages
output_chimeric = snakemake.output.chimeric
output_heatmap = snakemake.output.heatmap
log_file = snakemake.log[0]

# Set up logging
import logging
logging.basicConfig(
    filename=log_file,
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

def load_coverage_files(coverage_files, sample_names):
    """Load all coverage BED files into a dictionary of DataFrames."""
    coverage_data = {}
    
    for i, file_path in enumerate(coverage_files):
        sample = sample_names[i]
        logging.info(f"Loading coverage data for sample {sample} from {file_path}")
        
        try:
            # Load BED file with coverage data
            df = pd.read_csv(file_path, sep='\t', header=None, 
                            names=['contig', 'start', 'end', 'count', 'bases', 'length', 'coverage'])
            
            # Keep only required columns
            df = df[['contig', 'start', 'end', 'coverage']]
            
            # Add to dictionary
            coverage_data[sample] = df
            logging.info(f"Loaded {len(df)} coverage windows for sample {sample}")
        except Exception as e:
            logging.error(f"Error loading coverage data for sample {sample}: {e}")
            sys.exit(1)
    
    return coverage_data

def combine_coverage_data(coverage_data, min_contig_length):
    """Combine coverage data from all samples into a unified DataFrame for analysis."""
    # Get all unique contigs across all samples
    all_contigs = set()
    for sample, df in coverage_data.items():
        all_contigs.update(df['contig'].unique())
    
    # Get contig lengths for filtering
    contig_lengths = {}
    for sample, df in coverage_data.items():
        for contig in df['contig'].unique():
            contig_df = df[df['contig'] == contig]
            if contig in contig_lengths:
                contig_lengths[contig] = max(contig_lengths[contig], contig_df['end'].max())
            else:
                contig_lengths[contig] = contig_df['end'].max()
    
    # Filter contigs by length
    filtered_contigs = [contig for contig, length in contig_lengths.items() 
                        if length >= min_contig_length]
    
    logging.info(f"Found {len(filtered_contigs)} contigs >= {min_contig_length}bp for analysis")
    
    # Create a combined DataFrame
    combined_data = []
    
    for contig in filtered_contigs:
        windows = set()
        for sample, df in coverage_data.items():
            contig_df = df[df['contig'] == contig]
            for _, row in contig_df.iterrows():
                windows.add((row['start'], row['end']))
        
        windows = sorted(windows)
        
        for start, end in windows:
            row_data = {'contig': contig, 'start': start, 'end': end}
            
            for sample, df in coverage_data.items():
                try:
                    coverage = df[(df['contig'] == contig) & 
                                  (df['start'] == start) & 
                                  (df['end'] == end)]['coverage'].values
                    
                    row_data[sample] = coverage[0] if len(coverage) > 0 else 0
                except:
                    row_data[sample] = 0
            
            combined_data.append(row_data)
    
    if not combined_data:
        logging.error("No data found after combining coverage information")
        sys.exit(1)
    
    return pd.DataFrame(combined_data)

def detect_chimeric_contigs(df, min_coverage, min_windows, sample_names):
    """Detect potentially chimeric contigs by identifying regions where different
    samples contribute coverage to different parts of the same contig."""
    
    chimeric_contigs = []
    heatmap_data = []
    
    # Process each contig
    for contig in df['contig'].unique():
        contig_df = df[df['contig'] == contig].sort_values('start')
        
        if len(contig_df) < 5:  # Skip very short contigs
            continue
        
        # Create a matrix of coverage values for visualization
        matrix_rows = []
        for _, row in contig_df.iterrows():
            coverage_row = [row[sample] for sample in sample_names]
            matrix_rows.append(coverage_row)
        
        coverage_matrix = np.array(matrix_rows)
        
        # Determine which sample contributes the most to each window
        dominant_samples = []
        for i, row in enumerate(coverage_matrix):
            if np.max(row) >= min_coverage:
                dominant_sample = sample_names[np.argmax(row)]
                dominant_samples.append((i, dominant_sample, np.max(row)))
        
        # Check for changes in the dominant sample along the contig
        sample_regions = defaultdict(list)
        current_sample = None
        current_start = 0
        current_count = 0
        
        for i, (win_idx, sample, cov) in enumerate(dominant_samples):
            if current_sample is None:
                current_sample = sample
                current_start = win_idx
                current_count = 1
            elif sample == current_sample:
                current_count += 1
            else:
                if current_count >= min_windows:
                    sample_regions[current_sample].append((current_start, win_idx-1))
                current_sample = sample
                current_start = win_idx
                current_count = 1
        
        # Add the last region
        if current_sample is not None and current_count >= min_windows:
            sample_regions[current_sample].append((current_start, len(dominant_samples)-1))
        
        # Check if we have regions from different samples
        if len(sample_regions) > 1:
            region_info = []
            for sample, regions in sample_regions.items():
                for start, end in regions:
                    start_pos = contig_df.iloc[dominant_samples[start][0]]['start']
                    end_pos = contig_df.iloc[dominant_samples[end][0]]['end']
                    region_info.append(f"{sample}:{start_pos}-{end_pos}")
            
            if len(region_info) >= 2:
                chimeric_contigs.append((contig, "; ".join(region_info)))
                
                # Create heatmap data
                for i, row in contig_df.iterrows():
                    heatmap_row = [contig, row['start'], row['end']]
                    for sample in sample_names:
                        heatmap_row.append(row[sample])
                    heatmap_data.append(heatmap_row)
    
    logging.info(f"Detected {len(chimeric_contigs)} potentially chimeric contigs")
    return chimeric_contigs, heatmap_data

def main():
    """Main function to execute the analysis."""
    # Load coverage data from all samples
    coverage_data = load_coverage_files(coverage_files, sample_names)
    
    # Combine all coverage data
    combined_df = combine_coverage_data(coverage_data, min_contig_length)
    
    # Detect chimeric contigs
    chimeric_contigs, heatmap_data = detect_chimeric_contigs(
        combined_df, min_coverage, min_windows, sample_names
    )
    
    # Write results to output file
    with open(output_chimeric, 'w') as f:
        f.write("contig\tregions\n")
        for contig, regions in chimeric_contigs:
            f.write(f"{contig}\t{regions}\n")
    
    # Write heatmap data
    header = ['contig', 'start', 'end'] + sample_names
    heatmap_df = pd.DataFrame(heatmap_data, columns=header)
    heatmap_df.to_csv(output_heatmap, sep='\t', index=False)
    
    logging.info(f"Analysis complete. Results written to {output_chimeric} and {output_heatmap}")

if __name__ == "__main__":
    main()
