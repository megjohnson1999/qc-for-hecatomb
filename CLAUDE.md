# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This repository contains a Snakemake workflow for preprocessing metagenomic data focused on viral discovery. The workflow, named "qc-for-hecatomb", performs quality control, host removal, and assembly of high-throughput sequencing data to identify potential viral sequences.

Key components:
- Quality control and adapter removal from raw reads
- Primer and vector contamination removal
- Host genome masking and removal
- Read merging for overlapping pairs
- Metagenomic assembly using either co-assembly or individual assembly approaches
- Contig validation and chimeric contig detection

## Workflow Execution Commands

```bash
# Activate the Snakemake environment
source /ref/sahlab/software/miniforge3/bin/activate
conda activate snakemake_v8.20.3

# Run the workflow with default co-assembly strategy (single read directory)
cd workflow
snakemake --profile ../profile/slurm/ --config reads=/path/to/reads output=/path/to/output

# Run the workflow with individual assembly strategy
cd workflow
snakemake --profile ../profile/slurm/ --config reads=/path/to/reads output=/path/to/output assembly_strategy=individual

# Run the workflow with multiple read directories
cd workflow
snakemake --profile ../profile/slurm/ --config reads='["/path/to/reads1", "/path/to/reads2"]' output=/path/to/output
```

### Required Configuration Options
- `reads`: Directory or list of directories containing paired-end fastq files
  - Can be a single directory path or a list of paths in JSON format
  - Samples with the same name in different directories will use the first occurrence
- `output`: Directory where results will be written

### Optional Configuration Options
- `fastq_names_1`: Pattern for read 1 files (default: "{sample}_R1.fastq.gz")
- `fastq_names_2`: Pattern for read 2 files (default: "{sample}_R2.fastq.gz")
- `assembly_strategy`: Assembly approach ("coassembly" or "individual")

## Technical Architecture

The workflow is built with Snakemake v8+ and has several key components:

1. **Input Processing and Quality Control**
   - Raw read validation and statistics (check_fastq.smk)
   - Adapter and low-quality sequence removal with fastp
   - PrimerB and vector contamination removal with BBDuk

2. **Host Genome Processing**
   - Host genome masking using viral shreds to identify host regions with viral homology
   - Host read removal using minimap2 alignment to masked reference genome

3. **Read Processing**
   - Merging of overlapping reads with BBMerge
   - Processing of host-removed reads for assembly

4. **Assembly Strategies**
   - Co-assembly: Combines all samples using MEGAHIT
   - Individual assembly: Assembles each sample separately, then merges with Flye

5. **Contig Validation**
   - Mapping reads back to assembled contigs
   - Detection of potentially chimeric contigs
   - Coverage analysis and validation statistics

Each major workflow step runs in its own conda environment, defined in the workflow/envs/ directory.

## Computing Resources

The workflow is designed to run on a SLURM-based HPC environment. Key resource allocations:

- Default: 16GB RAM, 24 cores, 24-hour runtime
- Host processing: 64GB RAM
- Co-assembly: 300GB RAM, 66-hour runtime
- Flye assembly: 128GB RAM, 48-hour runtime

Resource configurations are defined in the profile/slurm/config.v8+.yaml file.