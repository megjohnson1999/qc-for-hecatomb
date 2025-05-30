# QC for Hecatomb Workflow Summary

## Overview
The QC for Hecatomb workflow is a comprehensive bioinformatics pipeline designed for quality control and preprocessing of metagenomic sequencing data, with a focus on viral sequence analysis. It supports both coassembly (combining all samples) and individual assembly strategies.

## Key Features
- **Multiple Input Directory Support**: Can process samples from one or more directories
- **Flexible Assembly Strategies**: Supports both coassembly and per-sample assembly approaches
- **Host Contamination Removal**: Removes human sequences using masked reference genomes
- **Viral Database Integration**: Uses viral genomes for masking host reference
- **Comprehensive QC**: Multiple quality control steps including adapter removal, primer trimming, and vector contamination removal

## Workflow Steps

### 1. Input Validation & Basic Statistics
- **fq_lint**: Validates FASTQ file format integrity
- **summary_stats**: Generates basic statistics on raw input data using seqkit
- **fastqc**: Performs quality control analysis on raw reads
- **multiqc**: Consolidates FastQC reports across all samples

### 2. Preprocessing & Quality Control
- **fastp**: Removes adapters, trims low-quality bases, filters by length and quality
  - Quality threshold: Q15
  - Minimum length: 90bp
  - Poly-G/X trimming enabled
- **Primer B Removal** (2 steps):
  - Step 1: Removes primer B sequences from the left (5' end)
  - Step 2: Removes reverse complement primer B from the right (3' end)
- **remove_vector_contamination**: Removes vector contamination using BBDuk
- **host_removal**: Maps reads to masked human reference genome and extracts unmapped reads

### 3. Read Merging
- **bbmerge**: Merges overlapping paired-end reads
  - Minimum overlap: 20bp
  - Produces merged reads and unmerged R1/R2 files

### 4. Assembly
Supports two strategies configured via `assembly_strategy`:

**Coassembly Mode** (default):
- Concatenates all samples' reads
- Performs single MEGAHIT assembly on combined data
- Minimum contig length: 1000bp

**Individual Assembly Mode**:
- Assembles each sample independently with MEGAHIT
- Merges individual assemblies using Flye in subassemblies mode
- Handles empty or failed assemblies gracefully

### 5. Contig Validation
- **index_contigs**: Creates minimap2 index of assembled contigs
- **map_reads_to_contigs**: Maps host-removed reads back to contigs
- **generate_per_sample_coverage**: Calculates coverage in 1kb windows
- **detect_chimeric_contigs**: Identifies potentially chimeric contigs based on coverage patterns
- **contig_validation_summary**: Generates comprehensive validation report

### 6. Host Genome Masking (preparation step)
- **vir_shred**: Shreds viral database into overlapping fragments
- **map_shreds**: Maps viral fragments to host genome
- **mask_host**: Masks regions of host genome similar to viral sequences

## Output Structure
```
output/
├── results/
│   ├── output/
│   │   ├── fastqc/           # FastQC reports
│   │   ├── qc/               # Intermediate QC files
│   │   ├── host_removed/     # Host-removed reads
│   │   ├── bbmerge/          # Merged/unmerged reads
│   │   └── assembly/         # Assembly results
│   └── stats/
│       ├── raw_input_data/   # Input statistics
│       ├── qc/               # QC statistics and plots
│       ├── assembly/         # Assembly statistics
│       └── contig_validation/ # Validation results
├── logs/                     # Process logs
├── benchmarks/               # Performance metrics
└── temp/                     # Temporary files
```

## Key Statistics Generated
- Read counts at each processing step
- Adapter/primer removal statistics
- Host contamination removal rates
- Assembly metrics (N50, total length, contig count)
- Read mapping rates to contigs
- Chimeric contig detection results

## Configuration
The workflow is configured through `config/config.yaml`:
- Input/output directories
- Sample naming patterns
- Assembly strategy selection
- Reference genome paths
- QC parameter thresholds

## Dependencies
Uses conda environments for tool management:
- fastp, bbmap suite (BBDuk, BBMerge, etc.)
- minimap2, samtools
- MEGAHIT, Flye assemblers
- FastQC, MultiQC
- R (for visualization)
- Python (pandas for statistics)