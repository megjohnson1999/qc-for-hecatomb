# How to run:
Run from an environment with Snakemake version 8+, [mamba](https://anaconda.org/conda-forge/mamba), and [snakemake-executor-plugin-slurm](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html)

```
source /ref/sahlab/software/miniforge3/bin/activate
conda activate snakemake_v8.20.3
```

```
cd workflow
snakemake --profile ../profile/slurm/ --config [options]
```

## Note:

If your fastq files have any suffixes other than _R1.fastq.gz and _R2.fastq.gz, need to specify this on the command line using the fastq_names_1 and fastq_names_2 options

# Options:

 - reads: specify path to directory where paired-end fastq reads are

 - output: specify path to directory where all outputs will be created

 - fastq_names_1: default is "{sample}_R1.fastq.gz"

 - fastq_names_2: default is "{sample}_R2.fastq.gz"


### Example command:

```
snakemake --profile ../profile/slurm/ --config reads=/scratch/sahlab/Megan/test_reads output=/scratch/sahlab/Megan/primer_removal_test.out
```

# Rules and databases:
- remove_5prime_primer
  - db: primerB.fa
 
- remove_3prime_contaminant
  - db: rc_primerB_ad6.fa

- remove_primer_free_adapter
  - db: nebnext_adapters.fa
  
- remove_adapter_free_primer
  - db: rc_primerB_ad6.fa
  
- remove_vector_contamination (remove contaminants like PhiX)
  - db: vector_contaminants.fa 

### See workflow/databases to look at all database files
Will need to be updated with new Illumina adapter or index sequences

# Outputs:

The output directory should contain a directory called "results" with a subdirectory called "stats". The stats directory contains files with statistics for the primers/contaminants removed from each sample.
Currently, the stats for each rule for each sample are all in separate files. I can change this if there's a different way that would make it easier to work with.

- see example output at /scratch/sahlab/Megan/rc2_subset/primer_removal_test.out/results/stats/
