# Purpose:
Preprocessing module for Hecatomb
- Remove adapters, primers, and vector contamination
- Merge overlapping read pairs
- Mask host reference genome for viruses
- Remove host contamination
- Report summary statistics for each step


Workflow:

![hecatomb-qc-0217](https://github.com/user-attachments/assets/ea29a2c3-2e05-4d0c-84f1-33beb60dfe74)

# How to run:
- Git clone this repository
- Run from an environment with Snakemake version 8+, [mamba](https://anaconda.org/conda-forge/mamba), and [snakemake-executor-plugin-slurm](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html)

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

 - reads: specify paths to directories where paired-end fastq reads are. Can be:
   - A single directory: e.g., `reads=/scratch/sahlab/Megan/test_reads`
   - Multiple directories: e.g., `reads='["/dir1", "/dir2", "/dir3"]'`
   - Note: For samples with the same name in different directories, the first encountered will be used

 - output: specify path to directory where all outputs will be created

 - fastq_names_1: default is "{sample}_R1.fastq.gz"

 - fastq_names_2: default is "{sample}_R2.fastq.gz"
 
 - assembly_strategy: Choice of assembly approach (default: "coassembly")
   - "coassembly": Performs a single assembly with all samples combined (default)
   - "individual": Assembles each sample separately and then merges assemblies using Flye


### Example commands:

```
# Default coassembly with a single read directory
snakemake --profile ../profile/slurm/ --config reads=/scratch/sahlab/Megan/test_reads output=/scratch/sahlab/Megan/pipeline_test.out

# Individual assembly with Flye merging
snakemake --profile ../profile/slurm/ --config reads=/scratch/sahlab/Megan/test_reads output=/scratch/sahlab/Megan/pipeline_test.out assembly_strategy=individual

# Using multiple read directories
snakemake --profile ../profile/slurm/ --config reads='["/scratch/sahlab/Megan/test_reads1", "/scratch/sahlab/Megan/test_reads2"]' output=/scratch/sahlab/Megan/pipeline_test.out
```
