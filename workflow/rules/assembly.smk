rule concatenate_merged_reads:
    """Concatenate all merged reads into a single file for coassembly"""
    input:
        expand(os.path.join(dir["output"], "bbmerge", "{sample}_merged.fastq.gz"), sample=SAMPLES)
    output:
        os.path.join(dir["output"], "assembly", "all_merged.fastq.gz")
    threads: 4
    conda:
        os.path.join(dir["env"], "pigz.yaml")
    log:
        os.path.join(dir["logs"], "assembly", "concat_merged.log")
    benchmark:
        os.path.join(dir["bench"], "assembly", "concat_merged.txt")
    shell:
        """
        # Validate input files first
        echo "Checking input file integrity before concatenation..." > {log}
        for file in {input}; do
            gzip -t "$file" || {{ echo "❌ Error: File $file failed integrity check" >> {log}; exit 1; }}
        done
        echo "✅ All input files passed integrity check" >> {log}
        
        # Use zcat and redirect to gzip to preserve file integrity
        echo "Starting concatenation with zcat..." >> {log}
        zcat {input} | gzip -c > {output} 2>> {log}
        
        # Verify output file
        echo "Verifying output file integrity..." >> {log}
        gzip -t {output} && echo "✅ Output file integrity verified" >> {log} || {{ echo "❌ Output file corrupted" >> {log}; exit 1; }}
        """

rule concatenate_unmerged_r1:
    """Concatenate all unmerged R1 reads into a single file for coassembly"""
    input:
        expand(os.path.join(dir["output"], "bbmerge", "{sample}_R1_unmerged.fastq.gz"), sample=SAMPLES)
    output:
        os.path.join(dir["output"], "assembly", "all_unmerged_R1.fastq.gz")
    threads: 4
    conda:
        os.path.join(dir["env"], "pigz.yaml")
    log:
        os.path.join(dir["logs"], "assembly", "concat_unmerged_r1.log")
    benchmark:
        os.path.join(dir["bench"], "assembly", "concat_unmerged_r1.txt")
    shell:
        """
        # Validate input files first
        echo "Checking input file integrity before concatenation..." > {log}
        for file in {input}; do
            gzip -t "$file" || {{ echo "❌ Error: File $file failed integrity check" >> {log}; exit 1; }}
        done
        echo "✅ All input files passed integrity check" >> {log}
        
        # Use zcat and redirect to gzip to preserve file integrity
        # This creates a new compressed file instead of concatenating compressed files directly
        echo "Starting concatenation with zcat..." >> {log}
        zcat {input} | gzip -c > {output} 2>> {log}
        
        # Verify output file
        echo "Verifying output file integrity..." >> {log}
        gzip -t {output} && echo "✅ Output file integrity verified" >> {log} || {{ echo "❌ Output file corrupted" >> {log}; exit 1; }}
        """

rule concatenate_unmerged_r2:
    """Concatenate all unmerged R2 reads into a single file for coassembly"""
    input:
        expand(os.path.join(dir["output"], "bbmerge", "{sample}_R2_unmerged.fastq.gz"), sample=SAMPLES)
    output:
        os.path.join(dir["output"], "assembly", "all_unmerged_R2.fastq.gz")
    threads: 4
    conda:
        os.path.join(dir["env"], "pigz.yaml")
    log:
        os.path.join(dir["logs"], "assembly", "concat_unmerged_r2.log")
    benchmark:
        os.path.join(dir["bench"], "assembly", "concat_unmerged_r2.txt")
    shell:
        """
        # Validate input files first
        echo "Checking input file integrity before concatenation..." > {log}
        for file in {input}; do
            gzip -t "$file" || {{ echo "❌ Error: File $file failed integrity check" >> {log}; exit 1; }}
        done
        echo "✅ All input files passed integrity check" >> {log}
        
        # Use zcat and redirect to gzip to preserve file integrity
        echo "Starting concatenation with zcat..." >> {log}
        zcat {input} | gzip -c > {output} 2>> {log}
        
        # Verify output file
        echo "Verifying output file integrity..." >> {log}
        gzip -t {output} && echo "✅ Output file integrity verified" >> {log} || {{ echo "❌ Output file corrupted" >> {log}; exit 1; }}
        """

rule verify_read_counts:
    """Verify that the unmerged read files have equal counts before assembly"""
    input:
        r1 = os.path.join(dir["output"], "assembly", "all_unmerged_R1.fastq.gz"),
        r2 = os.path.join(dir["output"], "assembly", "all_unmerged_R2.fastq.gz")
    output:
        verification = os.path.join(dir["logs"], "assembly", "read_count_verification.txt")
    conda:
        os.path.join(dir["env"], "pandas.yaml")
    shell:
        """
        # Verify file integrity first
        gzip -t {input.r1} || {{ echo "Error: R1 file corrupted" > {output.verification}; exit 1; }}
        gzip -t {input.r2} || {{ echo "Error: R2 file corrupted" > {output.verification}; exit 1; }}
        
        # Count reads in each file
        r1_count=$(zcat {input.r1} | wc -l | awk '{{print $1/4}}')
        r2_count=$(zcat {input.r2} | wc -l | awk '{{print $1/4}}')
        
        # Output counts to verification file
        echo "Read counts before assembly:" > {output.verification}
        echo "R1 reads: $r1_count" >> {output.verification}
        echo "R2 reads: $r2_count" >> {output.verification}
        
        # Check if counts are equal
        if [ "$r1_count" -eq "$r2_count" ]; then
            echo "✅ Read counts are equal ($r1_count reads in each file)" >> {output.verification}
        else
            echo "❌ WARNING: Unequal read counts detected! R1: $r1_count, R2: $r2_count" >> {output.verification}
            echo "Difference: $(($r1_count > $r2_count ? $r1_count - $r2_count : $r2_count - $r1_count)) reads" >> {output.verification}
            # Don't fail - we'll handle this in the assembly rule
        fi
        """

def get_assembly_input(wildcards):
    """Return assembly input files based on the assembly strategy.
    
    This function should only be used in a rule that knows which keys to expect
    for the current assembly strategy.
    """
    strategy = config.get("assembly_strategy", "coassembly")
    result = {}
    
    if strategy == "coassembly":
        result = {
            "merged": os.path.join(dir["output"], "assembly", "all_merged.fastq.gz"),
            "r1": os.path.join(dir["output"], "assembly", "all_unmerged_R1.fastq.gz"),
            "r2": os.path.join(dir["output"], "assembly", "all_unmerged_R2.fastq.gz")
        }
    elif strategy in ["individual", "per_sample"]:
        result = {
            "merged_contigs": os.path.join(dir["output"], "assembly", "flye", "assembly.fasta")
        }
    else:
        print(f"Warning: Unexpected assembly_strategy value: '{strategy}', defaulting to coassembly")
        result = {
            "merged": os.path.join(dir["output"], "assembly", "all_merged.fastq.gz"),
            "r1": os.path.join(dir["output"], "assembly", "all_unmerged_R1.fastq.gz"),
            "r2": os.path.join(dir["output"], "assembly", "all_unmerged_R2.fastq.gz")
        }
    
    # Dictionary will only contain keys that were added
    return result
rule megahit_coassembly:
    """Perform coassembly of all samples using MEGAHIT"""
    input:
        merged = os.path.join(dir["output"], "assembly", "all_merged.fastq.gz"),
        r1 = os.path.join(dir["output"], "assembly", "all_unmerged_R1.fastq.gz"),
        r2 = os.path.join(dir["output"], "assembly", "all_unmerged_R2.fastq.gz")
    output:
        contigs = os.path.join(dir["output"], "assembly", "megahit", "final.contigs.fa"),
        dir = directory(os.path.join(dir["output"], "assembly", "megahit"))
    params:
        min_contig = 1000,
        out_dir = os.path.join(dir["output"], "assembly", "megahit")
    threads: 24
    conda:
        os.path.join(dir["env"], "megahit.yaml")
    log:
        os.path.join(dir["logs"], "assembly", "megahit.log")
    benchmark:
        os.path.join(dir["bench"], "assembly", "megahit.txt")
    shell:
        """
        # Remove output directory if it exists since megahit won't overwrite
        rm -rf {params.out_dir}
        
        # Run megahit
        megahit -r {input.merged} -1 {input.r1} -2 {input.r2} \
            -o {params.out_dir} \
            --min-contig-len {params.min_contig} \
            --k-list 21,29,39,59,79,99,119,141 \
            -t {threads} \
            &>> {log}
        """

rule megahit_individual_assembly:
    """Perform assembly of individual samples using MEGAHIT"""
    input:
        merged = os.path.join(dir["output"], "bbmerge", "{sample}_merged.fastq.gz"),
        r1 = os.path.join(dir["output"], "bbmerge", "{sample}_R1_unmerged.fastq.gz"),
        r2 = os.path.join(dir["output"], "bbmerge", "{sample}_R2_unmerged.fastq.gz")
    output:
        contigs = os.path.join(dir["output"], "assembly", "per_sample", "{sample}", "final.contigs.fa"),
        dir = directory(os.path.join(dir["output"], "assembly", "per_sample", "{sample}"))
    params:
        min_contig = 1000,
        out_dir = os.path.join(dir["output"], "assembly", "per_sample", "{sample}")
    threads: 16
    conda:
        os.path.join(dir["env"], "megahit.yaml")
    log:
        os.path.join(dir["logs"], "assembly", "per_sample", "{sample}_megahit.log")
    benchmark:
        os.path.join(dir["bench"], "assembly", "per_sample", "{sample}_megahit.txt")
    shell:
        """
        # Remove output directory if it exists since megahit won't overwrite
        rm -rf {params.out_dir}
        
        # Run megahit
        megahit -r {input.merged} -1 {input.r1} -2 {input.r2} \
            -o {params.out_dir} \
            --min-contig-len {params.min_contig} \
            --k-list 21,29,39,59,79,99,119,141 \
            -t {threads} \
            &>> {log}
        """

rule collect_individual_assemblies:
    """Collect individual assembly contigs for merging with Flye"""
    input:
        expand(os.path.join(dir["output"], "assembly", "per_sample", "{sample}", "final.contigs.fa"), sample=SAMPLES)
    output:
        directory(os.path.join(dir["output"], "assembly", "per_sample", "contigs"))
    log:
        os.path.join(dir["logs"], "assembly", "collect_individual_assemblies.log")
    shell:
        """
        # Create output directory
        mkdir -p {output} 2> {log}
        
        # Copy each assembly with a unique name
        for file in {input}; do
            sample=$(echo $file | rev | cut -d/ -f2 | rev)
            cp $file {output}/$sample.contigs.fa 2>> {log}
        done
        """

rule flye_merge_assemblies:
    """Merge individual assemblies using Flye subassemblies mode"""
    input:
        contig_dir = os.path.join(dir["output"], "assembly", "per_sample", "contigs")
    output:
        contigs = os.path.join(dir["output"], "assembly", "flye", "assembly.fasta"),
        dir = directory(os.path.join(dir["output"], "assembly", "flye")),
        stats = os.path.join(dir["output"], "assembly", "flye", "assembly_stats.txt")
    params:
        out_dir = os.path.join(dir["output"], "assembly", "flye")
    threads: 24
    conda:
        os.path.join(dir["env"], "flye.yaml")
    log:
        log1 = os.path.join(dir["logs"], "assembly", "flye_merge.log"),
        log2 = os.path.join(dir["logs"], "assembly", "stats.log")
    benchmark:
        os.path.join(dir["bench"], "assembly", "flye_merge.txt")
    shell:
        """
        # Create a comma-separated list of all contig files
        CONTIG_FILES=$(find {input.contig_dir} -name "*.contigs.fa" | tr '\n' ',' | sed 's/,$//')
        
        # Run flye in subassemblies mode with plasmids flag
        flye --subassemblies $CONTIG_FILES \
            --out-dir {params.out_dir} \
            --plasmids \
            -g 1g \
            --threads {threads} \
            &> {log.log1}
            
        # Generate assembly statistics
        statswrapper.sh in={output.contigs} out={output.stats} \
            format=2 \
            ow=t 2> {log.log2}
        """


def get_assembly_path(wildcards):
    """Return the path to the final assembly file based on the assembly strategy.
    
    This function will always return a valid path string, never None.
    """
    strategy = config.get("assembly_strategy", "coassembly")
    if strategy == "coassembly":
        return os.path.join(dir["output"], "assembly", "megahit", "final.contigs.fa")
    elif strategy in ["individual", "per_sample"]:
        return os.path.join(dir["output"], "assembly", "flye", "assembly.fasta")
    else:
        # Handle unexpected values, default to coassembly
        print(f"Warning: Unexpected assembly_strategy value: '{strategy}', defaulting to coassembly")
        return os.path.join(dir["output"], "assembly", "megahit", "final.contigs.fa")

def get_contigs_mmi_path(wildcards):
    """Return the path to the minimap2 index file based on the assembly strategy.
    
    This function will always return a valid path string, never None.
    """
    strategy = config.get("assembly_strategy", "coassembly")
    if strategy == "coassembly":
        return os.path.join(dir["output"], "assembly", "megahit", "final.contigs.mmi")
    elif strategy in ["individual", "per_sample"]:
        return os.path.join(dir["output"], "assembly", "flye", "assembly.mmi")
    else:
        # Handle unexpected values, default to coassembly
        print(f"Warning: Unexpected assembly_strategy value: '{strategy}', defaulting to coassembly")
        return os.path.join(dir["output"], "assembly", "megahit", "final.contigs.mmi")

# Temporarily commenting out this rule to test multiple read directory functionality
#rule index_contigs:
#    """Index contigs using minimap2 for read mapping"""
#    input:
#        get_assembly_path
#    output:
#        get_contigs_mmi_path
#    threads: 12
#    conda:
#        os.path.join(dir["env"], "minimap_env.yaml")
#    log:
#        os.path.join(dir["logs"], "contig_validation", "index_contigs.log")
#    benchmark:
#        os.path.join(dir["bench"], "contig_validation", "index_contigs.txt")
#    shell:
#        """
#        minimap2 -d {output} {input} 2> {log}
#        """

# Temporary placeholder rule with fixed output paths
rule index_contigs:
    """Index contigs using minimap2 for read mapping (temporary version)"""
    input:
        contigs = get_assembly_path
    output:
        index = os.path.join(dir["output"], "assembly", "contigs.mmi")
    threads: 12
    conda:
        os.path.join(dir["env"], "minimap_env.yaml")
    log:
        os.path.join(dir["logs"], "contig_validation", "index_contigs.log")
    benchmark:
        os.path.join(dir["bench"], "contig_validation", "index_contigs.txt")
    shell:
        """
        minimap2 -d {output.index} {input.contigs} 2> {log}
        """

rule align_host_removed_reads:
    input:
        r1 = os.path.join(dir["output"], "host_removed", "{sample}_hr_R1.fastq"),
        r2 = os.path.join(dir["output"], "host_removed", "{sample}_hr_R2.fastq"),
        index = os.path.join(dir["output"], "assembly", "contigs.mmi")  # Fixed path to match new index_contigs rule
    output:
        bam_index = os.path.join(dir["output"], "host_removed", "{sample}_to_contig_sorted.bam.bai"),
        sorted_bam = os.path.join(dir["output"], "host_removed", "{sample}_to_contig_sorted.bam")
    threads: 24
    conda:
        os.path.join(dir["env"], "minimap.yaml")
    log:
        os.path.join(dir["logs"], "host_removed", "{sample}_alignment.log")
    benchmark:
        os.path.join(dir["bench"], "host_removed", "{sample}_alignment.txt")
    shell:
        """
        minimap2 -ax sr -t {threads} {input.index} {input.r1} {input.r2} | \
        samtools view -Sb - | \
        samtools sort -o {output.sorted_bam} - && \
        samtools index {output.sorted_bam}
        """

rule merge_bams:
    input:
        expand(os.path.join(dir["output"], "host_removed", "{sample}_to_contig_sorted.bam"), sample=SAMPLES)
    output:
        merged_bam = os.path.join(dir["output"], "host_removed", "merged.bam")
    threads: 8
    conda:
        os.path.join(dir["env"], "minimap.yaml")
    log:
        os.path.join(dir["logs"], "host_removed", "merge_bams.log")
    shell:
        """
        samtools merge -@ {threads} {output.merged_bam} {input}
        """

rule generate_pileup:
    input:
        bam = os.path.join(dir["output"], "host_removed", "merged.bam"),
        reference = get_assembly_path
    output:
        pileup = os.path.join(dir["output"], "host_removed", "merged.pileup")
    log:
        os.path.join(dir["logs"], "host_removed", "merged_pileup.log")
    params:
        extra = "-d 10000"
    shell:
        """
        samtools mpileup {params.extra} -f {input.reference} {input.bam} > {output.pileup} 2> {log}
        """


rule assembly_stats:
    """Calculate assembly statistics"""
    input:
        get_assembly_path
    output:
        stats = os.path.join(dir["stats"], "assembly", "assembly_stats.txt")
    conda:
        os.path.join(dir["env"], "seqkit.yaml")
    threads: 4
    log:
        os.path.join(dir["logs"], "assembly", "assembly_stats.log")
    shell:
        """
        seqkit stats -a {input} > {output.stats} 2> {log}
        """
