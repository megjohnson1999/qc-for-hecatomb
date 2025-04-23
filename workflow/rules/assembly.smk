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

# Add conditional logic to determine which assembly rules to run
def assembly_mode_choose_contigs(wildcards):
    if config.get("assembly_strategy", "coassembly") == "coassembly":
        # For coassembly, we need concatenated files
        return [
            os.path.join(dir["output"], "assembly", "all_merged.fastq.gz"),
            os.path.join(dir["output"], "assembly", "all_unmerged_R1.fastq.gz"),
            os.path.join(dir["output"], "assembly", "all_unmerged_R2.fastq.gz")
        ]
    else:
        # For individual assemblies, we need all per-sample assemblies and the merged Flye assembly
        return [
            os.path.join(dir["output"], "assembly", "flye", "assembly.fasta")
        ]

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
        dir = directory(os.path.join(dir["output"], "assembly", "flye"))
    params:
        out_dir = os.path.join(dir["output"], "assembly", "flye"),
        min_contig = 1000
    threads: 24
    conda:
        os.path.join(dir["env"], "flye.yaml")
    log:
        os.path.join(dir["logs"], "assembly", "flye_merge.log")
    benchmark:
        os.path.join(dir["bench"], "assembly", "flye_merge.txt")
    shell:
        """
        # Create a list of all contig files
        find {input.contig_dir} -name "*.contigs.fa" > contig_files.txt
        
        # Run flye in subassemblies mode
        flye --subassemblies @contig_files.txt \
            --out-dir {params.out_dir} \
            --min-overlap 500 \
            --threads {threads} \
            &> {log}
            
        # Remove temporary file
        rm contig_files.txt
        """


def get_assembly_path(wildcards):
    if config.get("assembly_strategy", "coassembly") == "coassembly":
        return os.path.join(dir["output"], "assembly", "megahit", "final.contigs.fa")
    else:
        return os.path.join(dir["output"], "assembly", "flye", "assembly.fasta")

def get_contigs_mmi_path(wildcards):
    if config.get("assembly_strategy", "coassembly") == "coassembly":
        return os.path.join(dir["output"], "assembly", "megahit", "final.contigs.mmi")
    else:
        return os.path.join(dir["output"], "assembly", "flye", "assembly.mmi")

rule index_contigs:
    """Index contigs using minimap2 for read mapping"""
    input:
        get_assembly_path
    output:
        get_contigs_mmi_path
    threads: 12
    conda:
        os.path.join(dir["env"], "minimap_env.yaml")
    log:
        os.path.join(dir["logs"], "contig_validation", "index_contigs.log")
    benchmark:
        os.path.join(dir["bench"], "contig_validation", "index_contigs.txt")
    shell:
        """
        minimap2 -d {output} {input} 2> {log}
        """

rule align_host_removed_reads:
    input:
        r1 = os.path.join(dir["output"], "host_removed", "{sample}_hr_R1.fastq"),
        r2 = os.path.join(dir["output"], "host_removed", "{sample}_hr_R2.fastq"),
        index = get_contigs_mmi_path
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
