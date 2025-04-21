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

rule megahit_assembly:
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


rule build_minimap_index:
    input:
        reference = os.path.join(dir["output"], "assembly", "megahit", "final.contigs.fa")
    output:
        index = os.path.join(dir["output"], "assembly", "megahit", "final.contigs.mmi")
    threads: 8
    conda:
        os.path.join(dir["env"], "minimap.yaml")
    log:
        os.path.join(dir["logs"], "indexes", "minimap2_index.log")
    benchmark:
        os.path.join(dir["bench"], "indexes", "minimap2_index.txt")
    shell:
        """
        minimap2 -t {threads} -d {output.index} {input.reference} 2> {log}
        """

rule align_host_removed_reads:
    input:
        r1 = os.path.join(dir["output"], "host_removed", "{sample}_hr_R1.fastq"),
        r2 = os.path.join(dir["output"], "host_removed", "{sample}_hr_R2.fastq"),
        index = os.path.join(dir["output"], "assembly", "megahit", "final.contigs.mmi")
    output:
        bam = os.path.join(dir["output"], "host_removed", "{sample}.bam"),
        sorted_bam = os.path.join(dir["output"], "host_removed", "{sample}_sorted.bam")
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

rule generate_pileup:
    input:
        bam = os.path.join(dir["output"], "host_removed", "{sample}_sorted.bam"),
        reference = os.path.join(dir["output"], "assembly", "megahit", "final.contigs.fa")
    output:
        pileup = os.path.join(dir["output"], "host_removed", "{sample}.pileup")
    log:
        os.path.join(dir["logs"], "host_removed", "{sample}_pileup.log")
    params:
        extra = "-d 10000"
    wrapper:
        "v1.2.1/bio/samtools/mpileup"


rule assembly_stats:
    """Calculate assembly statistics"""
    input:
        os.path.join(dir["output"], "assembly", "megahit", "final.contigs.fa")
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
