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
        cat {input} > {output} 2> {log}
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
        cat {input} > {output} 2> {log}
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
        cat {input} > {output} 2> {log}
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
        # Remove output directory if it exists since MEGAHIT won't overwrite
        rm -rf {params.out_dir}
        
        megahit --12 {input.merged} -1 {input.r1} -2 {input.r2} \
            -o {params.out_dir} \
            --min-contig-len {params.min_contig} \
            --k-list 21,29,39,59,79,99,119,141 \
            -t {threads} \
            &> {log}
        """

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
