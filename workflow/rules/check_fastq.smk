rule fq_lint:
    """Check that fastq files have valid format"""
    input:
        r1 = get_read1_path,
        r2 = get_read2_path,
    output:
        lint=os.path.join(dir["logs"],"lint","{sample}.lint")
    conda:
        os.path.join(dir["env"], "fq.yaml")
    shell:
        """
        fq lint {input.r1} {input.r2} > {output.lint} 
        """

rule summary_stats:
    """Get some basic stats about the input data"""
    input:
        r1_all = lambda wildcards: get_all_read_files()['r1_files'],
        r2_all = lambda wildcards: get_all_read_files()['r2_files'],
        lint = expand(os.path.join(dir["logs"],"lint","{sample}.lint"), sample=SAMPLES)
    output:
        os.path.join(dir["stats"], "raw_input_data", "basic_stats.txt")
    conda:
        os.path.join(dir["env"], "seqkit.yaml")
    shell:
        """
        seqkit stats -T {input.r1_all} {input.r2_all} > {output}
        """
        
rule fastqc:
    """Perform fastqc on the input data"""
    input:
        r1 = get_read1_path,
        r2 = get_read2_path,
    output:
        fastqc_r1 = os.path.join(dir["results"], "output", "fastqc", "{sample}_R1_fastqc.html"),
        fastqc_r2 = os.path.join(dir["results"], "output", "fastqc", "{sample}_R2_fastqc.html"),
    params:
        outdir=os.path.join(dir["results"], "output", "fastqc")
    conda:
        os.path.join(dir["env"], "fastqc.yaml")
    shell:
        """
        fastqc {input.r1} {input.r2} -o {params.outdir}
        """        

rule multiqc:
    """Get a report that consolidates the results for all samples"""
    input:
        expand(os.path.join(dir["results"], "output", "fastqc", "{sample}_R1_fastqc.html"), sample=SAMPLES),
        expand(os.path.join(dir["results"], "output", "fastqc", "{sample}_R2_fastqc.html"), sample=SAMPLES),
    output:
        os.path.join(dir["stats"], "raw_input_data", "multiqc_report.html")
    params:
        fastqc = os.path.join(dir["results"], "output", "fastqc"),
        outdir = os.path.join(dir["stats"], "raw_input_data"),
    conda:
        os.path.join(dir["env"], "multiqc.yaml")
    shell:
        """
        multiqc {params.fastqc} -o {params.outdir}
        #mv {params.fastqc}/multiqc_report.html {params.outdir} 
        """
 
