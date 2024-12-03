rule fq_lint:
    """Check that fastq files have valid format"""
    input:
        r1 = os.path.join(config["reads"], config["fastq_names_1"]),
        r2 = os.path.join(config["reads"], config["fastq_names_2"]),
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
        r1_all = expand(os.path.join(config["reads"], config["fastq_names_1"]), sample=SAMPLES),
        r2_all = expand(os.path.join(config["reads"], config["fastq_names_2"]), sample=SAMPLES),
        lint = expand(os.path.join(dir["logs"],"lint","{sample}.lint"), sample=SAMPLES)
    output:
        os.path.join(dir["stats"], "initial_summary.txt")
    conda:
        os.path.join(dir["env"], "seqkit.yaml")
    shell:
        """
        seqkit stats -T {input.r1_all} {input.r2_all} > {output}
        """
        
rule fastqc:
    """Perform fastqc on the input data"""
    input:
        r1 = expand(os.path.join(config["reads"], config["fastq_names_1"]), sample=SAMPLES),
        r2 = expand(os.path.join(config["reads"], config["fastq_names_2"]), sample=SAMPLES),
    output:
        fastqc_r1 = expand(os.path.join(dir["results"], "fastqc", "{sample}_R1_fastqc.html"), sample=SAMPLES),
        fastqc_r2 = expand(os.path.join(dir["results"], "fastqc", "{sample}_R2_fastqc.html"), sample=SAMPLES),
    params:
        outdir=os.path.join(dir["results"], "fastqc")
    conda:
        os.path.join(dir["env"], "fastqc.yaml")
    shell:
        """
        fastqc {input.r1} {input.r2} -o {params.outdir}
        """        

rule multiqc:
    """Get a report that consolidates the results for all samples"""
    input:
        expand(os.path.join(dir["results"], "fastqc", "{sample}_R1_fastqc.html"), sample=SAMPLES),
        expand(os.path.join(dir["results"], "fastqc", "{sample}_R2_fastqc.html"), sample=SAMPLES),
    output:
        os.path.join(dir["stats"], "multiqc_report.html")
    params:
        fastqc = os.path.join(dir["results"], "fastqc"),
        outdir = dir["stats"]
    conda:
        os.path.join(dir["env"], "multiqc.yaml")
    shell:
        """
        multiqc {params.fastqc} -o {params.outdir}
        #mv {params.fastqc}/multiqc_report.html {params.outdir} 
        """
