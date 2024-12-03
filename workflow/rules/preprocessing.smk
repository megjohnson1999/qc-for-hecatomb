rule fastp:
    """Run fastp to remove adapter sequences"""
    input:
        r1 = os.path.join(config["reads"], config["fastq_names_1"]),
        r2 = os.path.join(config["reads"], config["fastq_names_2"]),
    output:
        r1 = temp(os.path.join(dir["output"], "fastp", "{sample}_R1_trimmed.fastq.gz")),
        r2 = temp(os.path.join(dir["output"], "fastp", "{sample}_R2_trimmed.fastq.gz")),
        stats = os.path.join(dir["output"], "fastp", "{sample}.stats.json"),
        html = os.path.join(dir["output"], "fastp", "{sample}.stats.html"), 
    params:
        config["qc"]["fastp"]    
    threads: 12
    conda:
        os.path.join(dir["env"], "fastp.yaml")
    log:
        os.path.join(dir["logs"], "fastp", "fastp.{sample}.log")
    benchmark:
        os.path.join(dir["bench"], "fastp", "fastp.{sample}.txt")
    shell:
        """
        fastp \
        -i {input.r1} \
        -I {input.r2} \
        -o {output.r1} \
        -O {output.r2} \
        -j {output.stats} \
        -h {output.html} \
        --thread {threads} \
        {params} \
        &> {log}
        """

#rule(s) primerB:


rule remove_vector_contamination:
    input:
        r1 = os.path.join(dir["output"], "fastp", "{sample}_R1_trimmed.fastq.gz"),
        r2 = os.path.join(dir["output"], "fastp", "{sample}_R2_trimmed.fastq.gz"),
        contaminants = os.path.join(dir["db"], "vector_contaminants.fa")
    output:
        r1 = temp(os.path.join(dir["output"], "rm_vector_contamination", "{sample}_R1_rm_vc.fastq.gz")),
        r2 = temp(os.path.join(dir["output"], "rm_vector_contamination", "{sample}_R2_rm_vc.fastq.gz")),
        stats = os.path.join(dir["stats"], "rm_vector_contamination", "{sample}_rm_vc.stats"),
    params:
        config["qc"]["bbduk"]["rm_vc"]
    threads: 24
    conda:
        os.path.join(dir["env"], "bbmap.yaml")
    log:
        os.path.join(dir["logs"], "remove_vector_contamination", "{sample}_remove_vector_contamination.log")
    benchmark:
        os.path.join(dir["bench"], "remove_vector_contamination", "{sample}_remove_vector_contamination.txt")
    shell:
        """
        bbduk.sh \
        in={input.r1} \
        in2={input.r2} \
        ref={input.contaminants} \
        out={output.r1} \
        out2={output.r2} \
        stats={output.stats} \
        {params} \
        threads={threads} \
        2> {log}
        """

rule remove_low_quality:
    input:
        r1 = os.path.join(dir["output"], "rm_vector_contamination", "{sample}_R1_rm_vc.fastq.gz"),
        r2 = os.path.join(dir["output"], "rm_vector_contamination", "{sample}_R2_rm_vc.fastq.gz"),
    output:
        r1 = temp(os.path.join(dir["output"], "rm_low_quality", "{sample}_R1_rm_lq.fastq.gz")),
        r2 = temp(os.path.join(dir["output"], "rm_low_quality", "{sample}_R2_rm_lq.fastq.gz")),
        stats = os.path.join(dir["stats"], "rm_low_quality", "{sample}_low_quality.stats"),
    params:
        config["qc"]["bbduk"]["rm_lq"]
    threads: 24
    conda:
        os.path.join(dir["env"], "bbmap.yaml")
    log:
        os.path.join(dir["logs"], "remove_low_quality", "{sample}_remove_low_quality.log")
    benchmark:
        os.path.join(dir["bench"], "remove_low_quality", "{sample}_remove_low_quality.txt")
    shell:
        """
        bbduk.sh \
        in={input.r1} \
        in2={input.r2} \
        out={output.r1} \
        out2={output.r2} \
        stats={output.stats} \
        threads={threads} \
        {params} \
        2> {log}
        """


rule bbmerge:
    input:
        r1 = os.path.join(dir["output"], "rm_low_quality", "{sample}_R1_rm_lq.fastq.gz"),
        r2 = os.path.join(dir["output"], "rm_low_quality", "{sample}_R2_rm_lq.fastq.gz"),
    output:
        merged = os.path.join(dir["output"], "bbmerge", "{sample}_merged.fastq.gz"),
        unmerged1 = os.path.join(dir["output"], "bbmerge", "{sample}_R1_unmerged.fastq.gz"),
        unmerged2 = os.path.join(dir["output"], "bbmerge", "{sample}_R2_unmerged.fastq.gz"),
    params:
        config["qc"]["bbmerge"]
    threads: 24
    conda:
        os.path.join(dir["env"], "bbmap.yaml")
    log:
        os.path.join(dir["logs"], "bbmerge", "{sample}_bbmerge.log")
    benchmark:
        os.path.join(dir["bench"], "bbmerge", "{sample}_bbmerge.txt")
    shell:
        """
        bbmerge.sh \
        in1={input.r1} \
        in2={input.r2} \
        out={output.merged} \
        outu1={output.unmerged1} \
        outu2={output.unmerged2} \
        {params} \
        &> {log}
        """
