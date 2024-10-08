rule remove_5prime_primer:
    """Round A/B step 01: Remove 5' primer."""
    input:
        r1 = os.path.join(config["reads"], config["fastq_names_1"]),
        r2 = os.path.join(config["reads"], config["fastq_names_2"]),
        primers=os.path.join(dir["db"],"primerB.fa")
    output:
        r1=temp(os.path.join(dir["temp"],"{sample}.R1.s1.fastq")),
        r2=temp(os.path.join(dir["temp"],"{sample}.R2.s1.fastq")),
        stats=os.path.join(dir["stats"],"{sample}_5prime_stats")
    benchmark:
        os.path.join(dir["bench"],"remove_5prime_primer.{sample}.txt")
    log:
        os.path.join(dir["log"],"remove_5prime_primer.{sample}.log")
    threads:
        24
    params:
        params = config["qc"]["bbduk"]["rm_5p"],
    conda:
        os.path.join(dir["env"], "bbmap.yaml")
    group:
        "roundAB"
    shell:
        """
        bbduk.sh \
            in={input.r1} \
            in2={input.r2} \
            ref={input.primers} \
            out={output.r1} \
            out2={output.r2} \
            stats={output.stats} \
            threads={threads} \
            {params.params} \
            2> {log}
        """


rule remove_3prime_contaminant:
    """Round A/B step 02: Remove 3' read through contaminant."""
    input:
        r1=os.path.join(dir["temp"],"{sample}.R1.s1.fastq"),
        r2=os.path.join(dir["temp"],"{sample}.R2.s1.fastq"),
        primers=os.path.join(dir["db"],"rc_primerB_ad6.fa")
    output:
        r1=temp(os.path.join(dir["temp"],"{sample}.R1.s2.fastq")),
        r2=temp(os.path.join(dir["temp"],"{sample}.R2.s2.fastq")),
        stats=os.path.join(dir["stats"],"{sample}_3prime_stats")
    benchmark:
        os.path.join(dir["bench"],"remove_3prime_contaminant.{sample}.txt")
    log:
        os.path.join(dir["log"],"remove_3prime_contaminant.{sample}.log")
    threads:
        24
    params:
        params = config["qc"]["bbduk"]["rm_3rt"]
    conda:
        os.path.join(dir["env"], "bbmap.yaml")
    group:
        "roundAB"
    shell:
        """
        bbduk.sh \
            in={input.r1} \
            in2={input.r2} \
            ref={input.primers} \
            out={output.r1} \
            out2={output.r2} \
            stats={output.stats} \
            {params.params} \
            threads={threads} \
            2> {log}
        """


rule remove_primer_free_adapter:
    """Round A/B step 03: Remove primer free adapter (both orientations)."""
    input:
        r1=os.path.join(dir["temp"],"{sample}.R1.s2.fastq"),
        r2=os.path.join(dir["temp"],"{sample}.R2.s2.fastq"),
        primers=os.path.join(dir["db"],"nebnext_adapters.fa")
    output:
        r1=temp(os.path.join(dir["temp"],"{sample}.R1.s3.fastq")),
        r2=temp(os.path.join(dir["temp"],"{sample}.R2.s3.fastq")),
        stats=os.path.join(dir["stats"],"{sample}_primer_free_adapter_stats")
    benchmark:
        os.path.join(dir["bench"],"remove_primer_free_adapter.{sample}.txt")
    log:
        os.path.join(dir["log"],"remove_primer_free_adapter.{sample}.log")
    threads:
        24
    params:
        params = config["qc"]["bbduk"]["neb"]
    conda:
        os.path.join(dir["env"], "bbmap.yaml")
    group:
        "roundAB"
    shell:
        """
        bbduk.sh \
            in={input.r1} \
            in2={input.r2} \
            ref={input.primers} \
            out={output.r1} \
            out2={output.r2} \
            stats={output.stats} \
            {params.params} \
            threads={threads} \
            2> {log}
        """


rule remove_adapter_free_primer:
    """Round A/B step 04: Remove adapter free primer (both orientations)."""
    input:
        r1=os.path.join(dir["temp"],"{sample}.R1.s3.fastq"),
        r2=os.path.join(dir["temp"],"{sample}.R2.s3.fastq"),
        primers=os.path.join(dir["db"],"rc_primerB_ad6.fa")
    output:
        r1=temp(os.path.join(dir["temp"],"{sample}.R1.s4.fastq")),
        r2=temp(os.path.join(dir["temp"],"{sample}.R2.s4.fastq")),
        stats=os.path.join(dir["stats"],"{sample}_adapter_free_primer_stats")
    benchmark:
        os.path.join(dir["bench"],"remove_adapter_free_primer.{sample}.txt")
    log:
        os.path.join(dir["log"],"remove_adapter_free_primer.{sample}.log")
    threads:
        24
    params:
        params = config["qc"]["bbduk"]["rm_afp"]
    conda:
        os.path.join(dir["env"], "bbmap.yaml")
    group:
        "roundAB"
    shell:
        """
        bbduk.sh \
            in={input.r1} \
            in2={input.r2} \
            ref={input.primers} \
            out={output.r1} \
            out2={output.r2} \
            stats={output.stats} \
            {params.params} \
            threads={threads} \
            2> {log}
        """


rule remove_vector_contamination:
    """Round A/B step 05: Vector contamination removal (PhiX + NCBI UniVecDB)"""
    input:
        r1=os.path.join(dir["temp"],"{sample}.R1.s4.fastq"),
        r2=os.path.join(dir["temp"],"{sample}.R2.s4.fastq"),
        primers=os.path.join(dir["db"],"vector_contaminants.fa")
    output:
        r1=temp(os.path.join(dir["temp"],"{sample}.R1.s5.fastq")),
        r2=temp(os.path.join(dir["temp"],"{sample}.R2.s5.fastq")),
        stats=os.path.join(dir["stats"],"{sample}_vector_contamination_stats")
    benchmark:
        os.path.join(dir["bench"],"remove_vector_contamination.{sample}.txt")
    log:
        os.path.join(dir["log"],"remove_vector_contamination.{sample}.log")
    threads:
        24
    params:
        params = config["qc"]["bbduk"]["rm_vc"]
    conda:
        os.path.join(dir["env"], "bbmap.yaml")
    group:
        "roundAB"
    shell:
        """
        bbduk.sh \
            in={input.r1} \
            in2={input.r2} \
            ref={input.primers} \
            out={output.r1} \
            out2={output.r2} \
            stats={output.stats} \
            {params.params} \
            threads={threads} \
            2> {log}
        """


rule remove_low_quality:
    """Round A/B step 06: Remove remaining low-quality bases and short reads."""
    input:
        r1=os.path.join(dir["temp"],"{sample}.R1.s5.fastq"),
        r2=os.path.join(dir["temp"],"{sample}.R2.s5.fastq"),
    output:
        r1=temp(os.path.join(dir["temp"],"{sample}.R1.s6.fastq")),
        r2=temp(os.path.join(dir["temp"],"{sample}.R2.s6.fastq")),
        stats=os.path.join(dir["stats"],"{sample}_low_quality_stats")
    benchmark:
        os.path.join(dir["bench"],"remove_low_quality.{sample}.txt")
    log:
        os.path.join(dir["log"],"remove_low_quality.{sample}.log")
    threads:
        24
    params:
        params = config["qc"]["bbduk"]["rm_lq"]
    conda:
        os.path.join(dir["env"], "bbmap.yaml")
    group:
        "roundAB"
    shell:
        """
        bbduk.sh \
            in={input.r1} \
            in2={input.r2} \
            out={output.r1} \
            out2={output.r2} \
            stats={output.stats} \
            threads={threads} \
            {params.params} \
            2> {log}
        """


rule zip_roundAB:
    """Zip the final trimmed reads for Round A/B"""
    input:
        r1=os.path.join(dir["temp"],"{sample}.R1.s6.fastq"),
        r2=os.path.join(dir["temp"],"{sample}.R2.s6.fastq"),
    output:
        r1=temp(os.path.join(dir["output"],"{sample}_R1.fastq.gz")),
        r2=temp(os.path.join(dir["output"],"{sample}_R2.fastq.gz")),
    benchmark:
        os.path.join(dir["bench"],"zip_roundAB.{sample}.txt")
    log:
        os.path.join(dir["log"],"zip_roundAB.{sample}.log")
    threads:
        24
    params:
        compression = config["qc"]["compression"]
    conda:
        os.path.join(dir["env"], "pigz.yaml")
    group:
        "roundAB"
    shell:
        """
        pigz -p {threads} -{params.compression} -c {input.r1} > {output.r1} 2> {log}
        pigz -p {threads} -{params.compression} -c {input.r2} > {output.r2} 2> {log}
        """
