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


rule pB_step_1:
    # potential contamination type 1: Find 8 digits of primerB, trim to the left
    input:
        r1 = os.path.join(dir["output"], "fastp", "{sample}_R1_trimmed.fastq.gz"),
        r2 = os.path.join(dir["output"], "fastp", "{sample}_R2_trimmed.fastq.gz"), 
        ref = os.path.join(dir["db"], "last_8digits_primerB.fa"),
    output:
        out = os.path.join(dir["output"], "primerB", "step_1", "{sample}_R1.s1.out.fastq"),
        out2 = os.path.join(dir["output"], "primerB", "step_1", "{sample}_R2.s1.out.fastq"),
        stats = os.path.join(dir["output"], "primerB", "step_1", "{sample}_s1.stats"),
    params:
        config["qc"]["bbduk"]["s1"]
    threads: 8
    conda:
        os.path.join(dir["env"], "bbmap.yaml")
    log:
        os.path.join(dir["logs"], "primerB", "{sample}_left_pB.log")
    benchmark:
        os.path.join(dir["bench"], "primerB", "{sample}_left_pB.txt")
    shell:
        """
        bbduk.sh \
        in={input.r1} \
        in2={input.r2} \
        ref={input.ref} \
        out={output.out} \
        out2={output.out2} \
        stats={output.stats} \
        {params}
        """


rule pB_step_2:
    # potential contamination type 2: Find the first 8 digits of rc primerB, trim to the right
    input:
        r1 = os.path.join(dir["output"], "primerB", "step_1", "{sample}_R1.s1.out.fastq"),
        r2 = os.path.join(dir["output"], "primerB", "step_1", "{sample}_R2.s1.out.fastq"),
        ref = os.path.join(dir["db"], "first_8digits_primerB_rc.fa")
    output:
        out = os.path.join(dir["output"], "primerB", "step_2", "{sample}_R1.s2.out.fastq"),
        out2 = os.path.join(dir["output"], "primerB", "step_2", "{sample}_R2.s2.out.fastq"),
        stats = os.path.join(dir["output"], "primerB", "step_2", "{sample}_s2.stats"), 
    params:
        config["qc"]["bbduk"]["s2"]
    threads: 8
    conda:
        os.path.join(dir["env"], "bbmap.yaml")
    log:
        os.path.join(dir["logs"], "primerB", "{sample}_right_pB.log")
    benchmark:
        os.path.join(dir["bench"], "primerB", "{sample}_right_pB.txt")
    shell:
        """
        bbduk.sh \
        in={input.r1} \
        in2={input.r2} \
        ref={input.ref} \
        out={output.out} \
        out2={output.out2} \
        stats={output.stats} \
        {params}
        """


rule pB_step_3:
    # potential contamination type 3: Find the forward sequence primer and trim to the right
    input:
        r1 = os.path.join(dir["output"], "primerB", "step_2", "{sample}_R1.s2.out.fastq"),
        r2 = os.path.join(dir["output"], "primerB", "step_2", "{sample}_R2.s2.out.fastq"),
    output:
        out = os.path.join(dir["output"], "primerB", "step_3", "{sample}_R1.s3.out.fastq"),
        out2 = os.path.join(dir["output"], "primerB", "step_3", "{sample}_R2.s3.out.fastq"),
        stats = os.path.join(dir["output"], "primerB", "step_3", "{sample}_s3.stats"),
    params:
        config["qc"]["bbduk"]["s3_and_4"]
    threads: 8
    conda:
        os.path.join(dir["env"], "bbmap.yaml")
    log:
        os.path.join(dir["logs"], "primerB", "{sample}_pB_step_3.log")
    benchmark:
        os.path.join(dir["bench"], "primerB", "{sample}_pB_step_3.txt")
    shell:
        """
        bbduk.sh \
        in={input.r1} \
        in2={input.r2} \
        literal="CTGTCTCTTATACACATCT" \
        out={output.out} \
        out2={output.out2} \
        stats={output.stats} \
        {params}
        """


rule pB_step_4:
    # potential contamination type 4: Find the reverse sequence primer and trim to the right
    input:
        r1 = os.path.join(dir["output"], "primerB", "step_3", "{sample}_R1.s3.out.fastq"),
        r2 = os.path.join(dir["output"], "primerB", "step_3", "{sample}_R2.s3.out.fastq"),
    output:
        out = os.path.join(dir["output"], "primerB", "step_4", "{sample}_R1.s4.out.fastq"),
        out2 = os.path.join(dir["output"], "primerB", "step_4", "{sample}_R2.s4.out.fastq"),
        stats = os.path.join(dir["output"], "primerB", "step_4", "{sample}_s4.stats"),
    params:
        config["qc"]["bbduk"]["s3_and_4"]
    threads: 8
    conda:
        os.path.join(dir["env"], "bbmap.yaml")
    log:
        os.path.join(dir["logs"], "primerB", "{sample}_pB_step_4.log")
    benchmark:
        os.path.join(dir["bench"], "primerB", "{sample}_pB_step_4.txt")
    shell:
        """
        bbduk.sh \
        in={input.r1} \
        in2={input.r2} \
        literal="CTGTCTCTTCTACACATCT" \
        out={output.out} \
        out2={output.out2} \
        stats={output.stats} \
        {params}
        """


rule remove_vector_contamination:
    input:
        r1 = os.path.join(dir["output"], "primerB", "step_4", "{sample}_R1.s4.out.fastq"),
        r2 = os.path.join(dir["output"], "primerB", "step_4", "{sample}_R2.s4.out.fastq"),
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
