rule fastp:
    """Run fastp to remove adapter sequences"""
    input:
        r1 = os.path.join(config["reads"], config["fastq_names_1"]),
        r2 = os.path.join(config["reads"], config["fastq_names_2"]),
    output:
        r1 = temp(os.path.join(dir["output"], "qc", "fastp", "{sample}_R1_trimmed.fastq.gz")),
        r2 = temp(os.path.join(dir["output"], "qc", "fastp", "{sample}_R2_trimmed.fastq.gz")),
        stats = temp(os.path.join(dir["output"], "qc", "fastp", "{sample}.stats.json")),
        html = temp(os.path.join(dir["output"], "qc", "fastp", "{sample}.stats.html")), 
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


rule fastp_summary:
    input:
        expand(os.path.join(dir["output"], "qc", "fastp", "{sample}.stats.json"), sample = SAMPLES)
    output:
        os.path.join(dir["stats"], "qc", "fastp_stats.tsv")
    conda:
        os.path.join(dir["env"], "pandas.yaml")
    script:
        "../scripts/fastp_stats.py"


rule pB_step_1:
    # potential contamination type 1: Find 8 digits of primerB, trim to the left
    input:
        r1 = os.path.join(dir["output"], "qc", "fastp", "{sample}_R1_trimmed.fastq.gz"),
        r2 = os.path.join(dir["output"], "qc", "fastp", "{sample}_R2_trimmed.fastq.gz"), 
        ref = os.path.join(dir["db"], "primerB.fa"),
    output:
        out = temp(os.path.join(dir["output"], "qc", "step_1", "{sample}_R1.s1.out.fastq.gz")),
        out2 = temp(os.path.join(dir["output"], "qc", "step_1", "{sample}_R2.s1.out.fastq.gz")),
        stats = os.path.join(dir["stats"], "qc", "step_1", "{sample}_s1.stats"),
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
        r1 = os.path.join(dir["output"], "qc", "step_1", "{sample}_R1.s1.out.fastq.gz"),
        r2 = os.path.join(dir["output"], "qc", "step_1", "{sample}_R2.s1.out.fastq.gz"),
        ref = os.path.join(dir["db"], "primerB_rc.fa")
    output:
        out = temp(os.path.join(dir["output"], "qc", "step_2", "{sample}_R1.s2.out.fastq.gz")),
        out2 = temp(os.path.join(dir["output"], "qc", "step_2", "{sample}_R2.s2.out.fastq.gz")),
        stats = os.path.join(dir["stats"], "qc", "step_2", "{sample}_s2.stats"), 
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


rule remove_vector_contamination:
    input:
        r1 = os.path.join(dir["output"], "qc", "step_2", "{sample}_R1.s2.out.fastq.gz"),
        r2 = os.path.join(dir["output"], "qc", "step_2", "{sample}_R2.s2.out.fastq.gz"),
        contaminants = os.path.join(dir["db"], "vector_contaminants.fa")
    output:
        r1 = os.path.join(dir["output"], "qc", "rm_vector_contamination", "{sample}_R1_rm_vc.fastq.gz"),
        r2 = os.path.join(dir["output"], "qc", "rm_vector_contamination", "{sample}_R2_rm_vc.fastq.gz"),
        stats = os.path.join(dir["stats"], "qc", "rm_vector_contamination", "{sample}_rm_vc.stats"),
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
        

rule bbmerge:
    input:
        r1 = os.path.join(dir["output"], "qc", "rm_vector_contamination", "{sample}_R1_rm_vc.fastq.gz"),
        r2 = os.path.join(dir["output"], "qc", "rm_vector_contamination", "{sample}_R2_rm_vc.fastq.gz"),
    output:
        merged = os.path.join(dir["output"], "bbmerge", "{sample}_merged.fastq.gz"),
        unmerged1 = os.path.join(dir["output"], "bbmerge", "{sample}_R1_unmerged.fastq.gz"),
        unmerged2 = os.path.join(dir["output"], "bbmerge", "{sample}_R2_unmerged.fastq.gz"),
        hist = os.path.join(dir["stats"], "bbmerge", "{sample}_bbmerge.out")
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
        ihist={output.hist} \
        {params} \
        &> {log}
        """

rule index_host:
    input:
        config["host_ref_masked"]
    output:
        config["host_ref_indexed"]
    threads: 12
    conda:
        os.path.join(dir["env"], "minimap.yaml") 
    shell:
        """
        minimap2 -d {output} {input}
        """


rule host_removal:
    input:
        merged = os.path.join(dir["output"], "bbmerge", "{sample}_merged.fastq.gz"),
        unmerged1 = os.path.join(dir["output"], "bbmerge", "{sample}_R1_unmerged.fastq.gz"),
        unmerged2 = os.path.join(dir["output"], "bbmerge", "{sample}_R2_unmerged.fastq.gz"),
        index = config["host_ref_indexed"]
    output:
        merged_bam = temp(os.path.join(dir["output"], "host_removed", "{sample}_merged.bam")),
        merged_hr = os.path.join(dir["output"], "host_removed", "{sample}_merged_hr.fastq.gz"),
        unmerged_bam = temp(os.path.join(dir["output"], "host_removed", "{sample}_unmerged.bam")),
        unmerged_bam_sorted = temp(os.path.join(dir["output"], "host_removed", "{sample}_unmerged_sorted.bam")),
        unmerged1_hr = os.path.join(dir["output"], "host_removed", "{sample}_unmerged_hr_R1.fastq.gz"),
        unmerged2_hr = os.path.join(dir["output"], "host_removed", "{sample}_unmerged_hr_R2.fastq.gz"),
    threads: 24
    conda:
        os.path.join(dir["env"], "minimap.yaml") 
    log:
        os.path.join(dir["logs"], "host_removal", "{sample}_hr.log")
    benchmark:
        os.path.join(dir["bench"], "host_removal", "{sample}_hr.txt")
    shell:
        """
        # Merged reads
        minimap2 -ax sr {input.index} {input.merged} \
            | samtools view -bh -f 4 -F 256 \
            | samtools sort -o {output.merged_bam}
        samtools index {output.merged_bam}

        samtools fastq {output.merged_bam} | gzip -c > {output.merged_hr} 

        # Unmerged pairs
        minimap2 -ax sr {input.index} {input.unmerged1} {input.unmerged2} \
            | samtools view -bh \
            | samtools sort -o {output.unmerged_bam}
        samtools index {output.unmerged_bam} 

        samtools view -u -f 12 -F 256 {output.unmerged_bam} \
            | samtools sort -n -o {output.unmerged_bam_sorted}

        samtools fastq -1 {output.unmerged1_hr} -2 {output.unmerged2_hr} \
        -0 /dev/null -s /dev/null \
        -n {output.unmerged_bam_sorted}      

        """


rule primer_b_summary:
    input:
        step1 = expand(os.path.join(dir["stats"], "qc", "step_1", "{sample}_s1.stats"), sample = SAMPLES),
        step2 = expand(os.path.join(dir["stats"], "qc", "step_2", "{sample}_s2.stats"), sample = SAMPLES)
    output:
        os.path.join(dir["stats"], "qc", "primer_b_stats.tsv")
    conda:
        os.path.join(dir["env"], "pandas.yaml")
    script:
        "../scripts/primer_b_stats.py"

rule vector_stats_summary:
    input:
        expand(os.path.join(dir["stats"], "qc", "rm_vector_contamination", "{sample}_rm_vc.stats"), sample = SAMPLES)
    output:
        os.path.join(dir["stats"], "qc", "vector_stats.tsv")
    conda:
        os.path.join(dir["env"], "pandas.yaml")
    script:
        "../scripts/vector_stats.py"

rule bbmerge_summary:
    input:
        expand(os.path.join(dir["stats"], "bbmerge", "{sample}_bbmerge.out"), sample = SAMPLES)
    output:
        stats = os.path.join(dir["stats"], "qc", "bbmerge_stats.tsv"),
        hist = os.path.join(dir["stats"], "qc", "bbmerge_insert_hist.tsv")
    conda:
        os.path.join(dir["env"], "pandas.yaml")
    script:
        "../scripts/bbmerge_stats.py"

rule host_removal_summary:
    input:
        merged_bam = expand(os.path.join(dir["output"], "host_removed", "{sample}_merged.bam"), sample = SAMPLES),
        unmerged_bam = expand(os.path.join(dir["output"], "host_removed", "{sample}_unmerged.bam"), sample = SAMPLES)
    output:
        os.path.join(dir["stats"], "qc", "host_removal_stats.tsv")
    conda:
        os.path.join(dir["env"], "pandas.yaml")
    script:
        "../scripts/host_removal_stats.py"

rule preprocessing_plots:
    input:
        fastp = os.path.join(dir["stats"], "qc", "fastp_stats.tsv"),
        primer_b = os.path.join(dir["stats"], "qc", "primer_b_stats.tsv"),
        vector = os.path.join(dir["stats"], "qc", "vector_stats.tsv"),
        bbmerge = os.path.join(dir["stats"], "qc", "bbmerge_stats.tsv"),
        bbmerge_hist = os.path.join(dir["stats"], "qc", "bbmerge_insert_hist.tsv"),
        host = os.path.join(dir["stats"], "qc", "host_removal_stats.tsv")
    output:
        os.path.join(dir["stats"], "qc", "preprocessing_plots.html"),
    conda:
        os.path.join(dir["env"], "R_env.yaml")
    script:
        "../scripts/preprocessing_plots.Rmd"
