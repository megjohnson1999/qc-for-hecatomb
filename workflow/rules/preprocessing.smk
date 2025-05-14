rule fastp:
    """Run fastp to remove adapter sequences"""
    input:
        r1 = get_read1_path,
        r2 = get_read2_path,
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
        r1 = os.path.join(dir["output"], "host_removed", "{sample}_hr_R1.fastq"),
        r2 = os.path.join(dir["output"], "host_removed", "{sample}_hr_R2.fastq"),
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
        # Ensure input files exist and have content
        [ -s {input.r1} ] || {{ echo "Error: Input R1 file missing or empty" >> {log}; exit 1; }}
        [ -s {input.r2} ] || {{ echo "Error: Input R2 file missing or empty" >> {log}; exit 1; }}
        
        # Run bbmerge (with compression for output)
        bbmerge.sh \
        in1={input.r1} \
        in2={input.r2} \
        out={output.merged} \
        outu1={output.unmerged1} \
        outu2={output.unmerged2} \
        ihist={output.hist} \
        pigz=t \
        {params} \
        &> {log}
        
        # Verify that output files were created successfully
        gzip -t {output.merged} || {{ echo "Error: Merged output file corrupted" >> {log}; exit 1; }}
        gzip -t {output.unmerged1} || {{ echo "Error: Unmerged R1 output file corrupted" >> {log}; exit 1; }}
        gzip -t {output.unmerged2} || {{ echo "Error: Unmerged R2 output file corrupted" >> {log}; exit 1; }}
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
        r1 = os.path.join(dir["output"], "qc", "rm_vector_contamination", "{sample}_R1_rm_vc.fastq.gz"),
        r2 = os.path.join(dir["output"], "qc", "rm_vector_contamination", "{sample}_R2_rm_vc.fastq.gz"),
        index = config["host_ref_indexed"]
    output:
        bam = temp(os.path.join(dir["output"], "host_removed", "{sample}.bam")),
        sorted = temp(os.path.join(dir["output"], "host_removed", "{sample}_sorted.bam")),
        filtered = temp(os.path.join(dir["output"], "host_removed", "{sample}_filtered.bam")),
        r1_hr = os.path.join(dir["output"], "host_removed", "{sample}_hr_R1.fastq"),
        r2_hr = os.path.join(dir["output"], "host_removed", "{sample}_hr_R2.fastq"),
        temp_r1 = temp(os.path.join(dir["output"], "host_removed", "{sample}_temp_R1.fastq")),
        temp_r2 = temp(os.path.join(dir["output"], "host_removed", "{sample}_temp_R2.fastq")),
        stats = os.path.join(dir["logs"], "host_removal", "{sample}_stats.txt")
    threads: 24
    conda:
        os.path.join(dir["env"], "minimap.yaml") 
    log:
        os.path.join(dir["logs"], "host_removal", "{sample}_hr.log")
    benchmark:
        os.path.join(dir["bench"], "host_removal", "{sample}_hr.txt")
    shell:
        """
        # Create output directories if they don't exist
        mkdir -p $(dirname {output.bam})
        mkdir -p $(dirname {output.stats})
        mkdir -p $(dirname {log})
        
        # Check input file integrity thoroughly
        echo "Verifying input file integrity..." >> {log}
        gzip -t {input.r1} || {{ echo "Error: Input R1 file corrupted" >> {log}; exit 1; }}
        gzip -t {input.r2} || {{ echo "Error: Input R2 file corrupted" >> {log}; exit 1; }}
        echo "Input files verified successfully" >> {log}
        
        # Align paired reads to host genome - pipe to samtools to ensure proper BAM creation
        echo "Aligning reads to host genome..." >> {log}
        minimap2 -ax sr -t {threads} {input.index} {input.r1} {input.r2} 2>> {log} | samtools view -bh > {output.bam}
        
        # Check if output BAM file is valid
        if ! samtools quickcheck {output.bam}; then
            echo "Error: Alignment produced invalid BAM file" >> {log}
            exit 1
        fi
        
        # Print initial alignment stats
        echo "Initial alignment statistics:" > {output.stats}
        samtools flagstat {output.bam} >> {output.stats}
        
        # Filter reads where both in pair are unmapped (-f 12)
        # Also exclude secondary alignments (-F 256)
        echo "Filtering unmapped reads..." >> {log}
        samtools view -bh -f 12 -F 256 {output.bam} > {output.filtered}
        
        # Sort by name for paired extraction
        echo "Sorting by read name..." >> {log}
        samtools sort -n {output.filtered} > {output.sorted}
        
        # Print filtered alignment stats
        echo -e "\\nFiltered alignment statistics:" >> {output.stats}
        samtools flagstat {output.sorted} >> {output.stats}
        
        # Extract paired reads to temporary files first
        echo "Extracting unmapped reads to temporary files..." >> {log}
        samtools fastq -N -1 {output.temp_r1} -2 {output.temp_r2} -0 /dev/null -s /dev/null -n {output.sorted} 2>> {log}
        
        # Handle the case of empty or missing temporary files by creating empty fastq files
        # This can happen if all reads mapped to host (no reads to extract)
        if [ ! -s "{output.temp_r1}" ] || [ ! -s "{output.temp_r2}" ]; then
            echo "Warning: No unmapped reads found (all reads mapped to host)" >> {log}
            touch {output.temp_r1}
            touch {output.temp_r2}
        fi
        
        # Copy temporary files to final output locations
        echo "Copying to final output files..." >> {log}
        cat {output.temp_r1} > {output.r1_hr}
        cat {output.temp_r2} > {output.r2_hr}
        
        # Count reads in the output files
        r1_count=$(wc -l < {output.r1_hr} 2>/dev/null | awk '{{print int($1/4)}}' || echo 0)
        r2_count=$(wc -l < {output.r2_hr} 2>/dev/null | awk '{{print int($1/4)}}' || echo 0)
        
        echo -e "\\nRead counts in output files:" >> {output.stats}
        echo "R1 reads: $r1_count" >> {output.stats}
        echo "R2 reads: $r2_count" >> {output.stats}
        
        # Ensure read counts match, but don't fail if both are zero
        if [ "$r1_count" -ne "$r2_count" ]; then
            echo "ERROR: Unequal read counts in host-removed files!" >> {log}
            exit 1
        elif [ "$r1_count" -eq 0 ] && [ "$r2_count" -eq 0 ]; then
            echo "Warning: No reads remained after host removal for sample {wildcards.sample}" >> {log}
        fi
        
        echo "Host removal completed successfully for sample {wildcards.sample}" >> {log}
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
        r1 = expand(os.path.join(dir["output"], "host_removed", "{sample}_hr_R1.fastq"), sample=SAMPLES),
        r2 = expand(os.path.join(dir["output"], "host_removed", "{sample}_hr_R2.fastq"), sample=SAMPLES)
    output:
        os.path.join(dir["stats"], "qc", "host_removal_stats.tsv")
    conda:
        os.path.join(dir["env"], "pandas.yaml") # Ensure this points to the correct environment file
    script:
        "../scripts/host_removal_stats.py"

rule preprocessing_plots:
    input:
        raw_stats = os.path.join(dir["stats"], "raw_input_data", "basic_stats.txt"),
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
