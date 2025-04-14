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
        r2 = os.path.join(dir["output"], "assembly", "all_unmerged_R2.fastq.gz"),
        verification = os.path.join(dir["logs"], "assembly", "read_count_verification.txt")
    output:
        contigs = os.path.join(dir["output"], "assembly", "megahit", "final.contigs.fa"),
        dir = directory(os.path.join(dir["output"], "assembly", "megahit"))
    params:
        min_contig = 1000,
        out_dir = os.path.join(dir["output"], "assembly", "megahit"),
        temp_r1 = temp(os.path.join(dir["temp"], "balanced_R1.fastq.gz")),
        temp_r2 = temp(os.path.join(dir["temp"], "balanced_R2.fastq.gz"))
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
        
        # Check if read counts are equal
        r1_count=$(zcat {input.r1} | wc -l | awk '{{print $1/4}}')
        r2_count=$(zcat {input.r2} | wc -l | awk '{{print $1/4}}')
        
        if [ "$r1_count" -eq "$r2_count" ]; then
            echo "Read counts are equal, proceeding with assembly" >> {log}
            
            # Run megahit with original files
            megahit --12 {input.merged} -1 {input.r1} -2 {input.r2} \
                -o {params.out_dir} \
                --min-contig-len {params.min_contig} \
                --k-list 21,29,39,59,79,99,119,141 \
                -t {threads} \
                &>> {log}
        else
            echo "WARNING: Unequal read counts detected. Attempting to balance files..." >> {log}
            echo "R1 reads: $r1_count, R2 reads: $r2_count" >> {log}
            
            # Determine which file has more reads and needs resampling
            if [ "$r1_count" -gt "$r2_count" ]; then
                echo "R1 has more reads. Subsampling to match R2..." >> {log}
                
                # Use reformat.sh to subsample the larger file to match the smaller one
                reformat.sh in={input.r1} out={params.temp_r1} samplereadstarget=$r2_count sampleseed=42 &>> {log}
                cp {input.r2} {params.temp_r2}
            else
                echo "R2 has more reads. Subsampling to match R1..." >> {log}
                
                # Use reformat.sh to subsample the larger file to match the smaller one
                cp {input.r1} {params.temp_r1}
                reformat.sh in={input.r2} out={params.temp_r2} samplereadstarget=$r1_count sampleseed=42 &>> {log}
            fi
            
            # Verify balanced files
            balanced_r1=$(zcat {params.temp_r1} | wc -l | awk '{{print $1/4}}')
            balanced_r2=$(zcat {params.temp_r2} | wc -l | awk '{{print $1/4}}')
            echo "After balancing: R1 reads: $balanced_r1, R2 reads: $balanced_r2" >> {log}
            
            # Run megahit with balanced files
            megahit --12 {input.merged} -1 {params.temp_r1} -2 {params.temp_r2} \
                -o {params.out_dir} \
                --min-contig-len {params.min_contig} \
                --k-list 21,29,39,59,79,99,119,141 \
                -t {threads} \
                &>> {log}
        fi
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
