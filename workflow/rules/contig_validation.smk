rule index_contigs:
    """Index contigs using minimap2 for read mapping"""
    input:
        os.path.join(dir["output"], "assembly", "megahit", "final.contigs.fa")
    output:
        os.path.join(dir["output"], "assembly", "megahit", "final.contigs.mmi")
    threads: 12
    conda:
        os.path.join(dir["env"], "minimap_env.yaml")
    log:
        os.path.join(dir["logs"], "contig_validation", "index_contigs.log")
    benchmark:
        os.path.join(dir["bench"], "contig_validation", "index_contigs.txt")
    shell:
        """
        minimap2 -d {output} {input} 2> {log}
        """

rule map_reads_to_contigs:
    """Map host-removed reads back to contigs to evaluate assembly quality"""
    input:
        r1 = os.path.join(dir["output"], "host_removed", "{sample}_hr_R1.fastq"),
        r2 = os.path.join(dir["output"], "host_removed", "{sample}_hr_R2.fastq"),
        idx = os.path.join(dir["output"], "assembly", "megahit", "final.contigs.mmi")
    output:
        bam = os.path.join(dir["output"], "contig_validation", "mapping", "{sample}.bam"),
        sorted = os.path.join(dir["output"], "contig_validation", "mapping", "{sample}.sorted.bam"),
        stats = os.path.join(dir["stats"], "contig_validation", "mapping", "{sample}_mapping_stats.txt")
    threads: 24
    conda:
        os.path.join(dir["env"], "minimap_env.yaml")
    log:
        os.path.join(dir["logs"], "contig_validation", "map_reads", "{sample}.log")
    benchmark:
        os.path.join(dir["bench"], "contig_validation", "map_reads", "{sample}.txt")
    shell:
        """
        # Map reads to contigs
        minimap2 -ax sr -t {threads} {input.idx} {input.r1} {input.r2} | \
        samtools view -bh > {output.bam} 2> {log}
        
        # Sort BAM by coordinate
        samtools sort -@ {threads} -o {output.sorted} {output.bam} 2>> {log}
        
        # Index sorted BAM
        samtools index {output.sorted} 2>> {log}
        
        # Basic mapping statistics
        samtools flagstat {output.sorted} > {output.stats} 2>> {log}
        """

rule generate_per_sample_coverage:
    """Generate per-sample coverage in 1kb windows for each contig"""
    input:
        bam = os.path.join(dir["output"], "contig_validation", "mapping", "{sample}.sorted.bam"),
        contigs = os.path.join(dir["output"], "assembly", "megahit", "final.contigs.fa")
    output:
        coverage = os.path.join(dir["stats"], "contig_validation", "coverage", "{sample}_window_coverage.bed")
    params:
        window_size = 1000
    threads: 8
    conda:
        os.path.join(dir["env"], "minimap_env.yaml")
    log:
        os.path.join(dir["logs"], "contig_validation", "window_coverage", "{sample}.log")
    shell:
        """
        set -exo pipefail
    
        # Create output directory
        mkdir -p $(dirname {output.coverage})

        # Define sample-specific temporary files
        temp_seqkit_output="{wildcards.sample}_temp_seqkit_output.tab"
        temp_contigs="{wildcards.sample}_temp_contigs.bed"
        temp_windows="{wildcards.sample}_temp_windows.bed"
    
        # Get contig lengths correctly
        echo "Getting contig lengths" >> {log}
        seqkit fx2tab -nl {input.contigs} > $temp_seqkit_output 2>> {log}

        # Inspect output
        echo "First 10 lines of $temp_seqkit_output:" >> {log}
        head -n 10 $temp_seqkit_output >> {log}
    
        # Extract relevant columns and convert to BED format
        echo "Converting to BED format" >> {log}
        awk '{{print $1 "\\t0\\t" $NF}}' $temp_seqkit_output > $temp_contigs 2>> {log}  
        
        # Inspect output
        echo "First 10 lines of $temp_contigs:" >> {log}
        head -n 10 $temp_contigs >> {log}
    
        # Generate coverage windows
        echo "Generating coverage windows" >> {log}
        bedtools makewindows -b $temp_contigs -w {params.window_size} > $temp_windows 2>> {log}
    
        # Inspect output
        echo "First 10 lines of $temp_windows:" >> {log}
        head -n 10 $temp_windows >> {log}
    
        # Compute coverage
        echo "Computing coverage" >> {log}
        bedtools coverage -a $temp_windows -b {input.bam} -mean > {output.coverage} 2>> {log}
    
        # Clean up temporary files
        echo "Cleaning up" >> {log}
        rm $temp_seqkit_output $temp_contigs $temp_windows
    
        echo "Done" >> {log}
        """

rule detect_chimeric_contigs:
    """Detect chimeric contigs by analyzing per-sample coverage patterns along contigs"""
    input:
        coverages = expand(os.path.join(dir["stats"], "contig_validation", "coverage", "{sample}_window_coverage.bed"), sample=SAMPLES),
        contigs = os.path.join(dir["output"], "assembly", "megahit", "final.contigs.fa")
    output:
        chimeric = os.path.join(dir["stats"], "contig_validation", "chimeric", "chimeric_contigs.txt"),
        heatmap = os.path.join(dir["stats"], "contig_validation", "chimeric", "chimeric_contigs_heatmap.txt")
    params:
        min_contig_length = 3000,  # Only analyze contigs larger than this
        min_coverage = 5,          # Minimum coverage to consider a window covered by a sample
        min_windows = 3,           # Minimum windows in a row to consider a region
        sample_names = lambda wildcards: SAMPLES
    threads: 16
    conda:
        os.path.join(dir["env"], "pandas.yaml")
    log:
        os.path.join(dir["logs"], "contig_validation", "chimeric_detection.log")
    script:
        "../scripts/detect_chimeric_contigs.py"

rule contig_validation_summary:
    """Generate comprehensive summary of contig validation results"""
    input:
        mapping_stats = expand(os.path.join(dir["stats"], "contig_validation", "mapping", "{sample}_mapping_stats.txt"), sample=SAMPLES),
        chimeric = os.path.join(dir["stats"], "contig_validation", "chimeric", "chimeric_contigs.txt")
    output:
        summary = os.path.join(dir["stats"], "contig_validation", "contig_validation_summary.txt")
    threads: 4
    conda:
        os.path.join(dir["env"], "minimap_env.yaml")
    log:
        os.path.join(dir["logs"], "contig_validation", "summary.log")
    shell:
        """
        # Create summary report
        echo "Contig Validation Summary" > {output.summary}
        echo "=========================" >> {output.summary}
        echo "" >> {output.summary}
        
        # Add date
        date >> {output.summary}
        echo "" >> {output.summary}
        
        # Mapping statistics
        echo "Mapping Statistics:" >> {output.summary}
        echo "-------------------" >> {output.summary}
        for stat_file in {input.mapping_stats}; do
            sample=$(basename "$stat_file" _mapping_stats.txt)
            echo "Sample: $sample" >> {output.summary}
            grep "mapped (" "$stat_file" >> {output.summary}
            grep "properly paired" "$stat_file" >> {output.summary}
            echo "" >> {output.summary}
        done
        
        # Chimeric contig summary
        echo "Chimeric Contig Analysis:" >> {output.summary}
        echo "------------------------" >> {output.summary}
        echo "Total potential chimeric contigs detected: $(wc -l < {input.chimeric})" >> {output.summary}
        
        # Add first 10 chimeric contigs as examples
        if [ -s {input.chimeric} ]; then
            echo "" >> {output.summary}
            echo "Examples of chimeric contigs (showing contig ID and contributing samples):" >> {output.summary}
            head -10 {input.chimeric} >> {output.summary}
        fi
        
        echo "Contig validation summary completed at $(date)" >> {log}
        """
