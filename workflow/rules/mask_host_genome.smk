rule vir_shred:
    input:
        config["viral_db"]
    output:
        temp(os.path.join(dir["temp"], "virus_shreds.fa"))
    params:
        ol = 40,
        med = 80,
        var = 5,
        min = 20
    threads: 24
    conda:
        os.path.join(dir["env"], "bbmap.yaml")
    log:
        os.path.join(dir["logs"], "mask_reference", "shred.log")
    benchmark:
        os.path.join(dir["bench"], "mask_reference", "shred.txt")
    shell:
        """
        shred.sh \
        in={input} \
        out={output} \
        minlength={params.min} \
        overlap={params.ol} \
        median={params.med} \
        variance={params.var}
        """


rule map_shreds:
    input:
        shreds = os.path.join(dir["temp"], "virus_shreds.fa"),
        ref = config["host_ref"]
    output:
        os.path.join(dir["temp"], "mapped_to_host.sam")
    params:
        minid = 90,
        maxindel = 2
    threads: 24
    conda:
        os.path.join(dir["env"], "bbmap.yaml")
    log:
        os.path.join(dir["logs"], "mask_reference", "map.log")
    benchmark:
        os.path.join(dir["bench"], "mask_reference", "map.txt")
    shell:
        """
        bbmap.sh \
        ref={input.ref} \
        in={input.shreds} \
        ambig=all \
        outm={output} \
        minid={params.minid} \
        maxindel={params.maxindel} \
        ow=t \
        threads={threads}
        """


rule mask_host:
    input:
        sam = os.path.join(dir["temp"], "mapped_to_host.sam"),
        ref = config["host_ref"]
    output:
        config["host_ref_masked"]
    params:
        entropy = 0.5
    conda:
        os.path.join(dir["env"], "bbmap.yaml")
    threads: 24
    shell:
        """
        bbmask.sh \
        in={input.ref} \
        out={output} \
        entropy={params.entropy} \
        sam={input.sam} \
        ow=t \
        threads={threads}
        """
