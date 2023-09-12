rule bwa_mem:
    input:
        get_fastq_pe,
        idx=multiext(join(ref_path, "genome.fa"), "", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
        output=join("processed_data/{PROJECT_NAME}/", "alignment/{sample}.sam")
    params:
        extra=r"-R '@RG\tID:{wildcards.sample}\tSM:{wildcards.sample}'"
    conda:
        "envs/bwa_mem.yaml"
    threads: 32
    resources:
        machine_type = machines['high_mem']['name']
    shell:
        "bwa mem -t {threads} {params.extra} {input.idx[0]} {input.reads} > {output}"
        
rule build_bwa_index:
    input:
        ref=join(ref_path, "genome.fa"),
    output:
        idx=multiext(join(ref_path, "genome.fa"), "", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    params:
        algorithm="bwtsw",
    wrapper:
        "v2.6.0/bio/bwa/index"
