
rule bwa_mem:
    input:
        reads=get_fastq_pe,
        idx_genome=join(ref_path, "bwa", "genome.fa"),
        idx=multiext(join(ref_path, "bwa", "genome"), ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
        output=join("procdata/{PROJECT_NAME}/", "alignment/{sample}.sam")
    params:
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'"
    conda:
        "../envs/bwa_mem.yaml"
    threads: 16
    resources:
        machine_type = machines['high_mem']['name']
    shell:
        "ls -la $(dirname {input.idx_genome}); bwa mem -t {threads} $(dirname {input.idx_genome})/genome {input.reads} > {output}"
        
rule build_bwa_index:
    input:
        ref=join(ref_path, "bwa", "{genome}.fa"),
    output:
        idx=multiext(join(ref_path, "bwa", "{genome}"), ".amb", ".ann", ".bwt", ".pac", ".sa"),
    params:
        algorithm="bwtsw",
    log:
        "logs/bwa_index/{genome}.log"
    wrapper:
        "v2.6.0/bio/bwa/index"
