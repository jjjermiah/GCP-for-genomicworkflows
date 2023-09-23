rule kallisto_index:
    input:
        fasta=f"{ref_path}/transcriptome.fa",
    output:
        index=f"{ref_path}/transcriptome.idx",
    params:
        extra="",  # optional parameters
    log:
        "logs/kallisto_index_transcriptome.log",
    threads: 32
    conda:
        "../envs/kallisto.yaml"
    shell:
        """
        kallisto index \
            --threads {snakemake.threads} \
            -i {output.index} \
            {params.extra} \
            {input.fasta} \
            2> {log}
        """

rule kallisto_quant:
    input:
        fastq=get_fastq_pe,
        index=f"{ref_path}/transcriptome.idx",
    output:
        abundancefile="results/{PROJECT_NAME}/kallisto_quant/{ref_build}/{sample}_abundance.tsv",
    params:
        extra="",
    log:
        log="logs/{PROJECT_NAME}/{ref_build}_{sample}_kallisto_quant.log",
    threads: 15
    conda:
        "../envs/kallisto.yaml"
    shell:
        """
        kallisto quant \
            -i {input.index} \
            -o {output} \
            {params.extra} \
            {input.fastq} \
            2> {log}
        """