rule kallisto_index:
    input:
        fasta="{transcriptome}.fasta",
    output:
        index="{transcriptome}.idx",
    params:
        extra="",  # optional parameters
    log:
        "logs/kallisto_index_{transcriptome}.log",
    threads: 1
    wrapper:
        "v2.6.0/bio/kallisto/index"