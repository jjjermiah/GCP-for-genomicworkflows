rule star_index:
    input:
        fasta=join(ref_path, "genome.fa"),
        gtf=join(ref_path, "annotation.gtf"),
    output:
        directory(join(ref_path, "{genome}", "STAR_INDEX"))
    threads: 1
    params:
        extra="",
    log:
        "logs/star_index_{genome}.log",
    wrapper:
        "v2.6.0/bio/star/index"