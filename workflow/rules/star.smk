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


# rule create_index_star:
#     """
#         Create index for STAR alignments
#     """
#     input:
#         genome=join(ref_path, "genome.fa"),
#         gtf=join(ref_path, "annotation.gtf"),
#     output:
#         chromosome_info = join(ref_path, "STAR_INDEX", "chrNameLength.txt"),
#         chromosomes_names = join(ref_path, "STAR_INDEX", "chrName.txt"),
#     container:
#         "docker://quay.io/biocontainers/star:2.7.8a--h9ee0642_1"