rule star_index:
    input:
        fasta=join(ref_path, "genome.fa"),
        gtf=join(ref_path, "annotation.gtf"),
    output:
        directory(join(ref_path, "{genome}", "STAR_INDEX"))
    threads: 1
    conda:
        "../envs/star.yaml"
    params:
        extra="",
    log:
        "logs/star_index_{genome}.log",
    shell:
        """
        STAR \
        --runThreadN {threads} \
        --runMode genomeGenerate \
        --genomeFastaFiles {input.fasta} \
        --sjdbOverhang 99 \
        --sjdbGTFfile {input.gtf} \
        --genomeDir {output}
        """

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




#Generate Star genome index
## genomeSAindexNbases calculated to 9,557860944 (...min(14, log2(GenomeLength)/2 - 1))
# STAR --runThreadN 4 --runMode genomeGenerate --genomeSAindexNbases 9 --genomeDir /proj/uppstore2017134/stress/circ/star-genome9 \
# --genomeFastaFiles /*.fa --sjdbGTFfile NC_003112.gff

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
#     params:
#         cluster_log_path=config["cluster_log_dir"],
#         output_dir=lambda wildcards, output: os.path.dirname(output.chromosome_info),
#         outFileNamePrefix=os.path.join(
#             config["star_indexes"], "{organism}", "{index_size}", "STAR_index/STAR_"
#         ),
#         sjdbOverhang="{index_size}",
#         additional_params=parse_rule_config(
#             rule_config,
#             current_rule=current_rule,
#             immutable=(
#                 "--runMode",
#                 "--sjdbOverhang",
#                 "--genomeDir",
#                 "--genomeFastaFiles",
#                 "--outFileNamePrefix",
#                 "--sjdbGTFfile",
#             ),
#         ),
#     container:
#         "docker://quay.io/biocontainers/star:2.7.8a--h9ee0642_1"
#     conda:
#         "../envs/STAR.yaml"
#     threads: 8
#     resources:
#         mem_mb=lambda wildcards, attempt: 32000 * attempt,
#     shell:"""
#         STAR \
#         --runMode genomeGenerate \
#         --sjdbOverhang {params.sjdbOverhang} \
#         --genomeDir {params.output_dir} \
#         --genomeFastaFiles {input.genome} \
#         --runThreadN {threads} \
#         --outFileNamePrefix {params.outFileNamePrefix} \
#         --sjdbGTFfile {input.gtf}) \
#         {params.additional_params} \
#         1> {log.stdout} 2> {log.stderr}
#         """