# CIRI CIRCExplorer

rule CIRI2:
    input:
        sam=join("processed_data/{PROJECT_NAME}/", "alignment/{sample}.sam"),
        gtf=join(ref_path, "annotation.gtf"),
        idx=multiext(join(ref_path, "genome.fa"), "", ".amb", ".ann", ".bwt", ".pac", ".sa")
    output:
        CIRI2= "results/{PROJECT_NAME}/CIRI2/{sample}.tsv"
    threads: 20
    log: "log/{PROJECT_NAME}/CIRI2/{sample}.log"
    container:
        "docker://andremrsantos/ciri2:latest"
    resources:
        machine_type = machines['small_cpu']['name']
    shell:
        """CIRI2 \
        --thread_num {threads} \
        --in {input.sam} \
        --anno {input.gtf} \
        --out {output.CIRI2} \
        --ref_file {input.idx[0]} \
        --log {log}"""


# rule: 
#     input:

#     output:

#     shell:
#         """
#         STAR \
#         --runThreadN {threads} \
#         --genomeDir {params.index} \
#         --readFilesIn {input.reads1} {input.reads2} \
#         --readFilesCommand zcat \
#         --outFileNamePrefix {params.outFileNamePrefix} \
#         --outSAMattributes All \
#         --outStd BAM_Unsorted \
#         --outSAMtype BAM Unsorted \
#         --outSAMattrRGline ID:rnaseq_pipeline SM:{params.sample_id} \
#         {params.additional_params} \
#         > {output.bam};
#         """

# rule pe_map_genome_star:
#     """
#         Map to genome using STAR
#     """
#     input:
#         index=lambda wildcards: os.path.join(
#             config["star_indexes"],
#             get_sample("organism", search_id="index", search_value=wildcards.sample),
#             get_sample("index_size", search_id="index", search_value=wildcards.sample),
#             "STAR_index",
#             "chrNameLength.txt",
#         ),
#         reads1=os.path.join(
#             config["output_dir"],
#             "samples",
#             "{sample}",
#             "{sample}.fq1.pe.remove_polya.fastq.gz",
#         ),
#         reads2=os.path.join(
#             config["output_dir"],
#             "samples",
#             "{sample}",
#             "{sample}.fq2.pe.remove_polya.fastq.gz",
#         ),
#     output:
#         bam=os.path.join(
#             config["output_dir"],
#             "samples",
#             "{sample}",
#             "map_genome",
#             "{sample}.pe.Aligned.out.bam",
#         ),
#         logfile=os.path.join(
#             config["output_dir"],
#             "samples",
#             "{sample}",
#             "map_genome",
#             "{sample}.pe.Log.final.out",
#         ),
#     params:
#         cluster_log_path=config["cluster_log_dir"],
#         sample_id="{sample}",
#         index=lambda wildcards: os.path.abspath(
#             os.path.join(
#                 config["star_indexes"],
#                 get_sample(
#                     "organism", search_id="index", search_value=wildcards.sample
#                 ),
#                 get_sample(
#                     "index_size", search_id="index", search_value=wildcards.sample
#                 ),
#                 "STAR_index",
#             )
#         ),
#         outFileNamePrefix=lambda wildcards, output: output.bam.replace(
#             "Aligned.out.bam", ""
#         ),
#         additional_params=parse_rule_config(
#             rule_config,
#             current_rule=current_rule,
#             immutable=(
#                 "--genomeDir",
#                 "--readFilesIn",
#                 "--readFilesCommand",
#                 "--outFileNamePrefix",
#                 "--outSAMattributes",
#                 "--outStd",
#                 "--outSAMtype",
#                 "--outSAMattrRGline",
#             ),
#         ),
#     container:
#         "docker://quay.io/biocontainers/star:2.7.8a--h9ee0642_1"
#     conda:
#         os.path.join(workflow.basedir, "envs", "STAR.yaml")
#     threads: 12
#     resources:
#         mem_mb=lambda wildcards, attempt: 32000 * attempt,
#     log:
#         stderr=os.path.join(
#             config["log_dir"], "samples", "{sample}", current_rule + ".stderr.log"
#         ),
#     shell:
#         "(STAR \
#         --runThreadN {threads} \
#         --genomeDir {params.index} \
#         --readFilesIn {input.reads1} {input.reads2} \
#         --readFilesCommand zcat \
#         --outFileNamePrefix {params.outFileNamePrefix} \
#         --outSAMattributes All \
#         --outStd BAM_Unsorted \
#         --outSAMtype BAM Unsorted \
#         --outSAMattrRGline ID:rnaseq_pipeline SM:{params.sample_id} \
#         {params.additional_params} \
#         > {output.bam};) \
#         2> {log.stderr}"