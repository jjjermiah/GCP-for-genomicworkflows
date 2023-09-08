from os.path import join, basename,dirname
import re
import pandas as pd

configfile: "config.yaml"

# snakemake can implicitly use the --defaultremoteprovider and --defaultremoteprefix from command line for a given bucket
# but we use this for the data that might not be in that specific bucket
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()
GS = GSRemoteProvider()

# snakemake --kubernetes --default-remote-provider GS     --default-remote-prefix orcestra-archive/     --use-singularity --keep-going  -j2 

# TODO:: think of a better naming for my dockers. sratools:0.2 is specifically for fasterq-dump
sratoolkit_docker= "docker://jjjermiah/sratoolkit:0.2"
pigz_docker = "docker://jjjermiah/pigz:0.9"

#### CONFIGURE RESOURCES
# resources:`
small_cpu = "e2-standard-8"
med_cpu = "e2-standard-8"
high_cpu = "e2-standard-8"
high_mem = "n1-highmem-32" # 32 CPU, 208 GB RAM
# Make dictionaries of machines
machine_dict = {
    "small_cpu": small_cpu,
    "med_cpu": med_cpu,
    "high_cpu": high_cpu,
    "high_mem": high_mem
}

# PROJECT_NAME="gCSI"

REF_SPECIES = config["ref"]["SPECIES"]
REF_DATATYPE = config["ref"]["REF_DATATYPE"]
REF_BUILD = config["ref"]["REF_BUILD"]
REF_RELEASE = config["ref"]["REF_RELEASE"]

# end each path with / 


######### 
# Project specific data
# PROJECT_NAME = "CCLE"
SRA_METADATA_FILE = "metadata/sra_metadata.csv"
sra_metadata = pd.read_csv(SRA_METADATA_FILE)
sample_accessions = sra_metadata["run_accession"].tolist()
# sample_accessions = 'SRR8615581'
sample_accessions = sample_accessions[1]
# SRA_METADATA_FILE = "metadata/gCSI_metadata.csv"
# gCSI_metadata = pd.read_csv(gCSI_METADATA_FILE)
# sample_accessions = "586986_1"



include: "workflow/rules/ref_data.smk"
include: "workflow/rules/get_SRA_FASTQ.smk"


ref_path = f"reference_genomes/{REF_SPECIES}/{REF_RELEASE}/{REF_BUILD}/"
rule all: 
    input:
        expand("results/{PROJECT_NAME}/CIRI2/{sample}.tsv", sample=sample_accessions, PROJECT_NAME= "CCLE"),
        expand("results/{PROJECT_NAME}/CIRI2/{sample}.tsv", sample="586986_1", PROJECT_NAME= "gCSI"),
        # expand("{sample}_{split}_fastqc.done", sample=sample_accessions, split=[1,2]),
        # expand(join("processed_data/{PROJECT_NAME}/", "alignment/{sample}.sam"), sample=sample_accessions)

rule bowtie2_build:
    input:
        ref=join(ref_path, "genome.fa")
    output:
        multiext(join(ref_path, "genome.fa"), "", ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
    conda:
        "envs/bowtie2.yaml"
    threads: 
        8
    shell:
        "v2.6.0/bio/bowtie2/build"

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
        machine_type = high_mem
    shell:
        """CIRI2 \
        --thread_num {threads} \
        --in {input.sam} \
        --anno {input.gtf} \
        --out {output.CIRI2} \
        --ref_file {input.idx[0]} \
        --log {log}"""

def get_fastq_pe(wildcards):
    if wildcards.PROJECT_NAME == "CCLE":
        return [join("rawdata/{PROJECT_NAME}/", "FASTQ/{sample}_1.fastq.gz"), join("rawdata/{PROJECT_NAME}/", "FASTQ/{sample}_2.fastq.gz")]
    elif wildcards.PROJECT_NAME == "gCSI":
        return [join("rawdata/{PROJECT_NAME}/", "FASTQ/{sample}_1.rnaseq.fastq.gz"), join("rawdata/{PROJECT_NAME}/", "FASTQ/{sample}_2.rnaseq.fastq.gz")]

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
        machine_type = high_mem
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


# rule circexplorer_parse:
#     input: "align/{sample}_star/Chimeric.out.junction"
#     output: "circrna/circexplorer/{sample}_bsj.bed"
#     log: "log/{sample}_ce_parse.log"
#     threads: 4
#     conda: "envs/circexplorer.yaml"
#     shell: """
#     CIRCexplorer2 parse -t STAR -b {output} {input} > {log}
#     """

# rule circexplorer_annotate:
#     input: "circrna/circexplorer/{sample}_bsj.bed"
#     output: "circrna/circexplorer/{sample}_circRNA.tsv"
#     log: "log/{sample}_ce_annotate.log"
#     threads: 4
#     conda: "envs/circexplorer.yaml"
#     params:
#         ref = config["ref"]["bwa"],
#         pred = config["ref"]["pred"]
#     shell: """
#     CIRCexplorer2 annotate \
#     -g {params.ref} -r {params.pred} \
#     -b {input} -o {output} > {log}
#     """


# rule fastqc:
#     input:
#         fq=join("rawdata/{PROJECT_NAME}/", "FASTQ/{sample}_{split}.fastq.gz")
#     output:
#         html="QC/fastqc/{sample}_{split}.html",
#         zip="QC/fastqc/{sample}_{split}fastqc.zip", # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
#         tempfile=touch("{sample}_{split}_fastqc.done")
#     params:
#         extra = "--quiet"
#     threads: 1
#     log:
#         "logs/fastqc/{sample}_{split}.log"
#     resources:
#         mem_mb = 1024
#     wrapper:
#         "v2.6.0/bio/fastqc"


# Note on param extra: the extra  parameter is used to specify the read group information for the output SAM/BAM file. 
# It sets the read group ID (ID) and sample name (SM) to the value of the sample
# This information can be useful for downstream analysis or when working with multiple samples.
# Note on resources: https://hcc.unl.edu/docs/applications/app_specific/bioinformatics_tools/alignment_tools/bwa/samplening_bwa_commands/#useful-information