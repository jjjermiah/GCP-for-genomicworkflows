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
# from google.cloud import storage

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

PROJECT_NAME="CCLE"

REF_SPECIES = config["ref"]["SPECIES"]
REF_DATATYPE = config["ref"]["REF_DATATYPE"]
REF_BUILD = config["ref"]["REF_BUILD"]
REF_RELEASE = config["ref"]["REF_RELEASE"]

# end each path with / 
ref_path = f"reference_genomes/{REF_SPECIES}/{REF_RELEASE}/{REF_BUILD}/"
rawdata_path = f"rawdata/{PROJECT_NAME}/"
procdata_path = f"processed_data/{PROJECT_NAME}/"


######### 
# Project specific data

# # path to STORE SRA folders 
SRA_METADATA_FILE = "metadata/sra_metadata.csv"

# # Get the SRA metadata
sra_metadata = pd.read_csv(SRA_METADATA_FILE)

# # Get the SRA accessions
sra_accessions = sra_metadata["run_accession"].tolist()
# sra_accessions = 'SRR8615581'
sra_accessions = sra_accessions[1]


include: "workflow/rules/ref_data.smk"

rule all: 
    input:
        expand("results/CIRI2/{sample}.tsv", sample=sra_accessions),
        # expand("{sample}_{split}_fastqc.done", sample=sra_accessions, split=[1,2]),
        # expand(join(procdata_path, "alignment/{sample}.sam"), PATH_TO_DATA=PATH_TO_DATA, sample=sra_accessions)
        # expand(join(rawdata_path, "FASTQ/{sample}_{split}.fastq.gz", PATH_TO_DATA=PATH_TO_DATA, sample=sra_accessions, split=[1,2]),
        # expand("cachefiles/all_refseqs.csv", PATH_TO_DATA=PATH_TO_DATA),

rule bowtie2_build:
    input:
        ref=join(ref_path, "genome.fa")
    output:
        multiext(join(ref_path, "genome.fa"), "", ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
    conda:
        "envs/bowtie2.yaml"
    threads: 8
    shell:
        "v2.6.0/bio/bowtie2/build"

rule CIRI2:
    input:
        sam=join(procdata_path, "alignment/{sample}.sam"),
        gtf=join(ref_path, "annotation.gtf"),
        idx=multiext(join(ref_path, "genome.fa"), "", ".amb", ".ann", ".bwt", ".pac", ".sa")
    output:
        CIRI2= "results/CIRI2/{sample}.tsv"
    threads: 20
    log: "log/CIRI2/{sample}.log"
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

rule bwa_mem:
    input:
        reads=[
            join(rawdata_path, "FASTQ/{sample}_1.fastq.gz"), 
            join(rawdata_path, "FASTQ/{sample}_2.fastq.gz")
            ],
        idx=multiext(join(ref_path, "genome.fa"), "", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
        output=join(procdata_path, "alignment/{sample}.sam")
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
#         fq=join(rawdata_path, "FASTQ/{sample}_{split}.fastq.gz")
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