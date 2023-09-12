from os.path import join, basename,dirname
import re
import pandas as pd

configfile: "config.yaml"
machines = config['machine_type']
# snakemake can implicitly use the --defaultremoteprovider and --defaultremoteprefix from command line for a given bucket
# but we use this for the data that might not be in that specific bucket
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()
GS = GSRemoteProvider()

# snakemake --kubernetes --default-remote-provider GS     --default-remote-prefix orcestra-archive/     --use-singularity --keep-going  -j2 

# TODO:: think of a better naming for my dockers.
sratoolkit_docker= "docker://jjjermiah/sratoolkit:0.2"
pigz_docker = "docker://jjjermiah/pigz:0.9"


# PROJECT_NAME="gCSI"
reference_genome = config["ref"]
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


# Fastq file naming is not uniform between CCLE, gCSI
def get_fastq_pe(wildcards):
    if wildcards.PROJECT_NAME == "CCLE":
        return [join("rawdata/{PROJECT_NAME}/", "FASTQ/{sample}_1.fastq.gz"), join("rawdata/{PROJECT_NAME}/", "FASTQ/{sample}_2.fastq.gz")]
    elif wildcards.PROJECT_NAME == "gCSI":
        return [join("rawdata/{PROJECT_NAME}/", "FASTQ/{sample}_1.rnaseq.fastq.gz"), join("rawdata/{PROJECT_NAME}/", "FASTQ/{sample}_2.rnaseq.fastq.gz")]

# Path to store the reference genome data 
ref_path = f"reference_genomes/{reference_genome['SPECIES']}/{reference_genome['RELEASE']}/{reference_genome['BUILD']}/"

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




include: "workflow/rules/reference_genome.smk"
include: "workflow/rules/get_SRA_FASTQ.smk"
include: "workflow/rules/fastqc.smk"
include: "workflow/rules/circularRNA.smk"
include: "workflow/rules/bwa.smk"