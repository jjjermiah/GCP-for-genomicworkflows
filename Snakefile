from os.path import join, basename, dirname, splitext
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
julia_list = [
    "SRR8615281", "SRR8615282", "SRR8615300", "SRR8615338", "SRR8615368", "SRR8615369", "SRR8615376",
    "SRR8615398", "SRR8615459", "SRR8615475", "SRR8615504", "SRR8615506", "SRR8615561", "SRR8615564",
    "SRR8615580", "SRR8615581", "SRR8615671", "SRR8615709", "SRR8615749", "SRR8615753", "SRR8615788",
    "SRR8615803", "SRR8615817", "SRR8615854", "SRR8615856", "SRR8615896", "SRR8615924", "SRR8615933",
    "SRR8615934", "SRR8616011", "SRR8616014", "SRR8616032", "SRR8616033", "SRR8616044", "SRR8616058",
    "SRR8616142", "SRR8616177", "SRR8618304", "SRR8618306", "SRR8618307"
]
# sample_accessions = sample_accessions[1]
sample_accessions = julia_list
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
ref_path = f"reference_genomes/ENSEMBL/{reference_genome['SPECIES']}/{reference_genome['BUILD']}/release-{reference_genome['RELEASE']}/"

rule all: 
    input:
        expand("results/{PROJECT_NAME}/CIRI2/{sample}.tsv", sample=sample_accessions, PROJECT_NAME= "CCLE"),
        # expand("procdata/{PROJECT_NAME}/star/pe/{sample}/{sample}_pe_aligned.sam", sample=sample_accessions, PROJECT_NAME= "CCLE"),
        # expand("results/{PROJECT_NAME}/CIRI2/{sample}.tsv", sample="586986_1", PROJECT_NAME= "gCSI"),
        # expand("{sample}_{split}_fastqc.done", sample=sample_accessions, split=[1,2]),
        # expand(join("processed_data/{PROJECT_NAME}/", "alignment/{sample}.sam"), sample=sample_accessions)

include: "workflow/rules/reference_genome.smk"
include: "workflow/rules/sra_fastq.smk"
include: "workflow/rules/fastqc.smk"
include: "workflow/rules/circularRNA.smk"
include: "workflow/rules/bwa.smk"
include: "workflow/rules/star.smk"
include: "workflow/rules/bowtie.smk"