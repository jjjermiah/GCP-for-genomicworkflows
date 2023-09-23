from os.path import join, basename, dirname, splitext
import re
import pandas as pd

configfile: "config.yaml"
# snakemake can implicitly use the --defaultremoteprovider and --defaultremoteprefix from command line for a given bucket
# but we use this for the data that might not be in that specific bucket
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
GS = GSRemoteProvider()

# snakemake --kubernetes --default-remote-provider GS     --default-remote-prefix orcestra-archive/     --use-singularity --keep-going  -j2 

# TODO:: think of a better naming for my dockers.
sratoolkit_docker= "docker://jjjermiah/sratoolkit:0.2"
pigz_docker = "docker://jjjermiah/pigz:0.9"

# PROJECT_NAME="gCSI"

######### 
# Project specific data
# PROJECT_NAME = "CCLE"
SRA_METADATA_FILE = "metadata/sra_metadata.csv"
sra_metadata = pd.read_csv(SRA_METADATA_FILE)
sample_accessions = sra_metadata["run_accession"].tolist()

# METADATA from JULIA 
SRA_METADATA_FILE = "metadata/metadata_julia_SRA_df.csv"
sra_metadata = pd.read_csv(SRA_METADATA_FILE)
sample_accessions = sra_metadata["Run"].tolist()

#### gCSI 
# Define the path to the .tsv file
tsv_file = "metadata/gCSI_list.tsv"

# Read the .tsv file and extract the second column into a list
file_paths = []
with open(tsv_file, "r") as file:
    for line in file:
        columns = line.strip().split(" ")
        file_paths.append(columns[1])

# Get the basenames of the file paths
basenames = [os.path.basename(path) for path in file_paths]

# Remove everything after the first "_1" in each element of basenames
gCSI_accessions = [name.split("_1")[0] + "_1" for name in basenames]

machines = config["machine_type"]


# Path to store the reference genome data 
reference_genome = config["ref"]
reference_genome_source = "GENCODE"
reference_genome_build = "GRCh38"
gencode_release = 44
# ref_path = f"reference_genomes/{reference_genome_source}/{reference_genome['SPECIES']}/{reference_genome['BUILD']}/release-{reference_genome['RELEASE']}/"
ref_path = f"reference_genomes/{reference_genome_source}/{reference_genome['SPECIES']}/{reference_genome_build}/release-{gencode_release}/"
ref_build = f"{reference_genome_source}_{reference_genome_build}_v{gencode_release}"

# sample_accessions = sample_accessions [:2]
# sample_accessions = ['SRR8615504', 'SRR8615545' ]
rule all: 
    input:
        # expand("procdata/{PROJECT_NAME}/star/pe/{sample}/{sample}_pe_aligned.sam", sample=sample_accessions, PROJECT_NAME= "CCLE"),
        # expand("results/{PROJECT_NAME}/circRNA_finder/{sample}/{sample}_filteredJunctions.bed", PROJECT_NAME="CCLE", sample=sample_accessions),
        # expand("results/{PROJECT_NAME}/circRNA_finder/{sample}/{sample}_filteredJunctions.bed", PROJECT_NAME="gCSI", sample=gCSI_accessions),
        expand("results/{PROJECT_NAME}/kallisto_quant/{ref_build}/{sample}_abundance.tsv", PROJECT_NAME="CCLE", sample=sample_accessions[1], ref_build=ref_build),
        # expand("results/{PROJECT_NAME}/CIRI2/{sample}/{sample}.tsv", PROJECT_NAME= "CCLE", sample=sample_accessions),
        # expand("results/{PROJECT_NAME}/CIRI2/{sample}/{sample}.tsv", PROJECT_NAME= "gCSI", sample=gCSI_accessions),
        # idx=directory(join(ref_path, "star/index")),
        # gencode_annotation_file=f"reference_genomes/GENCODE/homo_sapiens/GRCh37/release-{gencode_release}/annotation.gtf",
        # gencode_genome_file=f"reference_genomes/GENCODE/homo_sapiens/GRCh37/release-{gencode_release}/genome.fa",
        # expand("{sample}_{split}_fastqc.done", sample=sample_accessions, split=[1,2]),


# Fastq file naming is not uniform between CCLE, gCSI
def get_fastq_pe(wildcards):
    if wildcards.PROJECT_NAME == "CCLE":
        return [join("rawdata/{PROJECT_NAME}/", "FASTQ/{sample}_1.fastq.gz"), join("rawdata/{PROJECT_NAME}/", "FASTQ/{sample}_2.fastq.gz")]
    elif wildcards.PROJECT_NAME == "gCSI":
        return [join("rawdata/{PROJECT_NAME}/", "FASTQ/{sample}_1.rnaseq.fastq.gz"), join("rawdata/{PROJECT_NAME}/", "FASTQ/{sample}_2.rnaseq.fastq.gz")]

include: "workflow/rules/reference_genome.smk"
# include: "workflow/rules/sra_fastq.smk"
include: "workflow/rules/kallisto.smk"
include: "workflow/rules/fastqc.smk"
include: "workflow/rules/circularRNA.smk"
include: "workflow/rules/bwa.smk"
include: "workflow/rules/star.smk"
include: "workflow/rules/bowtie.smk"