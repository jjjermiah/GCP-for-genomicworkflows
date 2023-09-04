from os.path import join, basename,dirname
import re
import pandas as pd
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
# from google.cloud import storage

# snakemake --kubernetes  -j2 --default-remote-provider GS     --default-remote-prefix orcestra-archive/     --use-singularity --keep-going 

# TODO:: think of a better naming for my dockers. sratools:0.2 is specifically for fasterq-dump
sratoolkit_docker= "docker://jjjermiah/sratoolkit:0.2"
sra_docker = "docker://jjjermiah/sratools:0.2"
pigz_docker = "docker://jjjermiah/pigz:0.9"

# # path to STORE SRA folders 
SRA_METADATA_FILE = "metadata/sra_metadata.csv"

PATH_TO_DATA = "rawdata/temp"

# # Get the SRA metadata
sra_metadata = pd.read_csv(SRA_METADATA_FILE)

# # Get the SRA accessions
sra_accessions = sra_metadata["run_accession"].tolist()
sra_accessions = 'SRR8615581'

#### CONFIGURE RESOURCES
# resources:
small_cpu = "e2-standard-8"
med_cpu = "e2-standard-8"
high_cpu = "e2-standard-8"
high_mem = "n1-highmem-32"
# Make dictionaries of machines
machine_dict = {
    "small_cpu": small_cpu,
    "med_cpu": med_cpu,
    "high_cpu": high_cpu,
    "high_mem": high_mem
}

rule all: 
    input:
        expand("{PATH_TO_DATA}/FASTQ/{run}_{split}.fastq.gz", PATH_TO_DATA=PATH_TO_DATA, run=sra_accessions, split=[1,2]) 
        
rule convertSRAtoFASTQ:
    input:
        expand("{PATH_TO_DATA}/SRA/{run}.tar.gz", PATH_TO_DATA=PATH_TO_DATA, run=sra_accessions)
    output:
        fq1="{PATH_TO_DATA}/FASTQ/{run}_1.fastq.gz",
        fq2="{PATH_TO_DATA}/FASTQ/{run}_2.fastq.gz"
    resources:
        machine_type = high_cpu
    threads:
        threads = 7
    container:
        sra_docker
    script:
        "scripts/convertSRAtoFASTQ.sh"        
        # "hi.txt"

rule download_sra:
    output:
        sra_zipped="{PATH_TO_DATA}/SRA/{run}.tar.gz",
        reflist="{PATH_TO_DATA}/SRA/{run}_refs.csv",
        filelist="{PATH_TO_DATA}/SRA/{run}_files.txt"
    container:
        sratoolkit_docker
    script:
        "scripts/getSRA.sh"



# rule download_sra:
#     output:
#         srafile="{PATH_TO_DATA}/{run}/{run}.sra",
#         reflist="{PATH_TO_DATA}/{run}/{run}_refs.csv",
#         vdbcachefile="{PATH_TO_DATA}/{run}/{run}.sra.vdbcache"
#     container:
#         "docker://pegi3s/sratoolkit:3.0.7"
#     script:
#         "scripts/getSRA.sh"

# docker pull pegi3s/sratoolkit:3.0.7