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
# sra_accessions = 'SRR8615581'
sra_accessions = sra_accessions[0:5]
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
print(sra_metadata.columns)
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
        sratoolkit_docker
    script:
        "scripts/convertSRAtoFASTQ.sh"        
        # "hi.txt"

def get_disk_mb(wildcards):
    # use wildcards.run to get the 'size_in_gb' from the sra_metadata
    size_in_gb = sra_metadata[sra_metadata['run_accession'] == wildcards.run]['size_in_GB'].values[0]

    # return the size in MB and add 10% for overhead
    return int((size_in_gb * 1000)*1.1)


rule download_sra:
    output:
        sra_zipped="{PATH_TO_DATA}/SRA/{run}.tar.gz",
        filelist="{PATH_TO_DATA}/SRA/{run}_files.txt"
    container:
        sratoolkit_docker
    threads:
        threads = 1
    resources:
        disk_mb = get_disk_mb
    script:
        "scripts/getSRA.sh"

### THIS SECTION IS REGARDING THE REFERENCE FILES

# input function for the rule aggregate to make sure all required cache files exist
def get_all_cache_files(wildcards):
    # decision based on content of output file
    # Important: use the method open() of the returned file!
    # This way, Snakemake is able to automatically download the file if it is generated in
    # a cloud environment without a shared filesystem.
    with checkpoints.create_all_refseqs.get().output[0].open() as f:
        # read in every line of the file and return it as a list
        return ["cachefiles/" + line.strip() for line in f]

checkpoint create_all_refseqs: 
    input:
        reflist=expand("{PATH_TO_DATA}/cachefiles/{run}_refs.csv", PATH_TO_DATA=PATH_TO_DATA, run=sra_accessions)
    output:
        "rawdata/RNA/refseqlists/all_refseqs.csv"
    resources:
        machine_type = med_cpu
    retries: 2
    script:
        "scripts/create_all_refseqs.py"

# the output here is to remove ambiguity
rule create_ref_seq_list:
    output:
       run_refs="{PATH_TO_DATA}/cachefiles/{run}_refs.csv"
    container:
        sratoolkit_docker
    retries: 5
    threads:
        1
    shell:
        "/opt/sratoolkit.3.0.7-ubuntu64/bin/align-info --ref ${wildcards.run} > ${output}"




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