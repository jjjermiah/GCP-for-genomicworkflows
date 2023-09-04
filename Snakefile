from os.path import join, basename,dirname
import re
import pandas as pd
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
# from google.cloud import storage

# snakemake --kubernetes  -j2 --default-remote-provider GS     --default-remote-prefix orcestra-archive/     --use-singularity --keep-going 

sra_docker = "docker://jjjermiah/sratools:0.2"
pigz_docker = "docker://jjjermiah/pigz:0.9"

SRA_METADATA_FILE = "metadata/sra_metadata.csv"

# # path to STORE SRA folders 
PATH_TO_SRA = "rawdata/temp/SRA"

# # Get the SRA metadata
sra_metadata = pd.read_csv(SRA_METADATA_FILE)

# # Get the SRA accessions
sra_accessions = sra_metadata["run_accession"].tolist()

sra_accessions = sra_accessions[2]
print(sra_metadata.iloc[2])

rule all:
    input:
        expand("{PATH_TO_SRA}/{run}/{run}.sra", PATH_TO_SRA=PATH_TO_SRA, run=sra_accessions) 
        # "hi.txt"

rule download_sra:
    output:
        srafile="{PATH_TO_SRA}/{run}/{run}.sra",
        reflist="{PATH_TO_SRA}/{run}/{run}_refs.csv",
        vdbcachefile="{PATH_TO_SRA}/{run}/{run}.sra.vdbcache"
    container:
        "docker://pegi3s/sratoolkit:3.0.7"
    script:
        "scripts/getSRA.sh"

# docker pull pegi3s/sratoolkit:3.0.7