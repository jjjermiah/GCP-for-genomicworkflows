
# NAME IDEAS:
GenomeCloud: A comprehensive repository for genome sequencing workflows on the Google Cloud Platform.
GCP-Genomics: Harnessing the power of Google Cloud Platform for genome sequencing and analysis.
GCP-GenomeSequencing: Centralized repository for genome sequencing workflows using Google Cloud services.
Genomic-CloudOps: Managing and streamlining genome sequencing processes with GCP and Snakemake.


This repository is (supposed to be) an all-in-one seq-processing hub utilizing google cloud tools and snakemake for efficient and cheap processing.

This work currently focusses on the pipelines associated with processing CCLE data by 
1) obtaining raw SRA data from ncbi
2) converting SRA to FASTQ
3) Quality control?
4) ??????


# Sunday September 3 
functionalities include:
    - querying NCBI data bases for which SRA metadata
    - once .sra files exist in GCS buckets, converting them to fastq files
    - some (limited) methods of compressing the fastq files
