

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


# Using Google Cloud Platform 
In this workflow, we utilize the following products from Google Cloud:
1) Google VM Instances
- As a cheap workspace (e2-micro is free) 
- Where we modify code and submit jobs
2) Google Kubernetes Clusters 
- To do a lot of the heavy processing 
- Snakemake has built-in functionality to use kubernetes clusters for jobs
3) Google Cloud Storage Buckets
- to store raw, processed, and meta data 
- also store our workspace for cheap! (for objects not pushed to github)
4) Google Dataflow
- simple task of using GZIP on fastq files

Looking to learn and incorporate:
5) Google BigQuery
Eventually: 
6) Google Cloud Functions | Cloud Run 
- adding metadata to google cloud storage bucket?
- performing quality control? 
- creating vizualizations 
- postprocessing/archival 
7) Batch
- with the deprecation of Google Cloud Life Sciences for the preferred Batch, Snakemake is working on a collaboration to bring Batch as an executor for Snakemake. 
    - this might provide added benefit as there might not be a need to setup any kubernetes cluster at all
    - all that is required is an idea of what resources each job/rule in snakemake will need, and GCP Batch will handle the rest
8) Google Workflows
https://console.cloud.google.com/workflows/
- seems cool
