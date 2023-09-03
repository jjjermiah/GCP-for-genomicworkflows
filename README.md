# Combining Google Cloud Platform tools and Genomic workflow pipelines 

Disclaimer: I am actively learning and exploring the Google Cloud platform, Bioinformatics tools, and Snakemake. 
As a beginner in these domains, I aim to share my journey and progress with others who may find it helpful. 
Please note that the content provided here may not always reflect the most advanced or expert-level knowledge. 
Instead, it represents my ongoing learning process as I figure some things out and test my understanding.  
I encourage feedback, any collaboration, and suggestions! 



# Using Google Cloud Platform 
In this repo, we utilize the following products from Google Cloud:
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
