# Note on job groups https://snakemake.readthedocs.io/en/stable/executing/grouping.html#job-grouping: 
# From Snakemake 7.11 on, Snakemake will request resources for groups by summing across jobs that can be run in parallel, and taking the max of jobs run in series. 



# def get_disk_mb(wildcards):
#     try:
#         # use wildcards.sample to get the 'size_in_gb' from the sra_metadata
#         # return the size in MB and add 10% for overhead
#         sra_metadata = pd.read_csv(join("../metadata/", "sra_metadata.csv"))
#         size_in_gb = sra_metadata[sra_metadata['sra_accession'] == wildcards.sample]['size_in_GB'].values[0]
#         return int((size_in_gb * 1000)*1.1)
#     except:
#         raise Exception("Please create a sra_metadata.csv using the jupyter notebook in metadata/ or disable the resources parameter in download_sra")

# rule download_sra:
#     output:
#         sra_zipped="{any_path}/SRA/SRR{sra_acc}.tar.gz"
#     container:
#         sratoolkit_docker
#     threads:
#         threads = 1
#     script:
#         "../scripts/getSRA.sh"
        
# Rules to specifically get SRA data
def get_sample_refseqs(wildcards):
   # use wildcards.sample to get the 'size_in_gb' from the sra_metadata
    with checkpoints.get_sra_ref_seqs.get(**wildcards).output[0].open() as f:
        refseqs = ["rawdata/cachefiles/" + line.strip() for line in f]
    # return the size in MB and add 10% for overhead
    return refseqs

# by forcing the SRR prefix to sra_acc, this will still work if the user
# inputs the full SRR prefix to any rule that requires an SRA accession
# it just prevents this rule from being used for non-SRA fastq files
rule sra_to_fastq:
    input:
        refs=get_sample_refseqs,
        sra="rawdata/{PROJECT_NAME}/SRA/{sample}.tar.gz"
    output:
        fq1="rawdata/{PROJECT_NAME}/FASTQ/{sample}_1.fastq.gz",
        fq2="rawdata/{PROJECT_NAME}/FASTQ/{sample}_2.fastq.gz"
    resources:
        machine_type = machines['high_mem']['name'],
    threads:
        threads = 32
    container:
        sratoolkit_docker
    log:
        stderr=os.path.join(
            config["log_dir"], "rawdata/{PROJECT_NAME}/FASTQ/{sample}_fasterq_dump.stderr.log"
        ),
        stdout=os.path.join(
            config["log_dir"], "rawdata/{PROJECT_NAME}/FASTQ/{sample}_fasterq_dump.stdout.log"
        ),
    script:
        "../scripts/convertSRAtoFASTQ.sh"        

checkpoint get_sra_ref_seqs:
    output:
        "rawdata/ref_lists/{sample}_refs.csv"
    container:
        "docker://jjjermiah/sratoolkit:0.2"
    retries: 5
    threads:
        1
    script:
        "../scripts/create_refseq_list.sh"    

rule download_refseqs:
    output:
        "rawdata/cachefiles/{refseq}"
    resources:
        machine_type = machines['med_cpu']['name']
    retries: 5
    threads:
        1
    shell:
        "wget https://sra-download.ncbi.nlm.nih.gov/traces/refseq/{wildcards.refseq} -O {output}"


# rule compress_fastq:
#     "Compress fastq inplace with pigz at best (9) compression level."
#     input:
#         fq="{any_path}/FASTQ/SRR{sra_acc}_{split}.fastq",
#     output:
#         fq="{any_path}/FASTQ/SRR{sra_acc}_{split}.fastq.gz"
#     threads:
#         threads = 8
#     container:
#         pigz_docker
#     log:
#         stderr=os.path.join(
#             config["log_dir"], "{any_path}/FASTQ/SRR","{sra_acc}", "compress_fastq_{split}.stderr.log"
#         ),
#         stdout=os.path.join(
#             config["log_dir"], "{any_path}/FASTQ/SRR", "{sra_acc}", "compress_fastq_{split}.stdout.log"
#         ),
#     shell:
#         """
#         pigz --best --processes {threads} {input.fq}; \
#         1> {log.stdout} 2> {log.stderr}; \
#         touch {output}; \
#         """