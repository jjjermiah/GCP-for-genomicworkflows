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

rule download_sra:
    output:
        sra_zipped="{any_path}/SRA/SRR{sra_acc}.tar.gz",
        filelist="{any_path}/SRA/SRR{sra_acc}_files.txt"
    container:
        sratoolkit_docker
    threads:
        threads = 1
    # resources:
    #     disk_mb = get_disk_mb
    script:
        "../scripts/getSRA.sh"
# Rules to specifically get SRA data
def get_sample_refseqs(wildcards):
   # use wildcards.sample to get the 'size_in_gb' from the sra_metadata
    with checkpoints.get_sra_ref_seqs.get(**wildcards).output[0].open() as f:
        # get first element of comma separated line
        # first_element=f.readline().split(",")[0]
        refseqs = ["tools/cachefiles/" + line.split(",")[0] for line in f]
    # return the size in MB and add 10% for overhead
    return refseqs

# by forcing the SRR prefix to sra_acc, this will still work if the user
# inputs the full SRR prefix to any rule that requires an SRA accession
# it just prevents this rule from being used for non-SRA fastq files
rule sra_to_fastq:
    input:
        get_sample_refseqs,
        "{any_path}/SRA/SRR{sra_acc}.tar.gz"
    output:
        fq1="{any_path}/FASTQ/SRR{sra_acc}_1.fastq",
        fq2="{any_path}/FASTQ/SRR{sra_acc}_2.fastq"
    resources:
        machine_type = machines['high_mem']['name'],
    group: "sra_fastq"
    threads:
        threads = 8
    container:
        sratoolkit_docker
    log:
        stderr=os.path.join(
            config["log_dir"], "{any_path}/FASTQ/SRR", "{sra_acc}", "fasterq_dump.stderr.log"
        ),
        stdout=os.path.join(
            config["log_dir"], "{any_path}/FASTQ/SRR", "{sra_acc}", "fasterq_dump.stdout.log"
        ),
    script:
        "../scripts/convertSRAtoFASTQ.sh"        
        # "hi.txt"

rule compress_fastq:
    "Compress fastq inplace with pigz at best (9) compression level."
    input:
        fq="{any_path}/FASTQ/SRR{sra_acc}_{split}.fastq",
    output:
        fq="{any_path}/FASTQ/SRR{sra_acc}_{split}.fastq.gz"
    group: "sra_fastq"
    threads:
        threads = 8
    container:
        pigz_docker
    log:
        stderr=os.path.join(
            config["log_dir"], "{any_path}/FASTQ/SRR","{sra_acc}", "compress_fastq_{split}.stderr.log"
        ),
        stdout=os.path.join(
            config["log_dir"], "{any_path}/FASTQ/SRR", "{sra_acc}", "compress_fastq_{split}.stdout.log"
        ),
    shell:
        """
        pigz --best --processes {threads} {input.fq}; \
        1> {log.stdout} 2> {log.stderr}; \
        touch {output}; \
        """

# the output here is to remove ambiguity
checkpoint get_sra_ref_seqs:
    output:
        "{any_path}/SRA/ref_lists/SRR{sra_acc}_refs.csv"
    container:
        "docker://ncbi/sra-tools"
    retries: 5
    threads:
        1
    shell:
        "/usr/local/bin/align-info {wildcards.sra_acc} | cut -d ',' -f1 >${snakemake_output[0]}"
    # script:
    #     "../scripts/create_refseq_list.sh"    

rule download_refseqs:
    output:
        "{any_path}/SRA/cachefiles/{refseq}"
    resources:
        machine_type = machines['med_cpu']['name']
    retries: 5
    threads:
        1
    shell:
        "wget https://sra-download.ncbi.nlm.nih.gov/traces/refseq/{wildcards.refseq} -O {output}"
