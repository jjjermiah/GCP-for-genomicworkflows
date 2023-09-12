# Rules to specifically get SRA data
def get_sample_refseqs(wildcards):
   # use wildcards.sample to get the 'size_in_gb' from the sra_metadata
    with checkpoints.create_ref_seq_list.get(**wildcards).output[0].open() as f:
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
        fq1="{any_path}/FASTQ/SRR{sra_acc}_1.fastq.gz",
        fq2="{any_path}/FASTQ/SRR{sra_acc}_2.fastq.gz"
    resources:
        machine_type = machines['med_cpu']['name']
    threads:
        threads = 7
    container:
        sratoolkit_docker
    script:
        "scripts/convertSRAtoFASTQ.sh"        
        # "hi.txt"

def get_disk_mb(wildcards):
    try:
        # use wildcards.sample to get the 'size_in_gb' from the sra_metadata
        # return the size in MB and add 10% for overhead
        sra_metadata = pd.read_csv("metadata/sra_metadata.csv")
        size_in_gb = sra_metadata[sra_metadata['sra_accession'] == wildcards.sample]['size_in_GB'].values[0]
        return int((size_in_gb * 1000)*1.1)
    except:
        raise Exception("Please create a sra_metadata.csv using the jupyter notebook in metadata/ or disable the resources parameter in download_sra")

rule download_sra:
    output:
        sra_zipped="{any_path}/SRA/SRR{sra_acc}.tar.gz",
        filelist="{any_path}/SRA/SRR{sra_acc}_files.txt"
    container:
        sratoolkit_docker
    threads:
        threads = 1
    resources:
        disk_mb = get_disk_mb
    script:
        "scripts/getSRA.sh"

# the output here is to remove ambiguity
checkpoint create_ref_seq_list:
    output:
        "{any_path}/SRA/ref_lists/SRR{sra_acc}_refs.csv"
    container:
        sratoolkit_docker
    retries: 5
    threads:
        1
    script:
        "scripts/createSRArefseq.sh"

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
