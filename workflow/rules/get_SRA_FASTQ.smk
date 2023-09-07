# Rules to specifically get SRA data



# def get_sample_refseqs(wildcards):
#     # use wildcards.sample to get the 'size_in_gb' from the sra_metadata
#     with checkpoints.create_ref_seq_list.get(sample=wildcards.sample).output[0].open() as f:
#         # get first element of comma separated line
#         # first_element=f.readline().split(",")[0]
#         refseqs = ["cachefiles/" + line.split(",")[0] for line in f]
#     # return the size in MB and add 10% for overhead
#     return refseqs

# rule convertSRAtoFASTQ:
#     input:
#         "SRA/{sample}.tar.gz",
#         get_sample_refseqs
#         # expand("cachefiles_list.done", PATH_TO_DATA=PATH_TO_DATA),
#     output:
#         fq1=join(rawdata_path, "FASTQ/{sample}_1.fastq.gz"),
#         fq2=join(rawdata_path, "FASTQ/{sample}_2.fastq.gz")
#     resources:
#         machine_type = high_cpu
#     threads:
#         threads = 7
#     container:
#         sratoolkit_docker
#     script:
#         "scripts/convertSRAtoFASTQ.sh"        
#         # "hi.txt"

# def get_disk_mb(wildcards):
#     # use wildcards.sample to get the 'size_in_gb' from the sra_metadata
#     size_in_gb = sra_metadata[sra_metadata['sample_accession'] == wildcards.sample]['size_in_GB'].values[0]

#     # return the size in MB and add 10% for overhead
#     return int((size_in_gb * 1000)*1.1)

# rule download_sra:
#     output:
#         sra_zipped="SRA/{sample}.tar.gz",
#         filelist="SRA/{sample}_files.txt"
#     container:
#         sratoolkit_docker
#     threads:
#         threads = 1
#     resources:
#         disk_mb = get_disk_mb
#     script:
#         "scripts/getSRA.sh"

# # the output here is to remove ambiguity
# checkpoint create_ref_seq_list:
#     output:
#        sample_refs="ref_lists/{sample}_refs.csv"
#     container:
#         sratoolkit_docker
#     retries: 5
#     threads:
#         1
#     script:
#         "scripts/createSRArefseq.sh"

# rule download_refseqs:
#     output:
#         "cachefiles/{refseq}"
#     resources:
#         machine_type = med_cpu
#     retries: 5
#     threads:
#         1
#     shell:
#         "wget https://sra-download.ncbi.nlm.nih.gov/traces/refseq/{wildcards.refseq} -O {output}"
