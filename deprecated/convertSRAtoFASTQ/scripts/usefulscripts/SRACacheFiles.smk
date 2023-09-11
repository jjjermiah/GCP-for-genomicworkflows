import pandas as pd
import_df = pd.read_csv("metadata/full_list_df.csv")
pathdf_ = import_df[['CellLine', 'Run', 'GCS_RNA_SRA_dirpaths']]


sra_docker = "docker://jjjermiah/sratools:0.5"



# SUBSET
df = pathdf_

df = df.iloc[:1]

cell_lines = df['CellLine'].values.tolist()
runs = df['Run'].values.tolist()
paths = df['GCS_RNA_SRA_dirpaths'].values.tolist()



checkpoint create_all_refseqs: 
    input:
        expand("RNA/refseqlists/{run}_refseqs.csv", run=runs)
    output:
        "SRA/all_refseqs.csv"
    shell:
        "touch scripts/create_all_refseqs.py"

# the output here is to remove ambiguity
checkpoint create_ref_seq_list:
    # input:
    #     "RNA/SRA/{run}/{run}.sra"
    output:
       "RNA/refseqlists/{run}_refseqs.csv",
    container:
        sra_docker
    threads:
        1
    script:
        "scripts/create_refseq_list.sh"


# input function for the rule aggregate to make sure all required cache files exist
def get_all_cache_files(wildcards):
    # decision based on content of output file
    # Important: use the method open() of the returned file!
    # This way, Snakemake is able to automatically download the file if it is generated in
    # a cloud environment without a shared filesystem.
    with checkpoints.create_all_refseqs.get(**wildcards).output[0].open() as f:
        # read in every line of the file and return it as a list
        return ["SRA/cachefiles/" + line.strip() for line in f]

# rule aggregate:
#     input:
#         get_all_cache_files
#     output:
#         touch("SRA/cachefiles_list.done")
#     shell:
#         "echo {input}"

# rule download_refseqs:
#     output:
#         "SRA/cachefiles/{refseq}"
#     threads:
#         1
#     shell:
#         # wget {wildcards.refseq} into {output}
#         "wget https://sra-download.ncbi.nlm.nih.gov/traces/refseq/{wildcards.refseq} -O {output}"

# def get_run_cache_files(wildcards):
#     with checkpoints.create_ref_seq_list.get(**wildcards).output[0].open() as f:
#         return ["SRA/cachefiles/" + line.strip() for line in f]


