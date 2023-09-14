


################################################################################################
################ ENSEMBL

# url_prefix = f"ftp://ftp.ensembl.org/pub/{branch}release-{release}/fasta/{species}/{datatype}/{species.capitalize()}.{spec}"
# reference_genome = config["ref"]
# # ref_path = f"reference_genomes/{reference_genome['SPECIES']}/release-{reference_genome['RELEASE']}/{reference_genome['BUILD']}/"

# rule get_ENSEMBL_genome:
#     output:
#         "{reference_build_spec}/genome.fa",
#     params:
#         species=reference_genome['SPECIES'],
#         datatype=reference_genome['DATATYPE'],
#         build=reference_genome['BUILD'],
#         release=reference_genome['RELEASE'],
#     log:
#         "{reference_build_spec}/logs/get_genome.log",
#     wrapper:
#         "v2.6.0/bio/reference/ensembl-sequence"

# rule get_ENSEMBL_annotation:
#     output:
#         "{reference_build_spec}/annotation.gtf",
#     params:
#         species=reference_genome['SPECIES'],
#         datatype=reference_genome['DATATYPE'],
#         build=reference_genome['BUILD'],
#         release=reference_genome['RELEASE'],
#     log:
#         "{reference_build_spec}/logs/get_annotation.log",
#     wrapper:
#         "v2.6.0/bio/reference/ensembl-annotation"

# rule get_annotation_gz:
#     output:
#         "{reference_build_spec}/annotation.gtf.gz",
#     params:
#         species=reference_genome['SPECIES'],
#         datatype=reference_genome['DATATYPE'],
#         build=reference_genome['BUILD'],
#         release=reference_genome['RELEASE'],
#         flavor="",  # optional, e.g. chr_patch_hapl_scaff, see Ensembl FTP.
#         # branch="plants",  # optional: specify branch
#     log:
#         "{reference_build_spec}/logs/get_annotation.log",
#     wrapper:
#         "v2.6.0/bio/reference/ensembl-annotation"


# ref_path = f"reference_genomes/ENSEMBL/{reference_genome['SPECIES']}/{reference_genome['BUILD']}/release-{reference_genome['RELEASE']}/"

def get_gencode_annotation(wildcards):
    if wildcards.ref_build == "GRCh37":
        if wildcards.species == "homo_sapiens" or "human" or "Gencode_human":
            ftp = f"ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{wildcards.gencode_release}/GRCh37_mapping/gencode.v{wildcards.gencode_release}lift37.annotation.gtf.gz"
    else:
        ftp = f"ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{wildcards.gencode_release}/gencode.v{wildcards.gencode_release}.annotation.gtf.gz"
    return FTP.remote(ftp, immediate_close=True, keep_local=True)

rule getGENCODEannotation:
    input:
        get_gencode_annotation
    output:
        gencode_annotation_file="reference_genomes/GENCODE/{species}/{ref_build}/release-{gencode_release}/annotation.gtf"
        # gencode_annotation_file="references/Gencode_human/gencode.v{gencode_release}.annotation.gtf"
    shell:
        "gzip -d -c {input} > {output}"

def get_gencode_genome(wildcards):
    if wildcards.ref_build == "GRCh37":
        if wildcards.species == "homo_sapiens" or "human" or "Gencode_human":
            ftp = f"ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{wildcards.gencode_release}/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz"
    else:
        ftp = f"ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{wildcards.gencode_release}/GRCh38.primary_assembly.genome.fa.gz"
    return FTP.remote(ftp, immediate_close=True, keep_local=True)

rule getGENCODEgenome:
    input:
        get_gencode_genome
    output:
        gencode_genome_file="reference_genomes/GENCODE/{species}/{ref_build}/release-{gencode_release}/genome.fa",
    shell:
        "gzip -d -c {input} > {output}"

# ## # Use the reference genomes from the google cloud bucket
# # Deprecating this. Looks like google-cloud-lifesciences api is being deprecated, not sure whats going to happen to this bucket. 
# # Useful if wishing to use data from genomics-public-data
# rule get_ref_genome:
#     input:
#         GS.remote(
#             expand("genomics-public-data/references/hg19/chr{chrom}.fa.gz",  # use the public google bucket from 
#             chrom=[str(i) for i in range(1, 23)] + ["X", "Y", "M"]
#                 )
#             )
#     params:
#         output_dir="hg19"
#     shell:
#         """
#         mkdir -p {params.output_dir}
#         zcat {input} > {output}
#         """