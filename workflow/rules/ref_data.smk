

################################################################################################
################ ENSEMBL

# url_prefix = f"ftp://ftp.ensembl.org/pub/{branch}release-{release}/fasta/{species}/{datatype}/{species.capitalize()}.{spec}"

rule get_genome:
    output:
        "{reference_build_spec}/genome.fa",
    params:
        species=REF_SPECIES,
        datatype=REF_DATATYPE,
        build=REF_BUILD,
        release=REF_RELEASE,
    log:
        "{reference_build_spec}/logs/get_genome.log",
    wrapper:
        "v2.6.0/bio/reference/ensembl-sequence"

rule get_annotation:
    output:
        "{reference_build_spec}/annotation.gtf",
    params:
        species=REF_SPECIES,
        release=REF_DATATYPE,
        build=REF_BUILD,
    log:
        "{reference_build_spec}/logs/get_annotation.log",
    wrapper:
        "v2.6.0/bio/reference/ensembl-annotation"

rule get_annotation_gz:
    output:
        "{reference_build_spec}/annotation.gtf.gz",
    params:
        species="homo_sapiens",
        release="105",
        build="GRCh37",
        flavor="",  # optional, e.g. chr_patch_hapl_scaff, see Ensembl FTP.
        # branch="plants",  # optional: specify branch
    log:
        "{reference_build_spec}/logs/get_annotation.log",
    wrapper:
        "v2.6.0/bio/reference/ensembl-annotation"



# def get_gencode_annotation(wildcards):
#     if reference_genome == "GRCh37":
#         ftp = f"ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{wildcards.gencode_release}/GRCh37_mapping/gencode.v{wildcards.gencode_release}lift37.annotation.gtf.gz"
#     else:
#         ftp = f"ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{wildcards.gencode_release}/gencode.v{wildcards.gencode_release}.annotation.gtf.gz"
#     return FTP.remote(ftp, immediate_close=True, keep_local=True)

# rule annotation:
#     input:
#         get_gencode_annotation
#     output:
#         gencode_annotation_file="references/Gencode_human/gencode.v{gencode_release}.annotation.gtf"
#     shell:
#         "gzip -d -c {input} > {output}"


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