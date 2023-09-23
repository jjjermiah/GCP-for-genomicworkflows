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

from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

def get_gencode_annotation(wildcards):
    if wildcards.ref_build == "GRCh37":
        ftp = f"ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{wildcards.gencode_release}/GRCh37_mapping/gencode.v{wildcards.gencode_release}lift37.annotation.gtf.gz"
    else:
        ftp = f"ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{wildcards.gencode_release}/gencode.v{wildcards.gencode_release}.annotation.gtf.gz"
    return HTTP.remote(ftp,  keep_local=True)

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
            ftp_genome = f"ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{wildcards.gencode_release}/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz"
    else:
        ftp_genome = f"ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{wildcards.gencode_release}/GRCh38.primary_assembly.genome.fa.gz"
    return HTTP.remote(ftp_genome, keep_local=True)

rule getGENCODEgenome:
    input:
        get_gencode_genome
    output:
        gencode_genome_file="reference_genomes/GENCODE/{species}/{ref_build}/release-{gencode_release}/genome.fa",
    shell:
        "gzip -d -c {input} > {output}"

def get_gencode_transcriptome(wildcards):
    if wildcards.ref_build == "GRCh37":
        ftp = f"ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{wildcards.gencode_release}/GRCh37_mapping/gencode.v{wildcards.gencode_release}lift37.transcripts.fa.gz"
    else:
        ftp = f"ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{wildcards.gencode_release}/gencode.v{wildcards.gencode_release}.transcripts.fa.gz"
    return HTTP.remote(ftp, keep_local=True)

rule getGENCODEtranscriptome:
    input:
        get_gencode_transcriptome
    output:
        gencode_genome_file="reference_genomes/GENCODE/{species}/{ref_build}/release-{gencode_release}/transcriptome.fa",
    shell:
        "gzip -d -c {input} > {output}"