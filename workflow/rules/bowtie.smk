
# rule bowtie2_build:
#     input:
#         ref=join(ref_path, "genome.fa")
#     output:
#         multiext(join(ref_path, "genome.fa"), "", ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
#     conda:
#         "envs/bowtie2.yaml"
#     threads: 
#         8
#     shell:
#         "v2.6.0/bio/bowtie2/build"