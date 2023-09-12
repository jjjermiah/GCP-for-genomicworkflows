
rule fastqc:
    input:
        fq=join("rawdata/{PROJECT_NAME}/", "FASTQ/{sample}_{split}.fastq.gz")
    output:
        html="QC/fastqc/{sample}_{split}.html",
        zip="QC/fastqc/{sample}_{split}fastqc.zip", # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
        tempfile=touch("{sample}_{split}_fastqc.done")
    params:
        extra = "--quiet"
    threads: 1
    log:
        "logs/fastqc/{sample}_{split}.log"
    resources:
        mem_mb = 1024
    wrapper:
        "v2.6.0/bio/fastqc"
