rule build_star_index:
    input:
        fasta=join(ref_path, "star", "genome.fa"),
        gtf=join(ref_path, "star", "annotation.gtf"),
    output:
        pipe(directory(join(ref_path, "star/index")))
    threads: 32
    resources:
        machine_type = "n1-highmem-32"
    conda:
        "../envs/star.yaml"
    params:
        extra="",
    log:
        "logs/star_index_genome.log",
    shell:
        """
        STAR \
        --runThreadN {threads} \
        --runMode genomeGenerate \
        --genomeFastaFiles {input.fasta} \
        --sjdbOverhang 99 \
        --sjdbGTFfile {input.gtf} \
        --genomeDir {output}
        """

rule star_pe_multi:
    input:
        fq1="rawdata/{PROJECT_NAME}/FASTQ/{sample}_1.fastq.gz",
        fq2="rawdata/{PROJECT_NAME}/FASTQ/{sample}_2.fastq.gz",
        # path to STAR reference genome index
        idx=directory(join(ref_path, "star/index"))
    output:
        # see STAR manual for additional output files
        aln="{procdata}/{PROJECT_NAME}/star/pe/{sample}/{sample}_pe_aligned.sam",
        log="{procdata}/{PROJECT_NAME}/logs/pe/{sample}/{sample}_Log.out",
        sj="{procdata}/{PROJECT_NAME}/star/pe/{sample}/{sample}_SJ.out.tab",
        chim_junc="{procdata}/{PROJECT_NAME}/star/pe/{sample}/{sample}_Chimeric.out.junction",
        unmapped=["{procdata}/{PROJECT_NAME}/star/pe/{sample}/{sample}_unmapped.1.fastq.gz","{procdata}/{PROJECT_NAME}/star/pe/{sample}/{sample}_unmapped.2.fastq.gz"],
    params:
        extra="--chimScoreMin 1 \
                --chimSegmentMin 20 \
                --alignIntronMax 1000000 \
                --alignTranscriptsPerReadNmax 100000 \
                --outFilterMismatchNoverReadLmax 0.02 \
                --outSAMtype BAM SortedByCoordinate \
                --outFilterMultimapNmax 2 \
                --chimOutType Junctions SeparateSAMold \
                --twopassMode Basic \
                --limitBAMsortRAM 500000000000"
    threads: 16
    wrapper:
        "v2.6.0/bio/star/align"
