rule star_index:
    input:
        fasta=join(ref_path, "genome.fa"),
        gtf=join(ref_path, "annotation.gtf"),
    output:
        directory(join(ref_path, "STAR_INDEX"))
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
        # use a list for multiple fastq files for one sample
        # usually technical replicates across lanes/flowcells
        # fq1=["reads/{sample}_R1.1.fastq", "reads/{sample}_R1.2.fastq"],
        # # paired end reads needs to be ordered so each item in the two lists match
        # fq2=["reads/{sample}_R2.1.fastq", "reads/{sample}_R2.2.fastq"],  #optional
        fq1="rawdata/{PROJECT_NAME}/FASTQ/{sample}_1.fastq.gz",
        fq2="rawdata/{PROJECT_NAME}/FASTQ/{sample}_2.fastq.gz",
        # path to STAR reference genome index
        idx=directory(join(ref_path, "STAR_INDEX"))
    output:
        # see STAR manual for additional output files
        aln="{procdata}/{PROJECT_NAME}/star/pe/{sample}/{sample}_pe_aligned.sam",
        log="{procdata}/{PROJECT_NAME}/logs/pe/{sample}/{sample}_Log.out",
        sj="{procdata}/{PROJECT_NAME}/star/pe/{sample}/{sample}_SJ.out.tab",
        chim_junc="{procdata}/{PROJECT_NAME}/star/pe/{sample}/{sample}_Chimeric.out.junction",
        unmapped=["{procdata}/{PROJECT_NAME}/star/pe/{sample}/{sample}_unmapped.1.fastq.gz","{procdata}/{PROJECT_NAME}/star/pe/{sample}/{sample}_unmapped.2.fastq.gz"],
    conda:
        "../envs/star.yaml"
    container:
        "docker://quay.io/biocontainers/star:2.7.8a--h9ee0642_1"
    params:
        # optional parameters
        extra="--chimScoreMin 1 \
                --chimSegmentMin 20 \
                --alignIntronMax 1000000 \
                --alignTranscriptsPerReadNmax 100000 \
                --outFilterMismatchNoverReadLmax 0.02 \
                --outSAMtype BAM SortedByCoordinate \
                --outFilterMultimapNmax 2 \
                --chimOutType Junctions SeparateSAMold \
                --twopassMode Basic"
    threads: 32
    script:
        "../scripts/map_star.py"
    # shell:
    #     """
    #     STAR \
    #     --runThreadN {threads} \
    #     --genomeDir {input.idx} \
    #     --readFilesIn {input.fq1} {input.fq2} \
    #     --readFilesCommand zcat \
    #     --outFileNamePrefix {params.outFileNamePrefix} \
    #     --outSAMattributes All \
    #     --outStd BAM_Unsorted \
    #     > {snakemake.output.aln}"
    #     " {log}"
    #     > {output.bam}; \
    #     2> {log.stderr}
    #     """

    # wrapper:
    #     "v2.6.0/bio/star/align"
    # """
    # STAR \
    # --genomeDir $genomeDir \
    # --readFilesCommand gunzip -c --readFilesIn $inFile1 $inFile2 \
    # --runThreadN $nThreads \
    # --chimSegmentMin $chimSegMin \
    # --chimScoreMin 1 \
    # --alignIntronMax $alignIntronMax \
    # --outFilterMismatchNoverReadLmax $maxMismatchFraction \
    # --alignTranscriptsPerReadNmax $alignTxPerReadMax \
    # --twopassMode Basic \
    # --outSAMtype BAM SortedByCoordinate \
    # --chimOutType Junctions SeparateSAMold \
    # --outFilterMultimapNmax 2 \
    # --outFileNamePrefix $outPrefix
    # """
