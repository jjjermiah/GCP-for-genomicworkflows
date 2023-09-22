# CIRI CIRCExplorer

rule CIRI2:
    input:
        sam=join("procdata/{PROJECT_NAME}/", "alignment/{sample}.sam"),
        gtf=join(ref_path, "bwa", "annotation.gtf"),
        idx=multiext(join(ref_path, "bwa", "genome"), ".fa", ".amb", ".ann", ".bwt", ".pac", ".sa")
    output:
        CIRI2= "results/{PROJECT_NAME}/CIRI2/{sample}.tsv"
    threads: 16
    log: "log/{PROJECT_NAME}/CIRI2/{sample}.log"
    container:
        "docker://andremrsantos/ciri2:latest"
    shell:
        """CIRI2 \
        --thread_num {threads} \
        --in {input.sam} \
        --anno {input.gtf} \
        --out {output.CIRI2} \
        --ref_file {input.idx[0]} \
        --log {log}"""

rule circRNA_finder:
    input:
        aln="procdata/{PROJECT_NAME}/star/pe/{sample}/{sample}_pe_aligned.sam",
        sj="procdata/{PROJECT_NAME}/star/pe/{sample}/{sample}_SJ.out.tab",
        chim_junc="procdata/{PROJECT_NAME}/star/pe/{sample}/{sample}_Chimeric.out.junction",
        # log="procdata/{PROJECT_NAME}/logs/pe/{sample}/{sample}_Log.out",
        # unmapped=[
        #         "procdata/{PROJECT_NAME}/star/pe/{sample}/{sample}_unmapped.1.fastq.gz",
        #         "procdata/{PROJECT_NAME}/star/pe/{sample}/{sample}_unmapped.2.fastq.gz"],
    output:
        filteredjunctions="results/{PROJECT_NAME}/circRNA_finder/{sample}/{sample}_filteredJunctions.bed",
        GT_AG_filteredjunctions="results/{PROJECT_NAME}/circRNA_finder/{sample}/{sample}_s_filteredJunctions.bed",
        fw_filteredjunctions="results/{PROJECT_NAME}/circRNA_finder/{sample}/{sample}_s_filteredJunctions_fw.bed",
        # sorted_bam="results/{PROJECT_NAME}/circRNA_finder/{sample}/{sample}_pe_aligned.sorted.bam",
        # sorted_indexed_bam="results/{PROJECT_NAME}/circRNA_finder/{sample}/{sample}_pe_aligned.sorted.bam.bai"
    conda:
        "../envs/circRNA_finder.yaml"
    script:
        "../scripts_circRNA_finder/runPostProcessStarAlignment.sh"
        
# a) _filteredJunctions.bed: 
    # A bed file with all circular junctions found by the pipeline. 
    # The score column indicates the number reads spanning each junction.

# b) _s_filteredJunctions.bed: 
    # A bed file with those juction in (a) that are flanked by GT-AG splice sites. 
    # The score column indicates the number reads spanning each junction.

# c) _s_filteredJunctions_fw.bed: 
    # A bed file with the same circular junctions as in file (b), 
    # but here the score column gives the average number of forward spliced reads at both splice sites around each circular junction.

# d) (Sorted and indexed) bam file with all chimeric reads identified by STAR. 
    # The circRNA junction spanning reads are a subset of these.

