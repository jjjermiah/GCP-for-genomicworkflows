#### INCOMPLETE
rule all:
    input:
        expand("{outDir}/{libname}.s_filteredJunctions_fw.bam", outDir=config["outDir"], libname=config["libnames"])

rule filter_circular_junctions:
    input:
        chimeric_file="{inDir}/{libname}Chimeric.out.junction",
        filter_script="filterCirc.awk"
    output:
        filtered_file="{outDir}/{libname}.filteredJunctions.txt"
    shell:
        "cat {input.chimeric_file} | awk -f {input.filter_script} | sort | uniq -c | sort -k1,1rn > {output.filtered_file}"

rule convert_to_bed:
    input:
        filtered_file=rules.filter_circular_junctions.output.filtered_file,
        min_len=config["minLen"]
    output:
        bed_file="{outDir}/{libname}.filteredJunctions.bed"
    shell:
        "{scriptDir}/starCirclesToBed.pl {input.filtered_file} {config[minLen]} > {output.bed_file}"

rule filter_by_splice_sites:
    input:
        bed_file=rules.convert_to_bed.output.bed_file,
        filter_script="filterSpliceSiteCircles.pl"
    output:
        filtered_bed_file="{outDir}/{libname}.s_filteredJunctions.bed"
    shell:
        "{scriptDir}/{input.filter_script} {input.bed_file} > {output.filtered_bed_file}"

rule count_forward_spliced_reads:
    input:
        splice_circles_bed=rules.filter_by_splice_sites.output.filtered_bed_file,
        junction_file="{inDir}/{libname}SJ.out.tab"
    output:
        fw_spliced_bed_file="{outDir}/{libname}.s_filteredJunctions_fw.bed"
    shell:
        "{scriptDir}/nrForwardSplicedReads.pl {input.splice_circles_bed} {input.junction_file} > {output.fw_spliced_bed_file}"

rule create_bam:
    input:
        sam_file="{inDir}/{libname}Chimeric.out.sam"
    output:
        junction_bam="{outDir}/{libname}.bam",
        sorted_bam="{outDir}/{libname}.sorted.bam",
        indexed_bam="{outDir}/{libname}.sorted.bam.bai"
    shell:
        "samtools view -bS -o {output.junction_bam} {input.sam_file} && "
        "samtools sort -o {output.sorted_bam} {output.junction_bam} && "
        "samtools index {output.sorted_bam}"

configfile: "config.yaml"
