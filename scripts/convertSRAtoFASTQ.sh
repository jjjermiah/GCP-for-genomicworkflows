#!/bin/bash

# output:
#     fq1="{PATH_TO_DATA}/FASTQ/{run}_1.fastq.gz",
#     fq2="{PATH_TO_DATA}/FASTQ/{run}_2.fastq.gz"

# Load files 
path_to_sra_zip=${snakemake_input[0]}
path_to_SRA_folder=$(dirname "$path_to_sra_zip")
sra_zip=$(basename "$path_to_sra_zip")

echo "${snakemake_wildcards['run']}"

echo "The snakemake input is: ${snakemake_input[@]}"

export OUTPUT_DIR="/workdir/$(dirname ${snakemake_output[0]})"
echo "The output directory is: $OUTPUT_DIR"
mkdir -p $OUTPUT_DIR

# decompress sra_zip. It already has the folder structure we need assuming that the sra_zip was created using the script getSRA.sh
# and that all the data we are storing is in the same bucket 
echo "decompressing $path_to_sra_zip"
tar -xzf $path_to_sra_zip && rm $path_to_sra_zip

echo "Changing directories to $path_to_SRA_folder"
cd $path_to_SRA_folder

echo "ls -la" 
ls -lah
# shell:
#     """
#     fasterq-dump {params.outdir} --outdir {params.outdir} \
#         --mem {resources.mem_mb}MB --threads {threads} \
#         --temp {resources.tmpdir} \
#         1> {log.stdout} 2> {log.stderr}; \
#     touch {output.flag}
#     """
# Fasterq-dump is a tool from the SRA toolkit that converts SRA files to FASTQ files
echo "Beginning Fasterq-dump"
fasterq-dump \
    --split-3 \
    --progress \
    --details \
    --threads ${snakemake[threads]} \
    ${snakemake_wildcards['run']} \
    1> ${snakemake_log['stdout']} 2> ${snakemake_log['stderr']} 

# The fastq files should be in the same directory that has the SRR_____ folder with all the files for the run 
echo "The current working directory is: $(pwd) and has the following fastq files:"
ls -lah | grep "fastq"


# TODO:: use the snakemake output as redirect instead of the hardcoded path
echo "Compressing fastq files"
echo "pigz -9 --processes ${snakemake[threads]} *.fastq "
# we want to compress it here using pigz multi-threaded compression 
pigz -9 --processes ${snakemake[threads]} *.fastq 

echo "Moving compressed fastq files to $OUTPUT_DIR"
# move it to the output directory 
mv *.fastq.gz $OUTPUT_DIR