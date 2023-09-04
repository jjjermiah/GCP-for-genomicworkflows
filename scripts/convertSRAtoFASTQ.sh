#!/bin/bash

# output:
#     fq1="{PATH_TO_DATA}/FASTQ/{run}_1.fastq.gz",
#     fq2="{PATH_TO_DATA}/FASTQ/{run}_2.fastq.gz"


echo "The snakemake input is: ${snakemake_input[@]}"


path_to_sra_zip=${snakemake_input[0]}
path_to_SRA_folder=$(dirname "$path_to_sra_zip")
sra_zip=$(basename "$path_to_sra_zip")

# decompress sra_zip. It already has the folder structure we need assuming that the sra_zip was created using the script getSRA.sh
# and that all the data we are storing is in the same bucket 
export OUTPUT_DIR="/workdir/$(dirname ${snakemake_output[0]})"
echo "The output directory is: $OUTPUT_DIR"
mkdir -p $OUTPUT_DIR

tar -I pigz -xf $path_to_sra_zip 
cd $path_to_SRA_folder

echo "Changing directories to $(dirname $path_to_SRA_folder)"
cd $(dirname $path_to_SRA_folder)
echo "ls -la" 
ls -1 $path_to_SRA_folder

# Fasterq-dump is a tool from the SRA toolkit that converts SRA files to FASTQ files
fasterq-dump \
    --split-3 \
    --progress \
    --details \
    --threads ${snakemake[threads]} \
    ${snakemake_wildcards}

# The fastq files should be in the same directory that has the SRR_____ folder with all the files for the run 
echo "The current working directory is: $(pwd) and has the following fastq files:"
ls -la | grep "fastq"

# we want to compress it here using pigz multi-threaded compression 
pigz -9 *.fastq 

# move it to the output directory 
mv *.fastq.gz $OUTPUT_DIR