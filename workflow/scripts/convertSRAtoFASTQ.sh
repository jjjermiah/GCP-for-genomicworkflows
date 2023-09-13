#!/bin/bash
# Load files 
path_to_sra_zip=${snakemake_input['sra']}
path_to_SRA_folder=$(dirname "$path_to_sra_zip")
sra_zip=$(basename "$path_to_sra_zip")

echo "Sample accession is: ${snakemake_wildcards['sample']}"

# echo "The snakemake input is: ${snakemake_input[@]}"

export OUTPUT_DIR="/workdir/$(dirname ${snakemake_output[0]})"
echo "The output directory is: $OUTPUT_DIR"
mkdir -p $OUTPUT_DIR
mkdir -p $(dirname ${snakemake_log['stderr']})
mkdir -p $(dirname ${snakemake_log['stdout']})


# decompress sra_zip. It already has the folder structure we need assuming that the sra_zip was created using the script getSRA.sh
# and that all the data we are storing is in the same bucket 
echo "decompressing $path_to_sra_zip"
tar -xzf $path_to_sra_zip && rm $path_to_sra_zip

echo "Changing directories to $path_to_SRA_folder"
cd $path_to_SRA_folder

echo "ls -lahR" 
ls -lahR

echo "Beginning Fasterq-dump"
fasterq-dump \
    --split-3 \
    --progress \
    --details \
    --threads ${snakemake[threads]} \
    ${snakemake_wildcards['sample']} 

# # The fastq files should be in the same directory that has the SRR_____ folder with all the files for the run 
echo "The current working directory is: $(pwd) and has the following fastq files:"
ls -lah | grep "fastq"

# # TODO:: use the snakemake output as redirect instead of the hardcoded path
# echo "Compressing fastq files"
echo "pigz -9 --processes ${snakemake[threads]} *.fastq "
# # we want to compress it here using pigz multi-threaded compression 
pigz -9 --processes ${snakemake[threads]} *.fastq 

mv *.fastq.gz $OUTPUT_DIR