
#!/bin/bash

export bucket=orcestra-archive

echo "The snakemake input is: ${snakemake_input[@]}"
echo "The snakemake input is: ${snakemake_input}"

path_to_sra=${snakemake_input[0]}
path_to_SRA_folder=$(dirname "$path_to_sra")

export path_to_cachefiles="orcestra-archive/rawdata/cachefiles"
# echo out the paths
echo "The path to the cachefiles is: $path_to_cachefiles "

echo "The path to the SRA folder is: $path_to_SRA_folder "

echo "Copying the cachefiles to the SRA folder"
cp -r $path_to_cachefiles/* $path_to_SRA_folder
echo "The files in the SRA folder are:"
ls -la $path_to_SRA_folder

echo "The wildcards are: ${snakemake_wildcards}"
echo "The output is: ${snakemake_output[0]} ${snakemake_output[1]}"

cp -r $path_to_cachefiles/* $path_to_SRA_folder

echo "Changing directories to $(dirname $path_to_SRA_folder)"
cd $(dirname $path_to_SRA_folder)
export OUTPUT_DIR="/workdir/$(dirname ${snakemake_output[0]})"
echo "The output directory is: $OUTPUT_DIR"

fasterq-dump \
    --split-3 \
    --progress \
    --details \
    --threads ${snakemake[threads]} \
    ${snakemake_wildcards}

ls -la 

mkdir -p $OUTPUT_DIR
mv *.fast* $OUTPUT_DIR
ls -la $OUTPUT_DIR

echo "fasterq-dump complete"
