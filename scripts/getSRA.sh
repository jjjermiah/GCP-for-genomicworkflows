#!/bin/bash
echo "test"

# echo "${snakemake_wildcards['run']}" > ${snakemake_output[0]}
export SRA_FILES=$(dirname ${snakemake_output[0]})
echo $SRA_FILES
export SRA_DIR=$(dirname $SRA_FILES)
echo $SRA_DIR

/opt/sratoolkit.3.0.7-ubuntu64/bin/prefetch ${snakemake_wildcards['run']} -O $SRA_DIR
/opt/sratoolkit.3.0.7-ubuntu64/bin/align-info --ref ${snakemake_wildcards['run']} > ${snakemake_output[1]}
# pigz --version 