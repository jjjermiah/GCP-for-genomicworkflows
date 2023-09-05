#!/bin/bash
echo "test"

# # echo "${snakemake_wildcards['run']}" > ${snakemake_output[0]}
# export SRA_FILES=$(dirname ${snakemake_output[0]})
# echo $SRA_FILES
# export SRA_DIR=$(dirname $SRA_FILES)
# echo $SRA_DIR

# export SRA_DIR=$(dirname $SRA_FILES)
# echo $SRA_DIR
echo "${snakemake_output[0]}"

export SRA_DIR=$(dirname ${snakemake_output[0]})
echo "$SRA_DIR"

# prefetch files to SRA_DIR
/opt/sratoolkit.3.0.7-ubuntu64/bin/prefetch ${snakemake_wildcards['run']} -O $SRA_DIR

# 
# echo "Getting Reference List"

# create reference list for alignment in case incorrect download or reference files cannot be accessed
# /opt/sratoolkit.3.0.7-ubuntu64/bin/align-info --ref ${snakemake_wildcards['run']} > ${snakemake_output[1]}
# echo "Done Getting Reference List"

echo "Getting List of Files Downloaded"
# create list of files downloaded
ls -1 $SRA_DIR/${snakemake_wildcards['run']} > ${snakemake_output[1]}

# remove all files in the SRA_DIR except the SRA files
echo "Removing all files in SRA_DIR except SRA files"
find $SRA_DIR/${snakemake_wildcards['run']} -type f ! -name "*${snakemake_wildcards['run']}*" -delete

# compress the folder with the SRA files
echo "Compressing SRA Files"
# pigz -9 -r $SRA_DIR/${snakemake_wildcards['run']} -c > ${snakemake_output[0]}
tar -I pigz -cf ${snakemake_output[0]} $SRA_DIR/${snakemake_wildcards['run']}/

# To decompress:
# tar -I pigz -xf compressed_folder.tar.gz -C /path/to/destination/directory
 