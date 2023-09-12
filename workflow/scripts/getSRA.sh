#!/bin/bash

export SRA_DIR=$(dirname ${snakemake_output[0]})

# prefetch files to SRA_DIR
/opt/sratoolkit.3.0.7-ubuntu64/bin/prefetch "SRR${snakemake_wildcards['sra_acc']}" -O $SRA_DIR

# # remove all files in the SRA_DIR except the SRA files
# echo "Removing all files in SRA_DIR except SRA files"
# find $SRA_DIR/"SRR${snakemake_wildcards['sra_acc']}" -type f ! -name "*${snakemake_wildcards['sra_acc']}*" -delete

# compress the folder with the SRA files
echo "Compressing SRA Files"
ls -laR $SRA_DIR/SRR${snakemake_wildcards['sra_acc']}/
# pigz -9 -r $SRA_DIR/${snakemake_wildcards['run']} -c > ${snakemake_output[0]}
tar -I pigz -cf ${snakemake_output[0]} $SRA_DIR/SRR${snakemake_wildcards['sra_acc']}/

# To decompress:
# tar -I pigz -xf compressed_folder.tar.gz -C /path/to/destination/directory
 