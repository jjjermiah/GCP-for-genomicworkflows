#!/bin/bash
set -euo pipefail


# get list of files from rule for sanity:


input_aln = ${snakemake_input['aln']}
input_log = ${snakemake_input['log']}
input_sj = ${snakemake_input['sj']}
input_chim_junc = ${snakemake_input['chim_junc']}
input_unmapped = ${snakemake_input['unmapped']}

output_filteredjunctions = ${snakemake_output['filteredjunctions']}
output_GT_AG_filteredjunctions = ${snakemake_output['GT_AG_filteredjunctions']}
output_fw_filteredjunctions = ${snakemake_output['fw_filteredjunctions']}
output_sorted_bam = ${snakemake_output['sorted_bam']}
output_sorted_indexed_bam = ${snakemake_output['sorted_indexed_bam']}

# get directory of input files:
input_dir = $(dirname ${snakemake_input['aln']})
output_dir = $(dirname ${snakemake_output['filteredjunctions']})

echo "The inputs are: --starDir ${input_dir} --minLen 200 --outDir  ${output_dir}"

circRNA_finder/postProcessStarAlignment.pl --starDir ${input_dir} --minLen 200 --outDir ${output_dir}
