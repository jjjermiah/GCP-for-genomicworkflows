#!/bin/bash
mkdir .ncbi
cp /root/.ncbi/user-settings.mkfg .ncbi/user-settings.mkfg

vdb-config --root --set /repository/remote/disabled="false"
vdb-config --root --set /repository/remote/main/CGI/resolver-cgi="https://trace.ncbi.nlm.nih.gov/Traces/names/names.fcgi"
vdb-config --root --set /repository/remote/main/SDL.2/resolver-cgi="https://locate.ncbi.nlm.nih.gov/sdl/2/retrieve"
vdb-config --root --set /repository/remote/protected/CGI/resolver-cgi="https://trace.ncbi.nlm.nih.gov/Traces/names/names.fcgi"
vdb-config --root --set /repository/remote/protected/SDL.2/resolver-cgi="https://locate.ncbi.nlm.nih.gov/sdl/2/retrieve"

/usr/local/bin/align-info "SRR${snakemake_wildcards['sra_acc']}" | cut -d ',' -f1 >${snakemake_output[0]}
# "/usr/local/bin/align-info {wildcards.sra_acc} | cut -d ',' -f1 >${snakemake_output[0]}"
