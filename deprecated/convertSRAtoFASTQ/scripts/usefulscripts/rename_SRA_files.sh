#!/bin/bash

# This script will rename the SRA files to the format: SRR<accession>.sra
# This script will also create a file called SRA_accessions.txt which contains the SRA accessions for each sample
# It will first create a list of all the SRA file in a given directory 

# In this example, we use a GCS bucket
# The SRA files are in the directory: ${PATH_TO_SRA_FOLDERS}/{SRR_Accession}/{SRR_Accession}
# The .vdbcache files also need to be renamed if they exist and they are in the directory: ${PATH_TO_SRA_FOLDERS}/{SRR_Accession}/{SRR_Accession}.vdbcache

PATH_TO_SRA_FOLDERS="gs://orcestra-archive/rawdata/RNA/SRA"

# use gsutil to list all the SRA files in the directory
gsutil ls ${PATH_TO_SRA_FOLDERS}/*/* > SRA_files.txt


# use awk to extract the SRR accession from the file path
# Assuming that the paths in the text file look like:
# ${PATH_TO_SRA_FOLDERS}/SRR8618310/SRR8618310
# we want to extract the SRR8618310 part
# we can do this by splitting the string by the '/' character and then taking the 7th element (note: // is counted as 2 elements)
awk -F'/' '{print $7}' SRA_files.txt > SRA_accessions.txt

# use a for loop to do the following:
# for every SRA_accession in the SRA_accessions.txt file
# if there exists a line in the SRA_files.txt file that is exactly:
# ${PATH_TO_SRA_FOLDERS}/{SRR_Accession}/{SRR_Accession}
# OR 
# ${PATH_TO_SRA_FOLDERS}/{SRR_Accession}/{SRR_Accession}.vdbcache
# then add a command that renames that file to a jobfile

# create an empty list called used_SRA
touch used_SRA.txt


for SRA_accession in $(cat SRA_accessions.txt)
do
    # check if the SRA_accession has already been used
    if cat used_SRA.txt | grep -Fx "${SRA_accession}"
    then
        continue
    fi
    # if it hasn't been used, then add it to the list
    echo ${SRA_accession} >> used_SRA.txt
    if cat  SRA_files.txt | grep -Fx "${PATH_TO_SRA_FOLDERS}/${SRA_accession}/${SRA_accession}"
    then
        echo "gsutil mv ${PATH_TO_SRA_FOLDERS}/${SRA_accession}/${SRA_accession} ${PATH_TO_SRA_FOLDERS}/${SRA_accession}/${SRA_accession}.sra" >> rename_SRA_jobfile.sh
    fi
    if cat  SRA_files.txt | grep -Fx "${PATH_TO_SRA_FOLDERS}/${SRA_accession}/${SRA_accession}.vdbcache" 
    then
        echo "gsutil mv ${PATH_TO_SRA_FOLDERS}/${SRA_accession}/${SRA_accession}.vdbcache ${PATH_TO_SRA_FOLDERS}/${SRA_accession}/${SRA_accession}.sra.vdbcache" >> rename_SRA_jobfile.sh
    fi
    echo "Done with ${SRA_accession}"
done

rm SRA_files.txt
rm SRA_accessions.txt
rm used_SRA.txt

# now there will be a job file called rename_SRA_jobfile.sh
# you can run this job anywhere, either normally, or parallel 
# I use parallel -j #numcores < rename_SRA_jobfile.sh


# then rename the file to SRR{SRA_accession}.sra
# if there exists a file called {SRA_accession}.vdbcache
# then rename the file to {SRA_accession}.sra.vdbcache
