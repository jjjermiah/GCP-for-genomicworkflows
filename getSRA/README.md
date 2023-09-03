# Obtaining SRA files from NCBI
Sequence Read Archive (SRA) files are a type of file format used in bioinformatics that store raw sequencing data. These files contain information about the sequence of nucleotides in a DNA or RNA sample, which can be used for various types of analysis and research. 
For the purposes of this work, the SRA files are precursors which will be converted into `fastq` files to be used as input for the analyses. 


As of Sept 1, 2023, there are a total of 27,712,956 SRA accessions according to the NIH-SRA BigQuery table (`nih-sra-datastore.sra.metadata`)

The workflows in this folder should be able to obtain any data from the SRA database.

# The first step is to retrieve the list of SRA accessions for the data you wish to use for analysis.
There are a few ways to do this
### 1) from the ncbi sra browser (details and steps explained: https://www.ncbi.nlm.nih.gov/sra/docs/srasearch/)
This link also has good insight on different ways to query the available data.

### 2) Using free tools/software:
    - using the SRA metadata information that NIH has stored in a Google BigQuery: `nih-sra-datastore.sra.metadata`
    - using open-source tools such as `pysradb` (https://github.com/saketkc/pysradb)
There are two jupyter notebooks in `scripts/`:

`scripts/getSRA_BigQuery.ipynb` is a more developed set of steps to use the bigquery API and constructed a SQL to query the BigQuery table. The results are used to create a `csv` and `json` containing the accession IDs and metadata of the SRA data interested in.

`scripts/getSRA_pysradb.ipynb` uses an open-source free tool to query the online API. I have not completed the notebook as of writing this.


# The second step is to use the metadata to obtain the SRA files 


