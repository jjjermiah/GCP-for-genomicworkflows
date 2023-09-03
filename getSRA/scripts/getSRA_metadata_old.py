import os
from google.cloud import storage
import pysradb

# The goal of this script is to get the metadata for each SRA file in the bucket
# and add it to the metadata of the SRA file in the bucket
# This was before reworking the pipeline to get metadata before obtaining the SRA files

# Keeping this script for reference
# Dont think it works right now anyways.


# Set the bucket name and path to SRA folders
BUCKET_NAME = "bucket-name"
PATH_TO_SRA_FOLDERS = "path/to/SRA/folders"

# Create a client to interact with the storage and SRAweb api 
storage_client = storage.Client()
db = pysradb.sraweb.SRAweb()

# Get a list of SRA files in the bucket folder
blobs = storage_client.list_blobs(BUCKET_NAME, prefix=PATH_TO_SRA_FOLDERS)

srafolders = dict()
for blob in blobs:
    if blob.name.endswith(".sra"):
        x = blob.name.split('/')[-1]
        print(x)
        metadata = db.sra_metadata(x.split('.')[0]).iloc[0].to_dict()
        print(metadata)


blob.metadata = pysradb_metadata

def set_blob_metadata(bucket_name, blob_name):
    """Set a blob's metadata."""
    # bucket_name = 'your-bucket-name'
    # blob_name = 'your-object-name'

    storage_client = storage.Client()
    bucket = storage_client.bucket(bucket_name)
    blob = bucket.get_blob(blob_name)
    metageneration_match_precondition = None

    # Optional: set a metageneration-match precondition to avoid potential race
    # conditions and data corruptions. The request to patch is aborted if the
    # object's metageneration does not match your precondition.
    metageneration_match_precondition = blob.metageneration

    metadata = {'color': 'Red', 'name': 'Test'}
    blob.metadata = metadata
    blob.patch(if_metageneration_match=metageneration_match_precondition)

    print(f"The metadata for the blob {blob.name} is {blob.metadata}")



# Extract the SRA accessions from the blob names
# sra_accessions = [os.path.basename(blob.name) for blob in blobs]
# sra_accessions = [os.path.basename(blob.name) for blob in blobs if blob.name.endswith('.sra')]
sra_files = [blob.name for blob in blobs if blob.name.endswith(".sra")]

dirnames, sra_files = zip(*[os.path.split(sra_file) for sra_file in sra_files])
sra_files = list(set(sra_files))

# sort the unique SRA accessions
sra_files.sort()

# Write the unique SRA accessions to a file
with open("SRA_accessions_unique.txt", "w") as file:
    file.write("\n".join(sra_files))

# Get the length of the unique SRA accessions
# num_unique_accessions = len(sra_files)
# print(f"Number of unique .sra files : {num_unique_accessions}")

# pysradb_metadata = db.sra_metadata([sra_file.split('.')[0] for sra_file in sra_files])
