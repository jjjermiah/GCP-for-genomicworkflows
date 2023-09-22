#!/bin/bash

# Read the file line by line
while IFS= read -r path; do
  # Run the gcloud command for each path
  echo " gcloud storage objects update '$path' --add-acl-grant=entity=AllUsers,role=READER"
  gcloud storage objects update "$path" --add-acl-grant=entity=AllUsers,role=READER
done < "CCLE_data_paths.txt"