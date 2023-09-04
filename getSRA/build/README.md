# Build docker image with gcsfuse

``` bash
export REGION='us-central1'
export BUCKET_NAME='orcestra-testcloudrun'

gsutil mb -l $REGION gs://$BUCKET_NAME
```


``` bash
export FS_SERVICE_ACC='fs-identity'
gcloud iam service-accounts create $FS_SERVICE_ACC
```


``` bash
export PROJECT_ID='orcestra-388613'

gcloud projects add-iam-policy-binding $PROJECT_ID \
     --member "serviceAccount:fs-identity@$PROJECT_ID.iam.gserviceaccount.com" \
     --role "roles/storage.objectAdmin" \
     --role "roles/storage.objects.list" 
```

``` bash
rm Dockerfile
cp Dockerfile.gcsfuse Dockerfile
```

``` bash
export RUN_NAME='filesystem-app'

gcloud run deploy $RUN_NAME --source . \
    --execution-environment gen2 \
    --allow-unauthenticated \
    --service-account 'orcestra@orcestra-388613.iam.gserviceaccount.com' \
    --update-env-vars BUCKET=$BUCKET_NAME
```

