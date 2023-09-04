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


 # Build image and deploy

#### Combine build and deploy
Might take some time to do both but for convenience since the build needs to copy the application code anyways 

``` bash
export BUCKET_NAME='orcestra-testcloudrun'
export RUN_NAME='filesystem-app'

gcloud run deploy --source . \
    --execution-environment gen2 \
    --allow-unauthenticated \
    --service-account 'orcestra@orcestra-388613.iam.gserviceaccount.com' \
    --update-env-vars BUCKET=$BUCKET_NAME

## Notable options to consider:
    --command= [COMMAND]    # entrypoint for the container image, otherwise image default Entrypoint is run
    --args=[ARG,...]        # comma-separated arguments passed to the command run by the container image
    --cpu=$CPU               # CPU limit in GKU cpu units, see gcloud run reploy --help for more details
    --min-instances=$MIN_INSTANCES
    --max-instances=$MAX_INSTANCES
    --set-env-vars=[KEY=VALUE,...]      # list of k-v pairs to set as environment variables, ALL EXISTING ENV VARIABLES WILL BE REMOVED
    --update-env-vars=[KEY=VALUE,...]   # same as above without removal
    --set-secrets=[KEY=VALUE,...]       # list of k-v pairs to set as secrets
    --update-secrets=[KEY=VALUE,...]    # ...
    --cpu-boost


```
TODO:: add in variable to set region. fine for now as the CLI asks you what region

# Something weird about how the gcloud run deploy 
I think something in my environment changes how the build works. 

#### Build Image using CI/CD
Using the Google Cloud Build tool to build on commit 
