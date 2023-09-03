
# Building the dockerfiles in this repo

``` bash
export DOCKER_REPO="jjjermiah"
export DOCKER_IMAGE_NAME="pigz"
export DOCKER_TAG="0.9"
export DOCKERFILE_NAME="Dockerfile.${DOCKER_IMAGE_NAME}"
docker build -t $DOCKER_REPO/$DOCKER_IMAGE_NAME:$DOCKER_TAG -f $DOCKERFILE_NAME .
```


# Pushing dockerfiles to dockerhub

``` bash
export 
docker push $DOCKER_REPO/$DOCKER_IMAGE_NAME:$DOCKER_TAG
```



# Pushing to artifact registry 
``` bash
docker tag jjjermiah/sratools:0.2 northamerica-northeast2-docker.pkg.dev/orcestra-388613/bhklab-docker-repo/sratools
docker push northamerica-northeast2-docker.pkg.dev/orcestra-388613/bhklab-docker-repo/sratools

```