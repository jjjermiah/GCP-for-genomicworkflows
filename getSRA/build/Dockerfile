FROM ubuntu:latest as build
ARG VERSION=current
RUN apt-get update && apt-get --quiet install --yes curl uuid-runtime && apt-get clean
RUN curl https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${VERSION}/sratoolkit.${VERSION}-centos_linux64-cloud.tar.gz | tar xz -C /
ENV PATH=/usr/local/ncbi/sra-tools/bin:${PATH}
RUN mkdir -p /root/.ncbi && \
    printf '/LIBS/IMAGE_GUID = "%s"\n' `uuidgen` > /root/.ncbi/user-settings.mkfg && \
    printf '/libs/cloud/report_instance_identity = "true"\n' >> /root/.ncbi/user-settings.mkfg \
    printf '/libs/cloud/accept_aws_charges = "true"\n/libs/cloud/accept_gcp_charges = "true"\n' >> /root/.ncbi/user-settings.mkfg

FROM python:3.11-buster
COPY --from=build /etc/ncbi /etc/ncbi
COPY --from=build /usr/local/ncbi /usr/local/ncbi
COPY --from=build /root/.ncbi /root/.ncbi
ENV PATH=/usr/local/ncbi/sra-tools/bin:${PATH}

# Install system dependencies
RUN set -e; \
    apt-get update -y && apt-get install -y \
    tini \
    lsb-release; \
    gcsFuseRepo=gcsfuse-`lsb_release -c -s`; \
    echo "deb http://packages.cloud.google.com/apt $gcsFuseRepo main" | \
    tee /etc/apt/sources.list.d/gcsfuse.list; \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | \
    apt-key add -; \
    apt-get update; \
    apt-get install -y gcsfuse pigz\
    && apt-get clean

# Set fallback mount directory
ENV MNT_DIR /mnt/gcs

# Copy local code to the container image.
ENV APP_HOME /app
WORKDIR $APP_HOME
COPY . ./

# Install production dependencies.
RUN pip install sra-downloader==1.0.6
RUN pip install Flask==2.1.0
RUN pip install gunicorn==20.1.0

# Ensure the script is executable
RUN chmod +x /app/gcsfuse_run.sh


RUN prefetch --version


# Use tini to manage zombie processes and signal forwarding
# https://github.com/krallin/tini
ENTRYPOINT ["/usr/bin/tini", "--"] 

# Pass the startup script as arguments to Tini
CMD ["/app/gcsfuse_run.sh"]