FROM quay.io/cellgeni/starsolo:base

# Install this repository
COPY . /opt/STARsolo
RUN cd /opt/STARsolo &&\
    rm data/test.tar.gz &&\
    tar -xzf data/whitelists.tar.gz &&\
    ./install.sh

# Set PATH to include all binaries
ENV STARSOLO_WL_DIR=/opt/STARsolo/data/whitelists
ENV PATH="/opt/STARsolo/bin:${PATH}"     

# Saving Dockerfile to a file for provenance
COPY Dockerfile /docker/
RUN chmod -R 755 /docker

# Default entrypoint: run the CLI
WORKDIR /workdir
ENTRYPOINT ["starsolo"]
CMD ["--help"]