# ---- STAR, built in its own stage so that only the binary reaches the image ----
# STAR_REPO/STAR_REF select the GitHub source tarball (a tag, branch or commit).
# The image stays on STAR 2.7.10a_alpha_220818 so that results match earlier
# batches; the cellgeni build adds the Velocyto fixes on top of it (AGENTS.md).
FROM ubuntu:20.04 AS star-build

ARG DEBIAN_FRONTEND=noninteractive
ARG STAR_REPO=alexdobin/STAR
ARG STAR_REF=2.7.10a_alpha_220818

RUN apt-get update && apt-get -y --no-install-recommends -qq install \
        build-essential ca-certificates wget zlib1g-dev && \
    rm -rf /var/lib/apt/lists/*

# -msse4.2 is needed for --clipAdapterType CellRanger4 (STAR#1218)
RUN wget -q -O star.tar.gz "https://github.com/${STAR_REPO}/archive/${STAR_REF}.tar.gz" && \
    mkdir -p /star && tar -xzf star.tar.gz -C /star --strip-components=1 && \
    make -C /star/source STAR CXXFLAGS_SIMD="-msse4.2" -j"$(nproc)" && \
    /star/source/STAR --version

# ---- Final image ---------------------------------------------------------------
FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

ARG STAR_REPO=alexdobin/STAR
ARG STAR_REF=2.7.10a_alpha_220818
ARG samtools_version=1.21
ARG bbmap_version=38.97
ARG seqtk_version=1.4

#Install OS packages
RUN apt-get update && apt-get -y --no-install-recommends -qq install \
    wget gcc build-essential software-properties-common libz-dev \
    git libncurses5-dev libbz2-dev liblzma-dev default-jre bsdmainutils pbzip2

#Install STAR (built in the star-build stage)
COPY --from=star-build /star/source/STAR /opt/STAR/bin/STAR

#Install seqtk
RUN wget --no-check-certificate https://github.com/lh3/seqtk/archive/refs/tags/v${seqtk_version}.tar.gz && \
    tar -xzf v${seqtk_version}.tar.gz -C /opt && \
    cd /opt/seqtk-${seqtk_version} && \
    make && \
    cd / && rm v${seqtk_version}.tar.gz

#Install samtools
RUN wget https://github.com/samtools/samtools/releases/download/${samtools_version}/samtools-${samtools_version}.tar.bz2 && \
    tar -xvf samtools-${samtools_version}.tar.bz2 -C /opt && \
    cd /opt/samtools-${samtools_version} && \
    ./configure && \
    make && \
    make install && \
    cd / && rm samtools-${samtools_version}.tar.bz2  

#Install BBMap
RUN wget https://sourceforge.net/projects/bbmap/files/BBMap_${bbmap_version}.tar.gz && \
    tar -xzf BBMap_${bbmap_version}.tar.gz -C /opt && \
    cd /opt/bbmap && \
    ./stats.sh in=resources/phix174_ill.ref.fa.gz && \
    cd / && rm BBMap_${bbmap_version}.tar.gz

# Install this repository
COPY . /opt/STARsolo
RUN cd /opt/STARsolo &&\
    tar -xzf data/whitelists.tar.gz &&\
    ./install.sh

# Set PATH to include all binaries
ENV STARSOLO_WL_DIR=/opt/STARsolo/data/whitelists
ENV PATH="/opt/STARsolo/bin:${PATH}:/opt/STAR/bin:/opt/seqtk-${seqtk_version}:/opt/bbmap"     

#Saving Software Versions to a file
RUN echo "STAR version: $(STAR --version) (${STAR_REPO}@${STAR_REF})" >> versions.txt && \
    echo "samtools version: ${samtools_version}" >> versions.txt && \
    echo "BBMap version: ${bbmap_version}" >> versions.txt && \
    echo "seqtk version: ${seqtk_version}" >> versions.txt 

COPY Dockerfile /docker/
RUN chmod -R 755 /docker



# Default entrypoint: run the CLI
WORKDIR /workdir
ENTRYPOINT ["starsolo"]
CMD ["--help"]