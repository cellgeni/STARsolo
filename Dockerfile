FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

ARG star_version=2.7.10a_alpha_220818
ARG samtools_version=1.21
ARG bbmap_version=38.97
ARG seqtk_version=1.4

#Install OS packages
RUN apt-get update && apt-get -y --no-install-recommends -qq install \
    wget gcc build-essential software-properties-common libz-dev \
    git libncurses5-dev libbz2-dev liblzma-dev default-jre bsdmainutils

#Install STAR
RUN wget --no-check-certificate https://github.com/alexdobin/STAR/archive/${star_version}.tar.gz && \
    tar -xzf ${star_version}.tar.gz -C /opt && \
    cd /opt/STAR-${star_version}/source && \
    make STAR CXXFLAGS_SIMD="-msse4.2" && \
    cd / && rm ${star_version}.tar.gz 

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
    rm data/test.tar.gz &&\
    tar -xzf data/whitelists.tar.gz &&\
    ./install.sh

# Set PATH to include all binaries
ENV STARSOLO_WL_DIR=/opt/STARsolo/data/whitelists
ENV PATH="/opt/STARsolo/bin:${PATH}:/opt/STAR-${star_version}/source:/opt/seqtk-${seqtk_version}:/opt/bbmap"     

#Saving Software Versions to a file
RUN echo "STAR version: ${star_version}" >> versions.txt && \
    echo "samtools version: ${samtools_version}" >> versions.txt && \
    echo "BBMap version: ${bbmap_version}" >> versions.txt && \
    echo "seqtk version: ${seqtk_version}" >> versions.txt 

COPY Dockerfile /docker/
RUN chmod -R 755 /docker

# Default entrypoint: run the CLI
ENTRYPOINT ["starsolo"]
CMD ["--help"]