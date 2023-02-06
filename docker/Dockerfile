FROM continuumio/miniconda3
MAINTAINER Fabio Cunial


# OS
RUN apt-get -y update \
	&& apt-get -y install git make cmake protobuf-compiler gcc g++ zlib1g-dev libcurl4-openssl-dev liblzma-dev libbz2-dev libdeflate-dev tree liblzma-dev bzip2 wget curl less \
	&& apt-get clean

# GSUTIL
ARG gcloud_version=405.0.1
RUN pip3 uninstall -y crcmod && pip3 install --no-cache-dir -U crcmod
RUN wget https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-cli-${gcloud_version}-linux-x86_64.tar.gz \
    && tar -xf google-cloud-cli-${gcloud_version}-linux-x86_64.tar.gz \
    && rm -f google-cloud-cli-${gcloud_version}-linux-x86_64.tar.gz \
    && yes | ./google-cloud-sdk/install.sh
ENV PATH=${work_dir}/google-cloud-sdk/bin:${PATH}

# BCFTOOLS
ARG bcftools_version=1.16
RUN wget https://github.com/samtools/bcftools/releases/download/${bcftools_version}/bcftools-${bcftools_version}.tar.bz2 \
    && tar xjf bcftools-${bcftools_version}.tar.bz2 \
    && rm bcftools-${bcftools_version}.tar.bz2 \
    && cd bcftools-${bcftools_version} \
    && ./configure --without-curses \
    && make \
    && make install \
    && cd ${work_dir} \
    && rm -rf bcftools-${bcftools_version} \
    && bcftools --help

# MINIMAP2
ARG minimap2_version=2.24
RUN wget https://github.com/lh3/minimap2/releases/download/v${minimap2_version}/minimap2-${minimap2_version}_x64-linux.tar.bz2 \
    && tar xjf minimap2-${minimap2_version}_x64-linux.tar.bz2 \
    && rm minimap2-${minimap2_version}_x64-linux.tar.bz2 \
    && cd minimap2-${minimap2_version}_x64-linux \
    && cp minimap2 /usr/local/bin/ \
    && minimap2 --help

# LRCALLER 0.2
RUN wget https://github.com/DecodeGenetics/LRcaller/releases/download/v0.2/LRcaller \
    && chmod +x LRcaller \
    && mv LRcaller /usr/local/bin/LRcaller \
    && LRcaller --help || echo "LRcaller returns an error code when called"
