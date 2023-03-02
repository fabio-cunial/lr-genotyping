FROM continuumio/miniconda3
MAINTAINER Fabio Cunial


# OS
RUN apt-get -qqy update --fix-missing \
    && apt-get -qqy dist-upgrade \
    && apt-get install -y --no-install-recommends \
        zlib1g-dev \
        liblzma-dev \
        libbz2-dev \
        libdeflate-dev \
        libssl-dev \
        libcurl4-openssl-dev \
        apt-transport-https \
        ca-certificates \
        tree \
        gnupg \
        ssh \
        time \
        curl \
        wget \
        autotools-dev \
        autoconf \
        automake \
        make \
        cmake \
        gcc \
        g++ \
        git \
        bc \
        xvfb \
        python3-pip

# GSUTIL
ARG gcloud_version=405.0.1
RUN pip3 uninstall -y crcmod && pip3 install --no-cache-dir -U crcmod
RUN wget https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-cli-${gcloud_version}-linux-x86_64.tar.gz \
    && tar -xf google-cloud-cli-${gcloud_version}-linux-x86_64.tar.gz \
    && rm -f google-cloud-cli-${gcloud_version}-linux-x86_64.tar.gz \
    && yes | ./google-cloud-sdk/install.sh
ENV PATH=/google-cloud-sdk/bin:${PATH}

# HTSLIB
ARG htslib_version=1.16
RUN wget https://github.com/samtools/htslib/releases/download/${htslib_version}/htslib-${htslib_version}.tar.bz2 \
    && tar xjf htslib-${htslib_version}.tar.bz2 \
    && rm htslib-${htslib_version}.tar.bz2 \
    && cd htslib-${htslib_version} \
    && ./configure \
    && make \
    && make install \
    && cd .. \
    && rm -rf htslib-${htslib_version} \
    && bgzip --help

# SAMTOOLS
ARG samtools_version=1.16.1
RUN wget https://github.com/samtools/samtools/releases/download/${samtools_version}/samtools-${samtools_version}.tar.bz2 \
    && tar xjf samtools-${samtools_version}.tar.bz2 \
    && rm samtools-${samtools_version}.tar.bz2 \
    && cd samtools-${samtools_version} \
    && ./configure --without-curses \
    && make \
    && make install \
    && cd ${work_dir} \
    && rm -rf samtools-${samtools_version} \
    && samtools --help

# BCFTOOLS
ARG bcftools_version=1.16
RUN wget https://github.com/samtools/bcftools/releases/download/${bcftools_version}/bcftools-${bcftools_version}.tar.bz2 \
    && tar xjf bcftools-${bcftools_version}.tar.bz2 \
    && rm bcftools-${bcftools_version}.tar.bz2 \
    && cd bcftools-${bcftools_version} \
    && ./configure --without-curses \
    && make \
    && make install \
    && cd .. \
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

# CUTESV
RUN git clone https://github.com/tjiangHIT/cuteSV.git
COPY cutesv_setup.py ./cuteSV/setup.py
RUN cd cuteSV/ \
    && python setup.py install \
    && cuteSV -h

# SNIFFLES2
ARG sniffles2_version=2.0.7
RUN conda update conda \
    && conda install -c bioconda -c conda-forge -c defaults sniffles=${sniffles2_version} \
    && conda clean --packages --tempfiles --logfiles \
    && sniffles --help

# IGV-REPORTS
RUN pip install igv-reports

# IGV-SNAPHOT-AUTOMATOR
RUN wget https://github.com/stevekm/IGV-snapshot-automator.git \
    && cd IGV-snapshot-automator \
    && make install

# ADDITIONAL FILES
COPY preprocess_vcf.py .
