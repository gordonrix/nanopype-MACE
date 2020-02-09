# \BUILD\---------------------------------------------------------------
#
#  CONTENTS      : Nanopype all-in-one Docker
#
#  DESCRIPTION   : Dockerfile
#
#  RESTRICTIONS  : none
#
#  REQUIRES      : none
#
# -----------------------------------------------------------------------
# Copyright (c) 2018-2020, Pay Giesselmann, Max Planck Institute for Molecular Genetics
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
# Written by Pay Giesselmann
# ---------------------------------------------------------------------------------

# BASE STAGE
FROM ubuntu:18.04 as base_stage

## system packages
RUN apt-get --yes update && apt-get install -y --no-install-recommends \
    wget gcc python3-dev python3-pip locales \
    ca-certificates lsb-release apt-transport-https

RUN update-ca-certificates
RUN locale-gen en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US:en
ENV LC_ALL en_US.UTF-8

## copy and configure nanopype
RUN mkdir -p /app
COPY . /app/
WORKDIR /app
RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install setuptools wheel
RUN python3 -m pip install -r requirements.txt




# BUILD STAGE
FROM base_stage as build_stage

MAINTAINER Pay Giesselmann <giesselmann@molgen.mpg.de>

## system packages
RUN apt-get --yes update && apt-get install -y --no-install-recommends wget \
    git gcc g++ ninja-build ca-certificates \
    binutils autoconf make cmake zlib1g-dev bzip2 libbz2-dev \
    liblzma-dev libncurses5-dev libcunit1-dev \
    python python3 python3-dev python3-pip

# run setup rules
RUN mkdir /build
RUN snakemake --snakefile rules/install.smk --config build_generator=Ninja --directory /build all
# HACK: fix broken symlink
RUN rm /build/src/guppy/lib/libz.so



# PACKAGE STAGE
FROM base_stage
RUN apt-get --yes update && \
apt-get install -y --no-install-recommends wget git gcc g++ \
    zlib1g-dev bzip2 libbz2-dev \
    liblzma-dev libncurses5-dev libcunit1 \
    libhdf5-100 libidn11 libgssapi-krb5-2 \
    ca-certificates \
    python python3 python3-dev python3-pip

## copy binaries from build stage
RUN mkdir -p /usr/local/bin
RUN mkdir -p /usr/local/lib/guppy
RUN mkdir -p /usr/local/data
COPY --from=build_stage /build/bin/* /usr/local/bin/
COPY --from=build_stage /build/lib/* /usr/local/lib/
COPY --from=build_stage /build/src/guppy/bin/* /usr/local/bin/
COPY --from=build_stage /build/src/guppy/lib/* /usr/local/lib/guppy/
COPY --from=build_stage /build/src/guppy/data/* /usr/local/data/
COPY --from=build_stage /usr/local/lib/python3.6/dist-packages /usr/local/lib/python3.6/dist-packages

# setup paths
RUN python3 -m pip install . --upgrade --install-option="--tools=/usr/local/bin"
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib/guppy

# create working directories
RUN mkdir -p /data/raw
RUN mkdir -p /data/import
RUN mkdir -p /data/processing
WORKDIR /data/processing
