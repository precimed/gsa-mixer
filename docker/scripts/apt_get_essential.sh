#!/usr/bin/env bash
set -euo pipefail

apt-get update && apt-get install -y --no-install-recommends apt-utils
apt-get update && apt-get install -y --no-install-recommends ca-certificates && \
   update-ca-certificates

# (!) Keep the list below sorted (!)
# Use available packages from the Ubuntu release selected in Dockerfile.
# Exact versions can disappear from the live Ubuntu repositories.
apt-get update && apt-get install -y --no-install-recommends \
   autoconf \
   automake \
   build-essential \
   bzip2 \
   cmake \
   curl \
   dos2unix \
   gdb \
   gfortran \
   git \
   less \
   libatlas-base-dev \
   libcurl4-openssl-dev \
   libgomp1 \
   libgsl-dev \
   libnss3 \
   libpcre2-dev \
   libxt-dev \
   pandoc \
   parallel \
   perl \
   pkg-config \
   python3 \
   python3-pytest \
   tar \
   tofrodos \
   unzip \
   vim \
   wget \
   zlib1g-dev

apt-get clean && rm -rf /var/lib/apt/lists/*
   
# /usr/bin/python must exist for bgenix, qctool
update-alternatives --install /usr/bin/python python /usr/bin/python3 10
update-alternatives --install /usr/bin/py.test py.test /usr/bin/py.test-3 10

