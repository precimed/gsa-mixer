#!/bin/sh
set -euo pipefail

version=20250825
wget --no-check-certificate https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_$version.zip && \
    unzip -j plink2_linux_x86_64_$version.zip && \
    rm -rf plink2_linux_x86_64_$version.zip

cp plink2 /bin