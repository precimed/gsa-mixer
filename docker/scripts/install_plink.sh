#!/bin/sh
set -euo pipefail

# plink
version=20250819
wget --no-check-certificate https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_$version.zip && \
    unzip -j plink_linux_x86_64_$version.zip && \
    rm -rf plink_linux_x86_64_$version.zip
cp plink /bin