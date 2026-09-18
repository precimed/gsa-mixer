#!/usr/bin/env bash
set -euo pipefail

version=20260914
archive="plink2_linux_x86_64_${version}.zip"
# Pinned Linux 64-bit build from https://www.cog-genomics.org/plink/2.0/
wget "https://s3.amazonaws.com/plink2-assets/alpha7/${archive}"
unzip -j "$archive"
rm "$archive"

cp plink2 /bin
plink2 --version
