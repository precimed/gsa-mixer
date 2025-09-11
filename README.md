## Introduction

**GSA-MiXeR** is a new technique for competitive gene-set analysis, which fits a model for gene-set heritability enrichments for complex human traits, thus allowing the quantification of partitioned heritability and fold enrichment for small gene-sets.

**Cross-trait MiXeR** is a statistical tool which quantifies polygenic overlap between complex traits irrespective of genetic correlation, using GWAS summary statistics. MiXeR results are presented as a Venn diagram of unique and shared polygenic components across traits.

This repository (https://github.com/precimed/gsa-mixer) contains source code for both of these tools. 
It is only relevant to developers interested to contribute pull requests to mixer code.

User documentation is provided in a separate repository, see here: https://github.com/precimed/mixer. Reference data is shared via https://github.com/comorment/mixer . These are the only two repositories that users of either of the MiXeR tools should interact with.

If you ancounter an issue, please submit a ticket: https://github.com/precimed/mixer/issues/new . In the past my response to new tickets was very bad. If you've already submitted a ticket but didn't get a response please give me a second chance to address your issue - just push a comment and tag me @ofrei , if your question is still relevant .

Additional instructions for users at the NORMENT centre are available in https://github.com/precimed/mixer_at_norment
If the link doesn't work please reach me out by e-mail to get access.

Kind regards,
Oleksandr Frei.

## Information for developers

MiXeR software has two parts: native C++ code (``src/`` folder), and python wrapper (``precimed/`` foder).

To build native code follow the steps from ``Dockefile``. To run C++ unit tests execute ``src/build/bin/bgmg-test``.

To run python unit tests type ``py.test``.

If you have built MiXeR's native code locally, use
```
export GSA_MIXER_ROOT=$HOME/github/precimed/gsa-mixer                   # adjust accordingly
export BGMG_SHARED_LIBRARY="$GSA_MIXER_ROOT/src/build/lib/libbgmg.so"
export MIXER_PY="python $GSA_MIXER_ROOT/precimed/mixer.py"
```
and run ``MIXER_PY`` as if it was pointing to a docker or singularity container. 

## Containers

The containers are based on the following [Dockerfile](Dockerfile), built using Github actions ([this workflow](.github/workflows/docker_build_push.yml)). We also include [scripts/from_docker_image.sh](scripts/from_docker_image.sh) shell script to convert locally built Docker container into singularity, which is only relevant if you're building these containers yourself.

## Releases

To release a new version, bump ``precimed/version.py``.

In case of changes to native C++ code, update ``VERSION`` in ``src/bgmg.h`` to match ``precimed/version.py``, otherwise let the ``src/bgmg.h`` lag behind.

For the version change, see [CHANGELOG](CHANGELOG.md)

## Misc

For all other questions refer to user documentation - https://github.com/precimed/mixer .