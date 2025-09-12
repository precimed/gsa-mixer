## Introduction

**GSA-MiXeR** is a new technique for competitive gene-set analysis, which fits a model for gene-set heritability enrichments for complex human traits, thus allowing the quantification of partitioned heritability and fold enrichment for small gene-sets.

**Cross-trait MiXeR** is a statistical tool which quantifies polygenic overlap between complex traits irrespective of genetic correlation, using GWAS summary statistics. MiXeR results are presented as a Venn diagram of unique and shared polygenic components across traits.

This repository (https://github.com/precimed/gsa-mixer) contains source code for both of these tools. 
It is only relevant to developers interested to contribute pull requests to mixer code.

User documentation is provided in a separate repository, see here: https://github.com/precimed/mixer. Reference data is shared via https://github.com/comorment/mixer . These are the only two repositories that users of either of the MiXeR tools should interact with.

If you ancounter an issue, please submit a ticket: https://github.com/precimed/mixer/issues/new . In the past my response to new tickets was very bad. If you've already submitted a ticket but didn't get a response please give me a second chance to address your issue - just push a comment and tag me @ofrei , if your question is still relevant .

Additional instructions for users at the NORMENT centre are available in https://github.com/precimed/mixer_at_norment
If the link doesn't work please reach me out by e-mail to get access.

For the detailed version changes, see [CHANGELOG](CHANGELOG.md)

Kind regards,
Oleksandr Frei.

## Information for developers

MiXeR software has two parts: native C++ code (``src/`` folder), and python wrapper (``precimed/`` foder).

To build native code follow the steps from ``Dockefile``. Or, check [Build from source - Linux](#build-from-source---linux) - but do this at your own risk of running into ``segmentation fault (core dumped)`` and all kind of other issues. For the end users the recommended way is to run via pre-compiled container.

To run C++ unit tests execute ``src/build/bin/bgmg-test``.

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


### Build from source - Linux

These are some legacy instructions - the exact steps depend  on your build environment. 
* If you work in HPC environment with modules system, you can load some existing combination of modules that include Boost libraries:
  ```
  module load CMake/3.15.3-GCCcore-8.3.0 Boost/1.73.0-GCCcore-8.3.0 Python/3.7.4-GCCcore-8.3.0 # TSD (gcc)
  module load Boost/1.71.0-GCC-8.3.0 Python/3.7.4-GCCcore-8.3.0 CMake/3.12.1                   # SAGA (gcc)  
  module load Boost/1.68.0-intel-2018b-Python-3.6.6 Python/3.6.6-intel-2018b CMake/3.12.1      # SAGA (intel)
  ```
* Alternatively, you may download and compile Boost libraries yourself:
  ```
  cd ~ && wget https://dl.bintray.com/boostorg/release/1.69.0/source/boost_1_69_0.tar.gz 
  tar -xzvf boost_1_69_0.tar.gz && cd boost_1_69_0
  ./bootstrap.sh --with-libraries=program_options,filesystem,system,date_time
  ./b2 --clean && ./b2 --j12 -a
  ```
* Clone and compile MiXeR repository
  ```
  cd ~ && git clone --recurse-submodules -j8 https://github.com/precimed/mixer.git
  mkdir mixer/src/build && cd mixer/src/build
  cmake .. && make bgmg -j16                                   # if you use GCC compiler
  CC=icc CXX=icpc cmake .. && make bgmg -j16                   # if you use Intel compiler
  cmake .. -DBOOST_ROOT=$HOME/boost_1_69_0 && make bgmg -j16   # if you use locally compiled boost
  ```
