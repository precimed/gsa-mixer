# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [2.0.1] - 2024-09-18

### Added

- add [CHANGELOG.md](CHANGELOG.md) (this file) describing changes to this project
- #15 ``mixer.py --help`` and ``mixer.py --version`` flags
- #3 implement Docker and Singularity builds via github actions; update [README.md](README.md) with new installation procedure
- [scripts/process_gsa_mixer_output.py](scripts/process_gsa_mixer_output.py) for formatting GSA-MiXeR results, together with sample input and output files for this script in [out_example/](out_example/) folder

### Changed

- Replace HRC with 1kG EUR reference panel in [README.md](README.md) and [scripts/GSA_MIXER.job](scripts/GSA_MIXER.job) instructions

### Fixed

- #4 add missing ``utils_obsolete.py`` file
- #6 #7 #8 #9 #10 address incomplete or inconsistent instructions in [README.md](README.md) file

### Removed

- ``scripts/MAGMA.job`` is removed; instructions to run MAGMA are now part of [scripts/GSA_MIXER.job](scripts/GSA_MIXER.job)

## [1.0.0] - 2024-02-03

- Initial release
