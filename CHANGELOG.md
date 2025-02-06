# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased] - https://github.com/precimed/gsa-mixer

### Added

* ``scripts/MIXER.job`` script with instructions reproducing univariate and bivariate MiXeR v1.3 analyses.

## [2.1.1] - 2025-02-05 - https://github.com/precimed/gsa-mixer

### Changed

* reduce search space in 'mixer.py fit1 diffevo-fast'

## [2.1.0] - 2024-09-25 - https://github.com/precimed/gsa-mixer

This version includes various fixes to support ``mixer.py [fit1,test1,fit2,test2]`` (univariate and cross-trait) analyses, plus a few new experimental options.

### Added

* new ``--kmax-pdf`` option allows to customize the number of sampling iteration used for QQ plot and power curves. In other scenarios (e.g. log-likelihood computation) the number of sampling iterations is controlled by ``--kmax`` parameter.
* ``--nckoef`` option allows user to specify coefficient between 0 and 1 for translating total number of causal variants to the number of causal variants explaining 90% of trait's heritability. By default this is now estimated from the data. Up until MiXeR v1.3 this was set to 0.319, a value specific to 1kG reference panel used for analysis.

### Changed

* The default settings for weighting SNPs in log-likelihood function in univariate (``fit1``, ``test1``) and cross-trait (``fit2``, ``test2``) analyses had changed from random prunning (``--randprune-r2 0.1 --randprune-n 64``) to weighting by inverse residual LD score, for consistency with GSA-MiXeR (``mixer.py plsa --gsa-base`` and ``--gsa-full``). In the cases where a maximum consistency with previous versions of univariate and cross-trait analyses is required, the following set of switches will configure SNP weights to work equivalently to the previous version: ``--disable-inverse-ld-score-weights --randprune-r2 0.1 --randprune-n 64 --hardprune-maf 0``. This should be added on top of ``--extract 1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.$REP.snps`` for ``fit1`` and ``fit2`` analyses. For new analyses it's recommended to keep the new weighting scheme, which is now the default setting.
* The inference scheme had a slight change due to new parametrization of the ``sig2_beta``; the old ``sig2_beta`` (an actual variance of causal effects) is now represented by the ratio of the ``sig2_beta`` over ``pi``. The reason for this re-parametrization is to facilitate inference of the parameters, so that trait's SNP-based heritability depends on the new ``sig2_beta``, but not on ``pi`` parameter. This implies a subtle change in the fit procedure (namely in ``diffevo`` and ``diffevo-fast`` steps of the ``fit1``), where the range of ``sig2_beta`` parameter is calibrated to search for trait's heritability between ``1e-5 <= h2 <= 10.0``.
* Old ``<out>.json`` files produced by previous versions of the MiXeR software will no longer work with ``MiXeR v2.0.0`` and later versions, as ``.json`` serialization includes many internal changes (e.g. ``UnivariateParams`` class was replaced with``AnnotUnivariateParams`` class).

### Fixed

* ``mixer_figures.py combine`` analysis applied to the output of ``mixer_figures.py test1`` will now also combine ``qqplot_bins`` sections of the ``json``; after this fix the subsequent ``mixer_figures.py one`` command will also produced a binned 3x3 figure with QQ plots
* [CHANGELOG.md](CHANGELOG.md) (this file) was greatly improed to cover full development history of the MiXeR and GSA-MiXeR packages.

### Experimental options (limited support)

* The ``mixer.py fit1`` now supports ``--load-params-file <out>.json`` option which can be used with ``<out>.json`` files produced by ``mixer.py plsa``. This is examplified by [PLSA_MIXER_UNIVAR.job](scripts/PLSA_MIXER_UNIVAR.job) script (and, subsequencly, [PLSA_MIXER_BIVAR.job](scripts/PLSA_MIXER_BIVAR.job) script for bivariate analysis). For this to work several new options (namely ``--annot-file``, ``--go-file``, ``--go-extend-bp``) were added to ``fit1``, ``fit2``, ``test1`` and ``test2`` analyses. These arguments must be kept consistent with arguments used in ``mixer.py plsa`` analysis used to produce the univariate model passed via ``--load-params-file``. The role of ``fit1`` (when used with ``--load-params-file``) is to estimate trait's polygenicity (which is only partly supported by ``mixer.py plsa`` analysis).
* New ``mixer.py fit1 --fit-sequence infinitesimal`` and ``mixer.py fit2 --fit-sequence infinitesimal`` allow to fit infinitesimal univariate and bivariate models.
* ``rg_sig2_factor`` output is included in the ``.csv`` tables from cross-trait analysis
* ``tag_ez2`` column included in ``.snps.csv`` tables, enabled by ``--make-snps-file``

## [2.0.1] - 2024-09-18 - https://github.com/precimed/gsa-mixer

### Added

* add [CHANGELOG.md](CHANGELOG.md) (this file) describing changes to this project
* #15 ``mixer.py --help`` and ``mixer.py --version`` flags
* #3 implement Docker and Singularity builds via github actions; update [README.md](README.md) with new installation procedure
* [scripts/process_gsa_mixer_output.py](scripts/process_gsa_mixer_output.py) for formatting GSA-MiXeR results, together with sample input and output files for this script in [out_example/](out_example/) folder

### Changed

* Replace HRC with 1kG EUR reference panel in [README.md](README.md) and [scripts/GSA_MIXER.job](scripts/GSA_MIXER.job) instructions

### Fixed

* #4 add missing ``utils_obsolete.py`` file
* #6 #7 #8 #9 #10 address incomplete or inconsistent instructions in [README.md](README.md) file

### Removed

* ``scripts/MAGMA.job`` is removed; instructions to run MAGMA are now part of [scripts/GSA_MIXER.job](scripts/GSA_MIXER.job)

## [1.0.0] - 2024-02-03 - https://github.com/precimed/gsa-mixer

Thsi is an initial release of GSA-MiXeR software, including the following changes as compared to non-GSA MiXeR v1.3 from https://github.com/precimed/mixer.

### Added

* New ``mixer.py plsa`` analysis allows fitting flexible univariate models (as described in [this](https://www.nature.com/articles/s41588-024-01771-1) paper), with model fit based on Adam algorithm (first-order gradient-based stochastic optimization). ``mixer.py plsa`` analysis has two specific configurations, triggered with ``mixer.py plsa --gsa-base`` and ``mixer.py plsa --gsa-full`` flags respectively, allowing to perform GSA-MiXeR baseline and full analses.
* New ``mixer.py split_sumstats`` option allows to split summary statistics into one file per chromosome (a required format for 
``mixer.py plsa`` analysis)
* new ``mixer.py plsa --load-baseline-params-file`` option allows to specify baseline model to use for fold enrichments (if requested by ``--go-file-test`` and/or ``--annot-file-test`` options); this is optional, if not specified enrichment will be calculated against model specified by ``--load-params-file``; when neither ``--load-params--file`` nor ``--load-baseline-params-file`` fold enrichments can not be calculated, and using ``[--go,--annot]-file-test`` will result in an error
* added support for ``effectallele``, ``otherallele``, ``rsid`` and ``pos`` column names for summary statistics, in addition to previous column names (``a1``, ``a2``, ``snp`` and ``bp``)
* the following options allow for a more flexible selection of SNP used for analysis, and for configuring their weights in log-likelihood function:
  * new ``--exclude-ranges`` option is added to all ``mixer.py`` analyses, allowing to exclude genomic regions based on CHR and BP (such as for example MHC or APOE)
  * new ``--allow-ambiguous-snps`` option allows to retain ambiguous SNPs for analysis (by default MiXeR still excludes them)
  * new ``--hardprune-maf``, ``--hardprune-r2``,  and ``--hardprune-subset`` options allow to exclude SNPs from fit procedure, offering a more convenient alternative to the ``--extract`` option used in the previous MiXeR version. ``--hardprune-subset`` selects random subset of SNPs of a certain size (as fraction of ``--bim-file`` SNPs when parameter is a float-point number between 0 and 1;, or as an absolute number of SNPs if parameter is an integer number above 1; this is done in a deterministic way, randomized based on the ``--seed`` argument); ``--hardprune-maf`` filters SNPs on minor allele frequency (allele frequencies are based on ``--ld-file``, i.e. from the genotype panel used to generate the LD reference); ``--hardprune-r2`` is used to randomly prune SNPs in high LD.
  * due to new default weighting scheme by inverse residual LD score a new flag ``--disable-inverse-ld-score-weights`` was added, allowing to disable inverse LD score weighting and thus restore previous behavior with weights set by random prunning; *residual LD score* means that LD score of each tag SNP is computed towards SNPs that have a well-defined z-score, and pass other filtering (``--extract``, ``--hardprune-maf``, etc).
  * new ``--save-weights`` option allows to save ``<out>.weights`` file with per-SNP weights generated using ``hardprune-`` and ``randprune-`` options. By default weights are saved for ``plsa`` analysis.
  * ``--weights-file`` loads weights  from a previously generated ``<out>.weights`` file. If weights are loaded, ``--save-weights`` argument is ignored. The default for ``plsa`` and ``fit1`` analysis is to use ``--weights-file auto`` which acts depending on ``--load-params-file``: if a previous params file is used, then the respective weights will be loaded. To overwrite ``--weights-file auto`` behavior set by default use ``--weights-file none``, which will disable loading weights. Whenever weights are not loaded, they are generated according to hardprune/randprune settings.
* ``--save-ldsc-reference`` now also exports the data in MATLAB format
* new ``--make-snps-file`` option (available for ``plsa``, ``fit1``, ``test1``, ``fit2``, ``test2`` analyses) allows to output per-SNP posterior effect size estimates
* new ``--use-complete-tag-indices``, ``--loadlib-file`` and ``--savelib-file`` options allows for faster loading of LD matrices
* new unit-tests are implemented (namely ``precimed/mixer-test/test_cli.py``)

### Changed

* the content of an existing ``<out>.log`` file is cleared prior to running ``mixer.py``

### Fixed

* Non-digit chr label codes should be ignored instead of crashing
* ``--help`` output is greatly improved for all ``mixer.py`` analyses and arguments

### Experimental options (limited support)

* ``mixer.py plsa --fit`` option allows to customize which parameters of the model should be fitted from the data. The choices for ``--fit`` are as follows: ``sig2-zeroL``, ``sig2-zeroA``, ``s``, ``l``, ``pi``, ``annot``, ``gene``. Multiple values can be provided to jointly fit several parameters. The ``--fit`` argument is required for ``mixer.py plsa`` analysis, except if ``--gsa-base`` and ``--gsa-full`` flags were used, in which case ``--fit`` argument will be ignroed, and fit sequence will be configured as needed for the GSA-MiXeR model. Other parameters, not listed in ``--fit``, will be constrained either based on ``--load-params-file``, or through ``--s-value``, ``--l-value``, ``--pi-value``, ``--sig2-zeroA`` and ``--sig2-zeroL`` values.

## [1.3.0] - 2020-08-31 - https://github.com/precimed/mixer

A comprehensive update to usage instructions of univariate and cross-trait MiXeR.
Confusingly, the main vehicle for releasing this version of MiXeR were singularity containers, shared separately in https://github.com/comorment/mixer repository, alongisde with usage instructions on [synthetic](https://github.com/comorment/mixer/blob/main/usecases/mixer_simu.md) and [real](uhttps://github.com/comorment/mixer/blob/main/usecases/mixer_real.md) data.
Further complicating things, the MiXeR software internally never had a v1.3.0 tag. If you're looking into the log files, they will still indicate ``MiXeR v1.2.0``. This is because the MiXeR software (python and C++ codes) used to perform ``MiXeR v1.3.0`` analysis was nearly identical to the codes used for ``v1.2.0``, however the README file with  instructions on how to run the software has changed, and a few new reference files were generated.

### Added

* 20 new reference files ``1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.repNN.snps`` (intended for ``--extract`` argument).
* ``mixer.py snps`` analyses allowing to generate a custom set of ``.snps`` files for the ``--extract`` argument.
* New example scripts are provided in [scripts](https://github.com/precimed/mixer/tree/master/scripts) folder, illustrating how to use a job array (with ``#SBATCH --array=1-20``) to execute MiXeR 20 times, and use ``--extract 1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep${SLURM_ARRAY_TASK_ID}.snps`` flag in ``mixer.py fit1`` and ``mixer.py fit2`` analyses, and also adjusting input/output file names accordingly.
* New ``mixer_figures.py combine`` option aggregates ``.json`` files produced by MiXeR, averaging individual parameter estimates, and calculating standard deviations.
* New option ``--statistic mean std`` is implemented for ``mixer_figures.py one`` and ``mixer_figures.py two`` commands. This option must be specified to generate figures based on "combined" ``.json`` produced by ``mixer_figures.py combine``.

### Changed

* If you previously generated your input summary statistics constraining them to HapMap3 SNPs, you need to re-generated them without constraining to HapMap3.
* With MiXeR v1.3 you should expect a slightly lower polygenicity estimate, mainly because of ``--maf 0.05`` constraint on the frequency of genetic variants used in the fit procedure. The rationale for the ``--maf 0.05`` recommendation is to still have a robust selection of SNPs in the fit procedure despite not having HapMap3 filter.
* With MiXeR v1.3 you should expect a 10 fold increase in CPU resources needed per ran (20 times more runs, but each run is ~600K SNPs which is half the size of HapMap3).


### Fixed

* In earlier version MiXeR estimates in bivariate analysis may, in some cases, depend on which of the two traits is selected as ``--trait1``.  Therefore an earlier recommendation was to run MiXeR twice, swapping the traits. This is no longer needed due to the averaging of the results across 20 iterations. Earlier the reason for instability was due to low power in GWAS, and AIC/BIC were able capture such an issue. In case of significant AIC and/or BIC values running MiXeR twice (with/without swapping the traits) hasn't ever been necessary.

## [1.2.0] - 2020-03-21 - https://github.com/precimed/mixer

### Added

* New ``mixer.py ld`` analysis generates files for ``--ld-file`` argument based on a custom genotype panel.
* ``mixer_figures.py`` is expanded to generate more data, including AIC and BIC values, a Dice coefficient for bivariate analysis, and negative log-likelihood plots.
* New optional flags ``--z1max`` and ``--z2max`` allow to specify thresholds for right-censoring.
* Added support for Intel compiler, with build instructions available for SAGA.
* Added ``mixer.py perf`` analysis to evaluate performance of MiXeR cost funciton (in terms or runtime / memory usage).


### Changed

* You are required to update your scripts using new reference files which are available on NIRD ``<SUMSTAT>/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld``. There is no need to update input files with summary statistics.
* ``mixer.py fit`` analysis is split into four commands: ``fit1``, ``test1`` (univariate analysis) and ``fit2``, ``test2`` (bivariate analysis).
  ``fit`` commands are time-consuming and estimate MiXeR model parameters, typically from a subset of SNPs such as HapMap3 (can be passed with ``--extract`` flag).
  ``test`` commands are much faster and estimate power plots, QQ-plots, and similar, typically on a full set of SNPs.
  All commands have custom defaults, so you no longer have to specify ``--fit-sequence``, ``--kmax``, ``--qq-plots``, ``--power-curves``, etc.
* MiXeR v1.2 uses precisely the same fit procedure as MiXeR v1.1, however it has a certain changes in how the LD matrix is stored internally, and also in how the cost function is evaluated. Particularly, the bivariate cost function based on sampling approach is now stateless, and it now consuming substentially less memory. The univariate cost function is based on convolution approach, as in v1.1, however it is now coupled with sampling for large z-scores to allow for ``--z1max`` and ``--z2max`` parameters.

### Removed

* MATLAB code is removed, from now you should run MiXeR via python interface.
* MiXeR no longer provides standard errors on bivariate parameter estimates.
  Univariate parameter estimates can still, in principle, be calculated with ``--ci-alpha 0.05`` option, but this is turned off default, and is not recommended.
  Inference about statistical significance of the parameters should be made in context of the model selection criterion: Akaike (AIC) or Bayesian (BIC).
* ``--plink-ld-bin0`` and ``--frq-file`` arguments are removed. Use ``--ld-file`` argument instead. 

## [1.1.0] - 2020-04-20 - https://github.com/precimed/mixer

Internal release of univariate and cross-trait MiXeR for Matlab and Python

## [1.0.0] - 2019-05-23 - https://github.com/precimed/mixer

 Initial public release of univariate and cross-trait MiXeR as Matlab software
