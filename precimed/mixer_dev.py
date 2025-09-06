#!/usr/bin/env python

import logging
import sys
import argparse

from version import VERSION as __version__, MASTHEAD
import gsa_mixer.cli
import mixer_dev.cli
import bivar_mixer.figures
import common.utils_cli

if __name__ == "__main__":
    logging.getLogger().setLevel(logging.INFO)

    parser = argparse.ArgumentParser(description=f"MiXeR v{__version__}: Univariate and Bivariate Causal Mixture for GWAS.")
    parser.add_argument('--version', action='version', version=f'MiXeR v{__version__}')

    parent_parser = argparse.ArgumentParser(add_help=False)
    common.utils_cli.parent_parser_add_argument_out_lib_log(parent_parser)
    
    parent_parser_figures = argparse.ArgumentParser(add_help=False)
    parent_parser_figures.add_argument('--argsfile', type=open, action=common.utils_cli.LoadFromFile, default=None, help="file with additional command-line arguments")
    parent_parser_figures.add_argument("--out", type=str, default="mixer", help="prefix for the output files")
    parent_parser_figures.add_argument('--ext', type=str, default=['png'], nargs='+', choices=['png', 'svg'], help="output extentions")
    parent_parser_figures.add_argument('--zmax', type=float, default=10, help="limit for z-vs-z density plots")
    parent_parser_figures.add_argument('--statistic', type=str, nargs='+', default=["point_estimate"], choices=["point_estimate", "mean", "median", "std", "min", "max"], help="Which statistic to show in the tables and on the Venn diagrams. Can have multiple values. In the case of venn diagram, the first value (typically 'point_estimate' or 'mean') indicate the size of the venn diagram; the second value (optional, typically 'std') allow to include error bars on the Venn diagramm.")
    parent_parser_figures.add_argument('--dev', default=True, action="store_true", help=argparse.SUPPRESS) 

    subparsers = parser.add_subparsers()

    gsa_mixer.cli.parser_ld_add_arguments(func=gsa_mixer.cli.execute_ld_parser, parser=subparsers.add_parser("ld", parents=[parent_parser], help='prepare files with linkage disequilibrium information'))
    gsa_mixer.cli.parser_plsa_add_arguments(func=gsa_mixer.cli.execute_plsa_parser, parser=subparsers.add_parser("plsa", parents=[parent_parser], help='fit PLSA model'))
    gsa_mixer.cli.parser_split_sumstats_add_arguments(func=gsa_mixer.cli.execute_split_sumstats_parser, parser=subparsers.add_parser("split_sumstats", parents=[parent_parser], help='split summary statistics file per-chromosome'))

    mixer_dev.cli.parser_snps_add_arguments(func=mixer_dev.cli.execute_snps_parser, parser=subparsers.add_parser("snps", parents=[parent_parser], help='generate random sets of SNPs'))
    mixer_dev.cli.parser_add_arguments(func=mixer_dev.cli.execute_fit1_or_test1_parser, parser=subparsers.add_parser("fit1", parents=[parent_parser], help='fit univariate MiXeR model'), analysis_type="fit1")
    mixer_dev.cli.parser_add_arguments(func=mixer_dev.cli.execute_fit1_or_test1_parser, parser=subparsers.add_parser("test1", parents=[parent_parser], help='test univariate MiXeR model'), analysis_type="test1")
    mixer_dev.cli.parser_add_arguments(func=mixer_dev.cli.execute_fit2_or_test2_parser, parser=subparsers.add_parser("fit2", parents=[parent_parser], help='fit bivariate MiXeR model'), analysis_type="fit2")
    mixer_dev.cli.parser_add_arguments(func=mixer_dev.cli.execute_fit2_or_test2_parser, parser=subparsers.add_parser("test2", parents=[parent_parser], help='test bivariate MiXeR model'), analysis_type="test2")
    mixer_dev.cli.parser_add_arguments(func=mixer_dev.cli.execute_save_matlab_reference_parser, parser=subparsers.add_parser("matlab", parents=[parent_parser], help='save ldsc reference'), analysis_type='matlab')
    mixer_dev.cli.parser_add_arguments(func=mixer_dev.cli.execute_perf_parser, parser=subparsers.add_parser("perf", parents=[parent_parser], help='run performance evaluation of the MiXeR'), analysis_type='perf')
    mixer_dev.cli.parser_add_arguments(func=mixer_dev.cli.execute_fixed_effects_parser, parser=subparsers.add_parser("fixed_effects", parents=[parent_parser], help='find fixed effects (beta) at the level of tagging variants'), analysis_type='fixed_effects')

    bivar_mixer.figures.parser_one_add_arguments(func=bivar_mixer.figures.execute_one_parser, parser=subparsers.add_parser("figures_one", parents=[parent_parser_figures], help='produce figures for univariate analysis'))
    bivar_mixer.figures.parser_two_add_arguments(func=bivar_mixer.figures.execute_two_parser, parser=subparsers.add_parser("figures_two", parents=[parent_parser_figures], help='produce figures for cross-trait analysis'))
    bivar_mixer.figures.parser_combine_add_arguments(func=bivar_mixer.figures.execute_combine_parser, parser=subparsers.add_parser("combine", parents=[parent_parser_figures], help='combine .json files MiXeR runs (e.g. with different --extract setting)'))

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args(sys.argv[1:])

    header = MASTHEAD
    header += f"Call: mixer_dev.py {sys.argv[1]} "
    for a in sys.argv[2:]:
        header += f'\\\n\t{a} ' if a.startswith('--') else f'{a} '

    if sys.argv[1] not in ['split_sumstats', 'figures_one', 'figures_two', 'combine']:
        # no need to rely on libbgmg just for logging
        common.utils_cli.initialize_libbgmg_logging(args, header)
    else:
        logging.info(header)

    args.func(args)
