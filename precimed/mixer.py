#!/usr/bin/env python

import logging
import sys
import argparse

from version import VERSION as __version__
import gsa_mixer.cli
import common.utils_cli

MASTHEAD = "***********************************************************************\n"
MASTHEAD += "* (c) 2016-2025 MiXeR software\n"
MASTHEAD += "* Norwegian Centre for Mental Disorders Research / University of Oslo\n"
MASTHEAD += "* Center for Multimodal Imaging and Genetics / UCSD\n"
MASTHEAD += "* GNU General Public License v3\n"
MASTHEAD += "***********************************************************************\n"

if __name__ == "__main__":
    logging.getLogger().setLevel(logging.INFO)

    parser = argparse.ArgumentParser(description=f"MiXeR v{__version__}: Univariate and Bivariate Causal Mixture for GWAS.")
    parser.add_argument('--version', action='version', version=f'MiXeR v{__version__}')

    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--argsfile', type=open, action=common.utils_cli.LoadFromFile, default=None, help=argparse.SUPPRESS) #"file with additional command-line arguments, e.g. those that are across all of your mixer.py runs (--lib, --bim-file and --ld-file)")
    parent_parser.add_argument("--out", type=str, default="mixer", help="prefix for the output files, such as <out>.json (default: %(default)s);")
    parent_parser.add_argument("--lib", type=str, default="libbgmg.so", help="path to libbgmg.so plugin (default: %(default)s); can be also specified via BGMG_SHARED_LIBRARY env variable.")
    parent_parser.add_argument("--log", type=str, default=None, help="file to output log (default: <out>.log); NB! if --log points to an existing file the new lines will be appended to it at the end of the file.")

    subparsers = parser.add_subparsers()

    if sys.argv[1] in ['ld', 'snps', 'plsa', 'split_sumstats']:
        gsa_mixer.cli.parser_ld_add_arguments(func=gsa_mixer.cli.execute_ld_parser, parser=subparsers.add_parser("ld", parents=[parent_parser], help='prepare files with linkage disequilibrium information'))
        gsa_mixer.cli.parser_snps_add_arguments(func=gsa_mixer.cli.execute_snps_parser, parser=subparsers.add_parser("snps", parents=[parent_parser], help='generate random sets of SNPs'))
        gsa_mixer.cli.parser_plsa_add_arguments(func=gsa_mixer.cli.execute_plsa_parser, parser=subparsers.add_parser("plsa", parents=[parent_parser], help='fit PLSA model'))
        gsa_mixer.cli.parser_split_sumstats_add_arguments(func=gsa_mixer.cli.execute_split_sumstats_parser, parser=subparsers.add_parser("split_sumstats", parents=[parent_parser], help='split summary statistics file per-chromosome'))
    elif sys.argv[1] in ['fit1', 'test1', 'fit2', 'test2', 'matlab', 'perf', 'fixed_effects']:
        raise NotImplementedError()
        bivar_mixer.cli.parser_add_arguments(func=bivar_mixer.cli.execute_fit1_or_test1_parser, parser=subparsers.add_parser("fit1", parents=[parent_parser], help='fit univariate MiXeR model'), analysis_type="fit1")
        bivar_mixer.cli.parser_add_arguments(func=bivar_mixer.cli.execute_fit1_or_test1_parser, parser=subparsers.add_parser("test1", parents=[parent_parser], help='test univariate MiXeR model'), analysis_type="test1")
        bivar_mixer.cli.parser_add_arguments(func=bivar_mixer.cli.execute_fit2_or_test2_parser, parser=subparsers.add_parser("fit2", parents=[parent_parser], help='fit bivariate MiXeR model'), analysis_type="fit2")
        bivar_mixer.cli.parser_add_arguments(func=bivar_mixer.cli.execute_fit2_or_test2_parser, parser=subparsers.add_parser("test2", parents=[parent_parser], help='test bivariate MiXeR model'), analysis_type="test2")
        bivar_mixer.cli.parser_add_arguments(func=bivar_mixer.cli.execute_save_matlab_reference_parser, parser=subparsers.add_parser("matlab", parents=[parent_parser], help='save ldsc reference'), analysis_type='matlab')
        bivar_mixer.cli.parser_add_arguments(func=bivar_mixer.cli.execute_perf_parser, parser=subparsers.add_parser("perf", parents=[parent_parser], help='run performance evaluation of the MiXeR'), analysis_type='perf')
        bivar_mixer.cli.parser_add_arguments(func=bivar_mixer.cli.execute_fixed_effects_parser, parser=subparsers.add_parser("fixed_effects", parents=[parent_parser], help='find fixed effects (beta) at the level of tagging variants'), analysis_type='fixed_effects')

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args(sys.argv[1:])

    header = MASTHEAD
    header += f"Call: mixer.py {sys.argv[1]} "
    for a in sys.argv[2:]:
        header += f'\\\n\t{a} ' if a.startswith('--') else f'{a} '
    logging.info(header)

    if sys.argv[1] not in ['split_sumstats']:
        # no need to rely on libbgmg just for logging
        common.utils_cli.initialize_libbgmg_logging(args, header)

    args.func(args)
