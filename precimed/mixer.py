#!/usr/bin/env python

import logging
import sys
import argparse

from version import VERSION as __version__, MASTHEAD
import gsa_mixer.cli
import bivar_mixer.cli
import common.utils_cli

if __name__ == "__main__":
    logging.getLogger().setLevel(logging.INFO)

    parser = argparse.ArgumentParser(description=f"MiXeR v{__version__}: Univariate and Bivariate Causal Mixture for GWAS.")
    parser.add_argument('--version', action='version', version=f'MiXeR v{__version__}')

    parent_parser = argparse.ArgumentParser(add_help=False)
    common.utils_cli.parent_parser_add_argument_out_lib_log(parent_parser)

    subparsers = parser.add_subparsers()

    gsa_mixer.cli.parser_ld_add_arguments(func=gsa_mixer.cli.execute_ld_parser, parser=subparsers.add_parser("ld", parents=[parent_parser], help='prepare files with linkage disequilibrium information'))
    gsa_mixer.cli.parser_plsa_add_arguments(func=gsa_mixer.cli.execute_plsa_parser, parser=subparsers.add_parser("plsa", parents=[parent_parser], help='fit PLSA model'))
    gsa_mixer.cli.parser_split_sumstats_add_arguments(func=gsa_mixer.cli.execute_split_sumstats_parser, parser=subparsers.add_parser("split_sumstats", parents=[parent_parser], help='split summary statistics file per-chromosome'))

    bivar_mixer.cli.parser_fit_or_test_add_arguments(func=bivar_mixer.cli.execute_fit1_or_test1_parser, parser=subparsers.add_parser("fit1", parents=[parent_parser], help='fit univariate MiXeR model'), do_fit=True, num_traits=1)
    bivar_mixer.cli.parser_fit_or_test_add_arguments(func=bivar_mixer.cli.execute_fit1_or_test1_parser, parser=subparsers.add_parser("test1", parents=[parent_parser], help='test univariate MiXeR model'), do_fit=False, num_traits=1)
    bivar_mixer.cli.parser_fit_or_test_add_arguments(func=bivar_mixer.cli.execute_fit2_or_test2_parser, parser=subparsers.add_parser("fit2", parents=[parent_parser], help='fit bivariate MiXeR model'), do_fit=True, num_traits=2)
    bivar_mixer.cli.parser_fit_or_test_add_arguments(func=bivar_mixer.cli.execute_fit2_or_test2_parser, parser=subparsers.add_parser("test2", parents=[parent_parser], help='test bivariate MiXeR model'), do_fit=False, num_traits=2)

    bivar_mixer.cli.parser_perf_add_arguments(func=bivar_mixer.cli.execute_perf_parser, parser=subparsers.add_parser("perf", parents=[parent_parser], help='run performance evaluation of the MiXeR'))
    bivar_mixer.cli.parser_snps_add_arguments(func=bivar_mixer.cli.execute_snps_parser, parser=subparsers.add_parser("snps", parents=[parent_parser], help='generate random sets of SNPs'))

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
