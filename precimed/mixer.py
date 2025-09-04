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

    if sys.argv[1] in ['ld', 'plsa', 'split_sumstats']:
        parser = argparse.ArgumentParser(description=f"MiXeR v{__version__}: Univariate and Bivariate Causal Mixture for GWAS.")
        parser.add_argument('--version', action='version', version=f'MiXeR v{__version__}')

        parent_parser = argparse.ArgumentParser(add_help=False)
        common.utils_cli.parent_parser_add_argument_out_lib_log(parent_parser)

        subparsers = parser.add_subparsers()

        gsa_mixer.cli.parser_ld_add_arguments(func=gsa_mixer.cli.execute_ld_parser, parser=subparsers.add_parser("ld", parents=[parent_parser], help='prepare files with linkage disequilibrium information'))
        gsa_mixer.cli.parser_plsa_add_arguments(func=gsa_mixer.cli.execute_plsa_parser, parser=subparsers.add_parser("plsa", parents=[parent_parser], help='fit PLSA model'))
        gsa_mixer.cli.parser_split_sumstats_add_arguments(func=gsa_mixer.cli.execute_split_sumstats_parser, parser=subparsers.add_parser("split_sumstats", parents=[parent_parser], help='split summary statistics file per-chromosome'))
    elif sys.argv[1] in ['fit1', 'test1', 'fit2', 'test2', 'perf', 'snps']:
        parser = bivar_mixer.cli.generate_args_parser(__version__)

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
