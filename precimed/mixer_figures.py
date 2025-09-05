#!/usr/bin/env python
'''
(c) 2016-2020 Oleksandr Frei, Alexey A. Shadrin, Dominic Holland
MiXeR software: Univariate and Bivariate Causal Mixture for GWAS
'''

import logging
import sys

import argparse

import common.utils_cli
import bivar_mixer.figures
from version import VERSION as __version__, MASTHEAD

if __name__ == "__main__":
    logging.getLogger().setLevel(logging.INFO)

    parser = argparse.ArgumentParser(description=f"MiXeR v{__version__}: visualization tools")
    parser.add_argument('--version', action='version', version=f'MiXeR v{__version__}')
    
    parent_parser_figures = argparse.ArgumentParser(add_help=False)
    parent_parser_figures.add_argument('--argsfile', type=open, action=common.utils_cli.LoadFromFile, default=None, help="file with additional command-line arguments")
    parent_parser_figures.add_argument("--out", type=str, default="mixer", help="prefix for the output files")
    parent_parser_figures.add_argument('--ext', type=str, default=['png'], nargs='+', choices=['png', 'svg'], help="output extentions")
    parent_parser_figures.add_argument('--zmax', type=float, default=10, help="limit for z-vs-z density plots")
    parent_parser_figures.add_argument('--statistic', type=str, nargs='+', default=["point_estimate"], choices=["point_estimate", "mean", "median", "std", "min", "max"], help="Which statistic to show in the tables and on the Venn diagrams. Can have multiple values. In the case of venn diagram, the first value (typically 'point_estimate' or 'mean') indicate the size of the venn diagram; the second value (optional, typically 'std') allow to include error bars on the Venn diagramm.")
    parent_parser_figures.add_argument('--dev', default=False, action="store_true", help=argparse.SUPPRESS) 

    subparsers = parser.add_subparsers()
    bivar_mixer.figures.parser_one_add_arguments(func=bivar_mixer.figures.execute_one_parser, parser=subparsers.add_parser("one", parents=[parent_parser_figures], help='produce figures for univariate analysis'))
    bivar_mixer.figures.parser_two_add_arguments(func=bivar_mixer.figures.execute_two_parser, parser=subparsers.add_parser("two", parents=[parent_parser_figures], help='produce figures for cross-trait analysis'))
    bivar_mixer.figures.parser_combine_add_arguments(func=bivar_mixer.figures.execute_combine_parser, parser=subparsers.add_parser("combine", parents=[parent_parser_figures], help='combine .json files MiXeR runs (e.g. with different --extract setting)'))

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args(sys.argv[1:])

    header = MASTHEAD
    header += f"Call: mixer_figures.py {sys.argv[1]} "
    for a in sys.argv[2:]:
        header += f'\\\n\t{a} ' if a.startswith('--') else f'{a} '
    logging.info(header)

    args.func(args)
