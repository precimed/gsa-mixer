#!/usr/bin/env python
'''
(c) 2016-2020 Oleksandr Frei, Alexey A. Shadrin, Dominic Holland
MiXeR software: Univariate and Bivariate Causal Mixture for GWAS
'''

import logging
import sys

import bivar_mixer.figures
from version import VERSION as __version__, MASTHEAD

if __name__ == "__main__":
    logging.getLogger().setLevel(logging.INFO)

    parser = bivar_mixer.figures.generate_args_parser(__version__)

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args(sys.argv[1:])

    header = MASTHEAD
    header += f"Call: mixer.py {sys.argv[1]} "
    for a in sys.argv[2:]:
        header += f'\\\n\t{a} ' if a.startswith('--') else f'{a} '
    logging.info(header)

    args.func(args)
