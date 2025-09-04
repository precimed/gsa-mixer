#!/usr/bin/env python

import logging
import sys

from version import VERSION as __version__
import gsa_mixer.cli
import bivar_mixer.cli

if __name__ == "__main__":
    logging.getLogger().setLevel(logging.INFO)
    if sys.argv[1] in ['ld', 'snps', 'plsa', 'split_sumstats']:
        parser = gsa_mixer.cli.generate_args_parser(__version__)
    else:
        parser = bivar_mixer.cli.generate_args_parser(__version__)
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args(sys.argv[1:])
    args.func(args)
