#!/usr/bin/env python

import logging
import sys

from version import VERSION as __version__
from mixer.cli import generate_args_parser

if __name__ == "__main__":
    logging.getLogger().setLevel(logging.INFO)
    parser = generate_args_parser(__version__)
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args(sys.argv[1:])
    args.func(args)
