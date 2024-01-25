#!/usr/bin/env python

import logging
import sys

from mixer.figures import parse_args

if __name__ == "__main__":
    logging.getLogger().setLevel(logging.INFO)
    args = parse_args(sys.argv[1:])
    args.func(args)
