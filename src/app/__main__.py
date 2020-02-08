#!/usr/bin/env python3

import sys
import argparse

from .app import runaggregate, runmatrix

def main():
    
    parser = argparse.ArgumentParser(description = "Aggregates signal for a number of regions from a BED file.")
    subparsers = parser.add_subparsers()
    aggregate = subparsers.add_parser("aggregate", help = "aggregate signal across regions")
    aggregate.add_argument("--bed-file", type = str, help = "Path to the BED file with the regions to aggregate.", required = True)
    aggregate.add_argument("--signal-file", type = str, help = "Path to a BigWig file with the signal to aggregate.", required = True)
    aggregate.add_argument("--output-file", type = str, help = "Path to write the output, in JSON format.", required = True)
    aggregate.add_argument("--extsize", type = int, help = "number of basepairs to expand each region; default 500.", default = 500)
    aggregate.add_argument("-j", type = int, help = "number of cores to use in parallel; default 8.", default = 8)
    aggregate.set_defaults(func = runaggregate)
    matrix = subparsers.add_parser("matrix", help = "produce a signal matrix for the given regions")
    matrix.add_argument("--bed-file", type = str, help = "Path to the BED file with the regions to aggregate.", required = True)
    matrix.add_argument("--signal-file", type = str, help = "Path to a BigWig file with the signal to aggregate.", required = True)
    matrix.add_argument("--output-file", type = str, help = "Path to write the output, in JSON format.", required = True)
    matrix.add_argument("--extsize", type = int, help = "number of basepairs to expand each region; default 500.", default = 500)
    matrix.add_argument("-j", type = int, help = "number of cores to use in parallel; default 8.", default = 8)
    matrix.set_defaults(func = runmatrix)

    args = parser.parse_args()
    args.func(args)

    return 0

if __name__ == "__main__":
    sys.exit(main())
