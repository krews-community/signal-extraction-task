#!/usr/bin/env python3

import sys
import argparse

from .app import runaggregate, runmatrix, runsequence, runzscore
from .aggregate import bedaggregate
from .sequence.twobit import TwoBitReader

def main():
    
    parser = argparse.ArgumentParser(description = "Aggregates signal for a number of regions from a BED file.")
    subparsers = parser.add_subparsers()
    
    aggregate = subparsers.add_parser("aggregate", help = "aggregate signal across regions")
    aggregate.add_argument("--bed-file", type = str, help = "Path to the BED file with the regions to aggregate.", required = True)
    aggregate.add_argument("--signal-file", type = str, help = "Path to a BigWig file with the signal to aggregate.", required = True)
    aggregate.add_argument("--output-file", type = str, help = "Path to write the output, in JSON format.", required = True)
    aggregate.add_argument("--extsize", type = int, help = "number of basepairs to expand each region; default 500.", default = 500)
    aggregate.add_argument("--start-index", type = int, help = "Index of the first element to aggregate (inclusive); default 1.", default = 0)
    aggregate.add_argument("--end-index", type = int, help = "Index of the last element to aggregate (not inclusive); defaults to the end of the list.", default = None)
    aggregate.add_argument("--resolution", type = int, help = "Bin size to use in basepairs; defaults to 1 (a single basepair)", default = 1)
    aggregate.add_argument("-j", type = int, help = "number of cores to use in parallel; default 8.", default = 8)
    aggregate.add_argument("--decimal-resolution", type = int, help = "Number of decimal places to keep in output.", default = 2)
    aggregate.add_argument("--grouped", action = "store_true", help = "If specified, groups output by the name field of each BED line.", default = False)
    aggregate.set_defaults(func = runaggregate)
    
    matrix = subparsers.add_parser("matrix", help = "produce a signal matrix for the given regions")
    matrix.add_argument("--bed-file", type = str, help = "Path to the BED file with the regions to aggregate.", required = True)
    matrix.add_argument("--signal-file", type = str, help = "Path to a BigWig file with the signal to aggregate.", required = True)
    matrix.add_argument("--output-file", type = str, help = "Path to write the output, in JSON format.", required = True)
    matrix.add_argument("--extsize", type = int, help = "number of basepairs to expand each region; default 500.", default = 500)
    matrix.add_argument("--start-index", type = int, help = "Index of the first element to aggregate (inclusive); default 1.", default = 0)
    matrix.add_argument("--end-index", type = int, help = "Index of the last element to aggregate (not inclusive); defaults to the end of the list.", default = None)
    matrix.add_argument("--resolution", type = int, help = "Bin size to use in basepairs; defaults to 1 (a single basepair)", default = 1)
    matrix.add_argument("-j", type = int, help = "number of cores to use in parallel; default 8.", default = 8)
    matrix.add_argument("--decimal-resolution", type = int, help = "Number of decimal places to keep in output.", default = 2)
    matrix.add_argument("--coordinate-map", action = "store_true", default = False, help = "if set, output JSON maps coordinates to values")
    matrix.add_argument("--streaming", action = "store_true", default = False, help = "if set, batches of results are streamed to an output file rather than kept in memory")
    matrix.add_argument("--random-access", action = "store_true", default = False, help = "if set, writes output in a format designed for seeking and requesting by index")
    matrix.set_defaults(func = runmatrix)

    sequence = subparsers.add_parser("sequence", help = "extract one hot encoded sequence for the given regions from a 2bit file")
    sequence.add_argument("--bed-file", type = str, help = "Path to the BED file with the regions for which to extract sequence.", required = True)
    sequence.add_argument("--two-bit-file", type = str, help = "Path to a 2bit file with the genomic sequence.", required = True)
    sequence.add_argument("--output-file", type = str, help = "Path to write the output, in JSON format.", required = True)
    sequence.add_argument("--extsize", type = int, help = "number of basepairs to expand each region; default 500.", default = 500)
    sequence.add_argument("-j", type = int, help = "number of cores to use in parallel; default 8.", default = 8)
    sequence.add_argument("--coordinate-map", action = "store_true", default = False, help = "if set, output JSON maps coordinates to values")
    sequence.add_argument("--streaming", action = "store_true", default = False, help = "if set, batches of results are streamed to an output file rather than kept in memory")
    sequence.add_argument("--random-access", action = "store_true", default = False, help = "if set, writes output in a format designed for seeking and requesting by index")
    sequence.set_defaults(func = runsequence)

    zscore = subparsers.add_parser("zscore", help = "computes Z-scores for aggregated signal across a list of regions")
    zscore.add_argument("--bed-file", type = str, help = "Path to the BED file with the regions to aggregate.", required = True)
    zscore.add_argument("--signal-file", type = str, help = "Path to a BigWig file with the signal to aggregate.", required = True)
    zscore.add_argument("--output-file", type = str, help = "Path to write the output, in BED format.", required = True)
    zscore.add_argument("--start-index", type = int, help = "Index of the first element to aggregate (inclusive); default 1.", default = 0)
    zscore.add_argument("--end-index", type = int, help = "Index of the last element to aggregate (not inclusive); defaults to the end of the list.", default = None)
    zscore.add_argument("--extsize", type = int, help = "if passed, extends each region by the given number of basepairs in each direction around the center points.", default = None)
    zscore.add_argument("--json", action = "store_true", default = False, help = "if set, writes output as a JSON array")
    zscore.add_argument("-j", type = int, help = "number of cores to use in parallel; default 8.", default = 8)
    zscore.set_defaults(func = runzscore)

    args = parser.parse_args()
    args.func(args)

    return 0

if __name__ == "__main__":
    sys.exit(main())
