#!/usr/bin/env python3

import ujson

from .aggregate import bedaggregate, bedAggregateByName
from .sequence.twobit import TwoBitReader

def runaggregate(args):
    if not args.grouped:
        values, _ = bedaggregate(
            args.signal_file, args.bed_file, args.extsize, args.j, args.start_index, args.end_index, args.resolution, args.decimal_resolution
        )
    else:
        values = bedAggregateByName(
            args.signal_file, args.bed_file, args.extsize, args.j, args.start_index, args.end_index, args.resolution, args.decimal_resolution
        )
    with open(args.output_file, 'w') as o:
        o.write(ujson.dumps(values) + '\n')

def runmatrix(args):
    def region(line):
        line = line.strip().split()
        return "%s:%s-%s" % (line[0], int(line[1]), int(line[2]))
    _, matrix = bedaggregate(
        args.signal_file, args.bed_file, args.extsize, args.j, args.start_index, args.end_index, args.resolution, args.decimal_resolution
    )
    if args.coordinate_map:
        with open(args.bed_file, 'r') as f:
            matrix = { region(f.readline()): x for x in matrix }
    with open(args.output_file, 'w') as o:
        o.write(ujson.dumps(matrix) + '\n')

def runsequence(args):
    def region(line):
        line = line.strip().split()
        m = int((int(line[1]) + int(line[2])) / 2)
        return line[0], int(line[1]) - args.extsize, int(line[2]) + args.extsize
    with TwoBitReader(args.two_bit_file) as t:
        with open(args.bed_file, 'r') as f:
            if not args.coordinate_map:
                results = [ t.read(*region(x)) for x in f ]
            else:
                results = { "%s:%s-%s" % region(x): t.read(*region(x)) for x in f }
    with open(args.output_file, 'w') as o:
        o.write(ujson.dumps(results) + '\n')
