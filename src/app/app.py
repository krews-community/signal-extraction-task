#!/usr/bin/env python3

import ujson

from .aggregate import bedaggregate, bedAggregateByName
from .sequence.twobit import TwoBitReader

def runaggregate(args):
    if not args.grouped:
        values, _ = bedaggregate(args.signal_file, args.bed_file, args.extsize, args.j, args.start_index, args.end_index, args.resolution)
    else:
        values = bedAggregateByName(args.signal_file, args.bed_file, args.extsize, args.j, args.start_index, args.end_index, args.resolution)
    with open(args.output_file, 'w') as o:
        o.write(ujson.dumps(values) + '\n')

def runmatrix(args):
    _, matrix = bedaggregate(args.signal_file, args.bed_file, args.extsize, args.j, args.start_index, args.end_index, args.resolution)
    with open(args.output_file, 'w') as o:
        o.write(ujson.dumps(matrix) + '\n')

def runsequence(args):
    def region(line):
        line = line.strip().split()
        return line[0], int(line[1]), int(line[2])
    with TwoBitReader(args.two_bit_file) as t:
        with open(args.bed_file, 'r') as f:
            results = [ t.read(*region(x)) for x in f ]
    with open(args.output_file, 'w') as o:
        o.write(ujson.dumps(results) + '\n')
