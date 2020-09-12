#!/usr/bin/env python3

import ujson

from .aggregate import bedaggregate, bedAggregateByName

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
