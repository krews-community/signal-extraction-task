#!/usr/bin/env python3

import ujson

from .aggregate import bedaggregate

def runaggregate(args):
    values, _ = bedaggregate(args.signal_file, args.bed_file, args.extsize, args.j, args.startindex, args.endindex)
    with open(args.output_file, 'w') as o:
        o.write(ujson.dumps(values) + '\n')

def runmatrix(args):
    _, matrix = bedaggregate(args.signal_file, args.bed_file, args.extsize, args.j, args.startindex, args.endindex)
    with open(args.output_file, 'w') as o:
        o.write(ujson.dumps(matrix) + '\n')
