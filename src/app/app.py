#!/usr/bin/env python3

import ujson

from .aggregate import bedaggregate, bedAggregateByName, summit, aggregate
from .sequence.twobit import TwoBitReader

def batch_size(args):
    try:
        return args.batch_size
    except:
        return 1000

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

def sregion(summit, extsize):
    c, s, _ = summit
    return "%s:%d-%d" % (c, s - extsize, s + extsize)

def runmatrix(args):
    if args.streaming:
        runmatrix_stream(args)
    else:
        runmatrix_all(args)

def runmatrix_stream(args):
    batch = []
    batchsize = batch_size(args)
    def write(cbatch, o, last = False):
        _, matrix = aggregate(
            args.signal_file, cbatch, args.extsize, args.j, args.start_index, args.end_index, args.resolution, args.decimal_resolution
        )
        if args.coordinate_map: matrix = { sregion(cbatch[i], args.extsize): x for i, x in enumerate(matrix) }
        o.write(ujson.dumps(matrix)[1:-1] + ("," if not last else ""))
    with open(args.output_file, 'w') as o:
        o.write("{" if args.coordinate_map else "[")
        with open(args.bed_file, 'r') as f:
            for line in f:
                if len(batch) == batchsize:
                    write(batch, o)
                    batch = []
                batch.append(summit(line))
        write(batch, o, last = True)
        o.write("}\n" if args.coordinate_map else "]\n")

def runmatrix_all(args):
    _, matrix = bedaggregate(
        args.signal_file, args.bed_file, args.extsize, args.j, args.start_index, args.end_index, args.resolution, args.decimal_resolution
    )
    if args.coordinate_map:
        with open(args.bed_file, 'r') as f:
            matrix = { sregion(summit(f.readline()), args.extsize): x for x in matrix }
    with open(args.output_file, 'w') as o:
        o.write(ujson.dumps(matrix) + '\n')

def seqregion(line, extsize):
    line = line.strip().split()
    m = int((int(line[1]) + int(line[2])) / 2)
    return line[0], int(line[1]) - extsize, int(line[2]) + extsize

def runsequence(args):
    if args.streaming:
        runsequence_stream(args)
    else:
        runsequence_all(args)

def runsequence_stream(args):
    batch = []
    batchsize = batch_size(args)
    def write(cbatch, t, o, last = False):
        if not args.coordinate_map:
            results = [ t.read(*seqregion(x, args.extsize)) for x in cbatch ]
        else:
            results = { "%s:%s-%s" % seqregion(x, args.extsize): t.read(*seqregion(x, args.extsize)) for x in cbatch }
        o.write(ujson.dumps(results)[1:-1] + ("," if not last else ""))
    with open(args.output_file, 'w') as o:
        o.write("{" if args.coordinate_map else "[")
        with TwoBitReader(args.two_bit_file) as t:
            with open(args.bed_file, 'r') as f:
                for line in f:
                    if len(batch) == batchsize:
                        write(batch, t, o)
                        batch = []
                    batch.append(line)
            write(batch, t, o, last = True)
        o.write("}\n" if args.coordinate_map else "]\n")

def runsequence_all(args):
    with TwoBitReader(args.two_bit_file) as t:
        with open(args.bed_file, 'r') as f:
            if not args.coordinate_map:
                results = [ t.read(*seqregion(x, args.extsize)) for x in f ]
            else:
                results = { "%s:%s-%s" % seqregion(x, args.extsize): t.read(*seqregion(x, args.extsize)) for x in f }
    with open(args.output_file, 'w') as o:
        o.write(ujson.dumps(results) + '\n')
