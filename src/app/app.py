#!/usr/bin/env python3

import ujson
import numpy
import math
from joblib import Parallel, delayed

from .aggregate import bedaggregate, bedAggregateByName, summit, aggregate
from .sequence.twobit import TwoBitReader
from .batch.batch import batch_size, flatten, BatchedFile

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

def tregion(line):
    x = line.strip().split()
    x[1] = int(x[1]); x[2] = int(x[2])
    if len(x) <= 3: return tuple(x + [ '.' ])
    return tuple(x[:4])

def runzscore(args):
    batchsize = batch_size(args)
    batches = []
    def mean(cbatch):
        cbatch = [ x for x in cbatch if x is not None ]
        _, matrix = aggregate(
            args.signal_file, cbatch, args.extsize, args.j, args.start_index, args.end_index, noextension = args.extsize is None
        )
        return [ sum(x) if len(x) > 0 else 0 for x in matrix ]
    with BatchedFile(args.bed_file, batchsize) as f:
        batches = [ [ summit(x) if args.extsize is not None else tregion(x) for x in batch ] for batch in f ]
    results = flatten( Parallel(n_jobs = args.j)(delayed(mean)(x) for x in batches) )
    mean = numpy.mean([ math.log(x) for x in results if x > 0 ])
    std = numpy.std([ math.log(x) for x in results if x > 0 ])
    minv = math.floor((min([ math.log(x) for x in results if x > 0 ]) - mean) / std)
    results = [ (math.log(x) - mean) / std if x > 0 else minv for x in results ]
    with open(args.bed_file, 'r') as f:
        with open(args.output_file, 'w') as o:
            i = 0
            for line in f:
                o.write("%s\t%.3f\n" % ('\t'.join(line.strip().split()[:4]), results[i]))
                i += 1

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
    first = True
    def write(cbatch, o, first = True):
        _, matrix = aggregate(
            args.signal_file, cbatch, args.extsize, args.j, args.start_index, args.end_index, args.resolution, args.decimal_resolution
        )
        if args.coordinate_map: matrix = { sregion(cbatch[i], args.extsize): x for i, x in enumerate(matrix) }
        o.write(("," if not first else "") + ujson.dumps(matrix)[1:-1])
        first = False
    with open(args.output_file, 'w') as o:
        o.write("{" if args.coordinate_map else "[")
        with BatchedFile(args.bed_file, batchsize) as f:
            for batch in f:
                write([ summit(x) for x in batch ], o, first)
                first = False
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
    return line[0], m - extsize, m + extsize

def runsequence(args):
    if args.streaming:
        runsequence_stream(args)
    else:
        runsequence_all(args)

def runsequence_stream(args):
    batch = []
    batchsize = batch_size(args)
    first = True
    def write(cbatch, t, o, last = False):
        if not args.coordinate_map:
            results = [ t.read(*seqregion(x, args.extsize)) for x in cbatch ]
        else:
            results = { "%s:%s-%s" % seqregion(x, args.extsize): t.read(*seqregion(x, args.extsize)) for x in cbatch }
        o.write(("," if not first else "") + ujson.dumps(results)[1:-1])
    with open(args.output_file, 'w') as o:
        o.write("{" if args.coordinate_map else "[")
        with TwoBitReader(args.two_bit_file) as t:
            with BatchedFile(args.bed_file, batchsize) as f:
                for batch in f:
                    write(batch, t, o, first)
                    first = False
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
