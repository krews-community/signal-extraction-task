#!/usr/bin/env python

import os
import numpy
import pyBigWig
import math
import gzip

from joblib import Parallel, delayed

def read(regions, divpoints, bigwig, i):
    """
    Reads signal for a given set of regions from a BigWig. Intended to be called in Parallel; data are read only
    for the subset of regions specified by the given index. For each region, if a strand is present and the region
    is on the minus strand, the order of the signal values is reversed for aggregation purposes.

    Args:
        regions (list): complete list of regions to read, each as a tuple of chromosome, start, end, strand
        divpoints (list): list of offsets in the region list at which each individual job should start
        bigwig (string): path to the BigWig file to read
        i (int): index in the divpoints list which this job should handle

    Returns:
        A list of value vectors, one per region beginning at index divpoints[i] in the regions list and ending
        at index divpoints[i + 1] - 1, inclusive.
    """
    bw = pyBigWig.open(bigwig)
    if bw is None:
        raise Exception("Error opening %s: no such file or directory." % bigwig)
    values = [
        bw.values(*tuple(region[:3])) if region[3] != '-' else list(reversed(bw.values(*tuple(region[:3]))))
        for region in regions[divpoints[i]:divpoints[i + 1]]
    ]
    return values

def valuematrix(bigwig, centers, extsize, j = 8):
    """
    Reads signal for a given set of regions from a BigWig in parallel. Each region is uniformly sized around the
    centers points passed in the "centers" parameter. For each region, if a strand is present and the region is on
    the minus strand, the order of the signal values is reversed for aggregation purposes.

    Args:
        bigwig (string): path to the BigWig to read
        centers (list): list of regions, each as a tuple of chromosome, center position, strand
        extsize (int): number of basepairs to read around the center point
        j (int): number of threads to use; default is 8

    Returns:
        A matrix of signal values for each region; the rows are the regions in their original order, and each column
        represents a single basepair.
    """
    bw = pyBigWig.open(bigwig)
    if bw is None:
        raise Exception("Error opening %s: no such file or directory." % bigwig)
    regions = [
        (chrom, x - extsize, x + extsize, strand)
        for chrom, x, strand in centers if x >= extsize and bw.chroms(chrom) is not None and x + extsize < bw.chroms(chrom)
    ]
    divpoints = [ math.floor(len(regions) * i / j) for i in range(j) ] if len(regions) > j else list(range(len(regions)))
    divpoints += [ divpoints[-1] + 1 ]
    readregions = Parallel(n_jobs = j)(delayed(read)(regions, divpoints, bigwig, i) for i in range(len(divpoints) - 1))
    retval = []
    for readregionset in readregions:
        retval += readregionset
    return retval

def aggregate(bigwig, centers, extsize, j = 8):
    """
    Aggregates signal around a given set of center points from the given BigWig file. For each region, if a strand
    is present and the region is on the minus strand, the order of the signal values is reversed.

    Args:
        bigwig (string): path to the BigWig to read.
        centers (list): list of regions, each as a tuple of chromosome, center position, strand
        extsize (int): number of basepairs to read around the center point; regions are extended by this amount in both directions
        j (int): number of threads to use; default is 8

    Returns:
        Tuple of aggregated and matrix-form results. The first element is a single vector of signal values, where each position
        is a single basepair and the value is the average of the signal values from each region at that basepair. The second
        element is the complete signal matrix, where each row is a region in the order of the "centers" parameter and each
        column is a basepair.
    """
    matrix = numpy.nan_to_num(valuematrix(bigwig, centers, extsize, j))
    return [ numpy.mean([ (x[i] if not numpy.isnan(x[i]) else 0.0) for x in matrix ]) for i in range(extsize * 2) ], [ [ float(x) for x in xx ] for xx in matrix ]

def bedaggregate(bigwig, bed, extsize, j = 8):
    """
    Aggregates signal around the center points of each region from a BED file, using signal from the given BigWig file.
    If the BED file has strand information, regions on the minus strand will be inverted before being aggregated; otherwise,
    all regions are assumed to be the same orientation.

    Args:
        bigwig (string): path to the BigWig to read
        bed (string): path to the BED file containing the regions to aggregate
        extsize (int): number of basepairs to read around the center point; regions are extended by this amount in both directions
        j (int): number of threads to use; default is 8

    Returns:
        Tuple of aggregated and matrix-form results. The first element is a single vector of signal values, where each position
        is a single basepair and the value is the average of the signal values from each region at that basepair. The second
        element is the complete signal matrix, where each row is a region in the order of the "centers" parameter and each
        column is a basepair.
    """
    if not os.path.exists(bed):
        raise Exception("Error opening %s: no such file or directory." % bed)
    with (gzip.open if bed.endswith(".gz") else open)(bed, 'rt') as f:
        centers = [(
            l.split('\t')[0],
            int((int(l.split('\t')[2].strip()) + int(l.split('\t')[1])) / 2),
            l.split('\t')[3].strip() if len(l.split('\t')) >= 4 else '.'
        ) for l in f ]
    return aggregate(bigwig, centers, extsize, j)
