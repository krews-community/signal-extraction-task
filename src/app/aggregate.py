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
        regions (list): complete list of regions to read, each as a tuple of chromosome, start, end, strand or None
        divpoints (list): list of offsets in the region list at which each individual job should start
        bigwig (string): path to the BigWig file to read
        i (int): index in the divpoints list which this job should handle

    Returns:
        A list of value vectors, one per region beginning at index divpoints[i] in the regions list and ending
        at index divpoints[i + 1] - 1, inclusive. Regions which are None are assumed to be missing and return an
        empty vector.
    """
    bw = pyBigWig.open(bigwig)
    if bw is None:
        raise Exception("Error opening %s: no such file or directory." % bigwig)
    values = [
        ( bw.values(*tuple(region[:3])) if region[3] != '-' else list(reversed(bw.values(*tuple(region[:3])))) ) if region is not None else []
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
        represents a single basepair. Regions which are out of range for the given bigWig are represented by empty vectors.
    """
    bw = pyBigWig.open(bigwig)
    if bw is None:
        raise Exception("Error opening %s: no such file or directory." % bigwig)
    regions = [
        (chrom, x - extsize, x + extsize, strand) if x >= extsize and bw.chroms(chrom) is not None and x + extsize < bw.chroms(chrom) else None
        for chrom, x, strand in centers
    ]
    divpoints = [ math.floor(len(regions) * i / j) for i in range(j) ] if len(regions) > j else list(range(len(regions)))
    divpoints.append(len(regions))
    if j == 1: divpoints = [ 0, len(regions) ]
    readregions = Parallel(n_jobs = j)(delayed(read)(regions, divpoints, bigwig, i) for i in range(len(divpoints) - 1))
    retval = []
    for readregionset in readregions:
        retval += readregionset
    return retval

def condense(a, r = 1, dr = 2):
    if r == 1: return a
    rv = [ 0. for x in range(int(len(a) / r)) ]
    for i in range(len(a)):
        rv[int(math.floor(i / r))] += a[i]
    return [ round(x, dr) for x in rv ]

def aggregate(bigwig, centers, extsize, j = 8, startindex = 0, endindex = None, resolution = 1, decimal_resolution = 2):
    """
    Aggregates signal around a given set of center points from the given BigWig file. For each region, if a strand
    is present and the region is on the minus strand, the order of the signal values is reversed.

    Args:
        bigwig (string): path to the BigWig to read.
        centers (list): list of regions, each as a tuple of chromosome, center position, strand
        extsize (int): number of basepairs to read around the center point; regions are extended by this amount in both directions
        j (int): number of threads to use; default is 8
        startindex (int): first index to include in the aggregate; default is 0
        endindex (int): last index to aggregate (not inclusive); default is None, indicating aggregation should continue to the end of the list.
        resolution (int): if set, returns bins which represent the average signal across this number of basepairs
        decimal_resolution (int): rounds values in the matrix and aggregate vector to the given number of decimal places

    Returns:
        Tuple of aggregated and matrix-form results. The first element is a single vector of signal values, where each position
        is a bin with a size determined by the resolution parameter and the value is the average of the signal values from each
        region at that bin. The second element is the complete signal matrix, where each row is a region in the order of the "centers"
        parameter and each column is a bin. Regions which are missing from or out of range for the bigWig are represented by vectors
        of zeroes.
    """
    if endindex is None: endindex = len(centers)
    if endindex <= startindex: return [], []
    matrix = valuematrix(bigwig, centers[startindex : endindex], extsize, j)
    for i in range(len(matrix)):
        if len(matrix[i]) == 0: matrix[i] = [ 0. for _ in range(extsize * 2) ]
    matrix = [ condense([ float(x) for x in xx ], resolution) for xx in numpy.nan_to_num(matrix) ]
    aggregate = [
        round(numpy.mean([ (x[i] if not numpy.isnan(x[i]) else 0.) for x in matrix ]), decimal_resolution)
        for i in range(len(matrix[0]))
    ]
    return aggregate, matrix

def bedaggregate(bigwig, bed, extsize, j = 8, startindex = 0, endindex = None, resolution = 1, decimal_resolution = 2):
    """
    Aggregates signal around the center points of each region from a BED file, using signal from the given BigWig file.
    If the BED file has strand information, regions on the minus strand will be inverted before being aggregated; otherwise,
    all regions are assumed to be the same orientation.

    Args:
        bigwig (string): path to the BigWig to read
        bed (string): path to the BED file containing the regions to aggregate
        extsize (int): number of basepairs to read around the center point; regions are extended by this amount in both directions
        j (int): number of threads to use; default is 8
        startindex (int): first index to aggregate (inclusive); default is 0
        endindex (int): last index to aggregate (not inclusive); default is None, indicating aggregation should continue to the end of the list
        resolution (int): if set, returns bins which represent the average signal across this number of basepairs
        decimal_resolution (int): rounds values in the matrix and aggregate vector to the given number of decimal places

    Returns:
        Tuple of aggregated and matrix-form results. The first element is a single vector of signal values, where each position
        is a bin with a size determined by the resolution parameter and the value is the average of the signal values from each
        region at that basepair. The second element is the complete signal matrix, where each row is a region in the order of the
        "centers" parameter and each column is a bin. Regions which are missing from or out of range for the bigWig are represented
        by vectors of zeroes.
    """
    if not os.path.exists(bed):
        raise Exception("Error opening %s: no such file or directory." % bed)
    with (gzip.open if bed.endswith(".gz") else open)(bed, 'rt') as f:
        centers = [(
            l.split('\t')[0],
            int((int(l.split('\t')[2].strip()) + int(l.split('\t')[1])) / 2),
            l.split('\t')[3].strip() if len(l.split('\t')) >= 4 else '.'
        ) for l in f ]
    return aggregate(bigwig, centers, extsize, j, startindex, endindex, resolution, decimal_resolution)

def bedAggregateByName(bigwig, bed, extsize, j = 8, startindex = 0, endindex = None, resolution = 1, decimal_resolution = 2):
    """
    Aggregates signal around the center points of each region from a BED file, using signal from the given BigWig file.
    If the BED file has strand information, regions on the minus strand will be inverted before being aggregated; otherwise,
    all regions are assumed to be the same orientation. Results are grouped by the names in the fourth field of the BED file;
    strand information, if present, must be given in the fifth field.

    Args:
        bigwig (string): path to the BigWig to read
        bed (string): path to the BED file containing the regions to aggregate
        extsize (int): number of basepairs to read around the center point; regions are extended by this amount in both directions
        j (int): number of threads to use; default is 8
        startindex (int): first index to aggregate (inclusive); default is 0
        endindex (int): last index to aggregate (not inclusive); default is None, indicating aggregation should continue to the end of the list
        resolution (int): if set, returns bins which represent the average signal across this number of basepairs
        decimal_resolution (int): rounds values in the matrix and aggregate vector to the given number of decimal places

    Returns:
        Dictionary of aggregated results. Keys are names from the fourth BED field. Values are vectors of signal values, where each
        position is a bin with a size determined by the resolution parameter and the value is the average of the signal values from each
        region at that basepair.
    """
    if not os.path.exists(bed):
        raise Exception("Error opening %s: no such file or directory." % bed)
    centers = {}
    with (gzip.open if bed.endswith(".gz") else open)(bed, 'rt') as f:
        for line in f:
            line = line.strip().split()
            if line[3] not in centers: centers[line[3]] = []
            centers[line[3]].append((
                line[0],
                int((int(line[2]) + int(line[1])) / 2),
                line[4] if len(line) >= 5 else '.'
            ))
    return { k: aggregate(bigwig, v, extsize, j, startindex, endindex, resolution, decimal_resolution)[0] for k, v in centers.items() }
