#!/usr/bin/env python

def batch_size(args):
    try:
        return args.batch_size
    except:
        return 1000

def flatten(l):
    r = []
    for x in l:
        r += x
    return r

class BatchedFile:
    
    def __init__(self, file, batchsize):
        self.file = file
        self.batchsize = batchsize
    
    def __enter__(self):
        self.handle = open(self.file, 'r')
        return self
    
    def __iter__(self):
        return self
    
    def __next__(self):
        results = [ x for x in [ self.handle.readline().strip() for _ in range(self.batchsize) ] if len(x.strip().split()) >= 3 ]
        if len(results) == 0: raise StopIteration
        return results

    def __exit__(self, *args):
        self.handle.close()
