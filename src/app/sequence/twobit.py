#!/usr/bin/env python3

import twobitreader
from typing import List

from .onehot import onehot, ONEHOT

class TwoBitReader:

    @staticmethod
    def pad(a, l):
        return a + [ ONEHOT['n'] for _ in range(l - len(a)) ]

    def __init__(self, path: str):
        self.path = path
    
    def __enter__(self):
        self.twobit = twobitreader.TwoBitFile(self.path)
        return self
    
    def __exit__(self, *args):
        self.twobit.close()

    def read(self, chromosome: str, start: int, end: int) -> List[int]:
        try:
            if start < 0:
                return [ ONEHOT['n'] for _ in range(start, 0) ] + onehot(self.twobit[chromosome][:end])
            return TwoBitReader.pad( onehot(self.twobit[chromosome][start:end]), end - start )
        except:
            return [ ONEHOT['n'] for _ in range(end - start) ]
