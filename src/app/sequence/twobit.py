#!/usr/bin/env python3

import twobitreader
from typing import List

from .onehot import onehot

class TwoBitReader:

    def __init__(self, path: str):
        self.path = path
    
    def __enter__(self):
        self.twobit = twobitreader.TwoBitFile(self.path)
        return self
    
    def __exit__(self, *args):
        self.twobit.close()

    def read(self, chromosome: str, start: int, end: int) -> List[int]:
        return onehot(self.twobit[chromosome][start:end])
