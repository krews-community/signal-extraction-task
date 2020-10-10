#!/usr/bin/env python3

from typing import List

ONEHOT = {
    'a': [ 1, 0, 0, 0 ],
    'c': [ 0, 1, 0, 0 ],
    'g': [ 0, 0, 1, 0 ],
    't': [ 0, 0, 0, 1 ],
    'n': [ 0, 0, 0, 0 ]
}

def onehot(sequence: str) -> List[int]:
    return [
        ONEHOT[x.lower()] for x in sequence.replace(" ", "")
    ]
