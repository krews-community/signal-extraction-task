#!/usr/bin/env python3

import os
import tempfile
import unittest
import hashlib

from app.app import runaggregate, runmatrix

class TestInput:
    
    def __init__(self, testbed = "test.bed", startindex = 0, endindex = None, resolution = 1):
        self.signal_file = os.path.join(os.path.dirname(__file__), "resources", "test.bigWig")
        self.bed_file = os.path.join(os.path.dirname(__file__), "resources", testbed)
        self.extsize = 5
        self.j = 1
        self.start_index = startindex
        self.end_index = endindex
        self.resolution = resolution

    def __enter__(self):
        self.output = tempfile.NamedTemporaryFile()
        self.output_file = self.output.name
        return self

    def __exit__(self, exc_type, exc_value, tb):
        if exc_type is not None:
            raise
        self.output.close()

class TestApp(unittest.TestCase):
        
    def test_runaggregate(self):
        with TestInput() as test:
            runaggregate(test)
            self.assertEqual(hashlib.md5(test.output.read()).hexdigest(), "ef3de94cf04e83f8b924aa1275692ef7")

    def test_runaggregate_10(self):
        with TestInput(resolution = 2) as test:
            runaggregate(test)
            self.assertEqual(hashlib.md5(test.output.read()).hexdigest(), "069e6d1dd9aed62c90319af1ae3fe65d")

    def test_runmatrix(self):
        with TestInput() as test:
            runmatrix(test)
            self.assertEqual(hashlib.md5(test.output.read()).hexdigest(), "b9b14755e8fd0c29342fb417aa0db286")

    def test_runmatrix_10(self):
        with TestInput(resolution = 2) as test:
            runmatrix(test)
            self.assertEqual(hashlib.md5(test.output.read()).hexdigest(), "04224186a7d33f8b5379a287d0a48cec")

    def test_runaggregate_startindex(self):
        with TestInput(testbed = "test.offsets.bed", startindex = 1) as test:
            runaggregate(test)
            self.assertEqual(hashlib.md5(test.output.read()).hexdigest(), "3679fd43377e574c0281e51247e0f019")

    def test_runaggregate_endindex(self):
        with TestInput(testbed = "test.offsets.bed", endindex = 2) as test:
            runaggregate(test)
            self.assertEqual(hashlib.md5(test.output.read()).hexdigest(), "5e0f8d2c34d3c796d46ad70b204f9068")
