#!/usr/bin/env python3

import os
import tempfile
import unittest
import hashlib

from app.app import runaggregate, runmatrix, runsequence

class TestInput:
    
    def __init__(self, testbed = "test.bed", startindex = 0, endindex = None, resolution = 1, decimal_resolution = 2, grouped = False):
        self.signal_file = os.path.join(os.path.dirname(__file__), "resources", "test.bigWig")
        self.two_bit_file = os.path.join(os.path.dirname(__file__), "resources", "chrTest.2bit")
        self.bed_file = os.path.join(os.path.dirname(__file__), "resources", testbed)
        self.extsize = 5
        self.j = 1
        self.start_index = startindex
        self.end_index = endindex
        self.resolution = resolution
        self.decimal_resolution = decimal_resolution
        self.grouped = grouped

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
            self.assertEqual(hashlib.md5(test.output.read()).hexdigest(), "fd76fa57c7a036b6798a3f6b4ad3a944")

    def test_runaggregate_endindex(self):
        with TestInput(testbed = "test.offsets.bed", endindex = 2) as test:
            runaggregate(test)
            self.assertEqual(hashlib.md5(test.output.read()).hexdigest(), "5e0f8d2c34d3c796d46ad70b204f9068")

    def test_runaggregate_grouped(self):
        with TestInput(testbed = "test.group.bed", grouped = True) as test:
            runaggregate(test)
            self.assertEqual(hashlib.md5(test.output.read()).hexdigest(), "05355371a11e2df119eab0ebadd99dd0")

    def test_runsequence(self):
        with TestInput(testbed = "test.chrTest.bed", grouped = True) as test:
            runsequence(test)
            self.assertEqual(hashlib.md5(test.output.read()).hexdigest(), "a4271d6a803b6801bbbbe48e62d6a826")
