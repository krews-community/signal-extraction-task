#!/usr/bin/env python3

import os
import tempfile
import unittest
import hashlib

from app.app import runaggregate, runmatrix, runsequence, runzscore

class TestInput:
    
    def __init__(
        self, testbed = "test.bed", startindex = 0, endindex = None, resolution = 1, decimal_resolution = 2, grouped = False,
        coordinate_map = False, extsize = 5, streaming = False
    ):
        self.signal_file = os.path.join(os.path.dirname(__file__), "resources", "test.bigWig")
        self.two_bit_file = os.path.join(os.path.dirname(__file__), "resources", "chrTest.2bit")
        self.bed_file = os.path.join(os.path.dirname(__file__), "resources", testbed)
        self.extsize = extsize
        self.j = 1
        self.start_index = startindex
        self.end_index = endindex
        self.resolution = resolution
        self.decimal_resolution = decimal_resolution
        self.grouped = grouped
        self.coordinate_map = coordinate_map
        self.batch_size = 3
        self.streaming = streaming

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
    
    def test_runzscore(self):
        with TestInput(testbed = "test.chrTest.signal.bed") as test:
            runzscore(test)
            self.assertEqual(hashlib.md5(test.output.read()).hexdigest(), "a9309e8fa95fb36d6651842f9bc086ea")
    
    def test_runzscore_no_extension(self):
        with TestInput(testbed = "test.chrTest.signal.bed", extsize = None) as test:
            runzscore(test)
            self.assertEqual(hashlib.md5(test.output.read()).hexdigest(), "0e3397fce1282ab6b4ad26955c2264b9")

    def test_runaggregate_10(self):
        with TestInput(resolution = 2) as test:
            runaggregate(test)
            self.assertEqual(hashlib.md5(test.output.read()).hexdigest(), "069e6d1dd9aed62c90319af1ae3fe65d")

    def test_runmatrix(self):
        with TestInput() as test:
            runmatrix(test)
            self.assertEqual(hashlib.md5(test.output.read()).hexdigest(), "b9b14755e8fd0c29342fb417aa0db286")
        
    def test_runmatrix_coordinate_map(self):
        with TestInput(coordinate_map = True) as test:
            runmatrix(test)
            self.assertEqual(hashlib.md5(test.output.read()).hexdigest(), "53a11a17db1d51d5646542783fa551ac")

    def test_runmatrix_streamed(self):
        with TestInput(streaming = True) as test:
            runmatrix(test)
            streamed = hashlib.md5(test.output.read()).hexdigest()
        with TestInput() as test:
            runmatrix(test)
            self.assertEqual(hashlib.md5(test.output.read()).hexdigest(), streamed)
        
    def test_runmatrix_coordinate_map_streamed(self):
        with TestInput(coordinate_map = True, streaming = True) as test:
            runmatrix(test)
            streamed = hashlib.md5(test.output.read()).hexdigest()
        with TestInput(coordinate_map = True) as test:
            runmatrix(test)
            self.assertEqual(hashlib.md5(test.output.read()).hexdigest(), streamed)

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
        with TestInput(testbed = "test.chrTest.bed") as test:
            runsequence(test)
            self.assertEqual(hashlib.md5(test.output.read()).hexdigest(), "6efa138a53b74e9cdf7489367cde7832")

    def test_runsequence_coordinate_map(self):
        with TestInput(testbed = "test.chrTest.bed", coordinate_map = True, extsize = 7) as test:
            runsequence(test)
            self.assertEqual(hashlib.md5(test.output.read()).hexdigest(), "b6b4aff95070cf5a11a0e4d1617b4eb6")

    def test_runsequence_stream(self):
        with TestInput(testbed = "test.chrTest.bed", extsize = 7, streaming = True) as test:
            runsequence(test)
            streamed = hashlib.md5(test.output.read()).hexdigest()
        with TestInput(testbed = "test.chrTest.bed", extsize = 7) as test:
            runsequence(test)
            self.assertEqual(hashlib.md5(test.output.read()).hexdigest(), streamed)

    def test_runsequence_coordinate_map_stream(self):
        with TestInput(testbed = "test.chrTest.bed", coordinate_map = True, extsize = 7, streaming = True) as test:
            runsequence(test)
            streamed = hashlib.md5(test.output.read()).hexdigest()
        with TestInput(testbed = "test.chrTest.bed", coordinate_map = True, extsize = 7) as test:
            runsequence(test)
            self.assertEqual(hashlib.md5(test.output.read()).hexdigest(), streamed)
