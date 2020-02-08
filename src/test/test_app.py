#!/usr/bin/env python3

import os
import tempfile
import unittest
import hashlib

from app.app import runaggregate, runmatrix

class TestInput:
    
    def __init__(self):
        self.signal_file = os.path.join(os.path.dirname(__file__), "resources", "test.bigWig")
        self.bed_file = os.path.join(os.path.dirname(__file__), "resources", "test.bed")
        self.extsize = 5
        self.j = 1

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

    def test_runmatrix(self):
        with TestInput() as test:
            runmatrix(test)
            self.assertEqual(hashlib.md5(test.output.read()).hexdigest(), "b9b14755e8fd0c29342fb417aa0db286")
