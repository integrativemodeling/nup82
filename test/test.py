#!/usr/bin/env python

import unittest
import os
import sys
import subprocess

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), '..'))

class Tests(unittest.TestCase):

    def test_simple(self):
        """Test model building"""
        os.chdir(os.path.join(TOPDIR, 'template'))
        p = subprocess.check_call(["python", "1_modeling_initial_random.py",
                     "-r", "1000", "-out", "output",
                     "-em2d", "../data/em2d/2.pgm", "-weight", "10000.0"])
        # todo: assert outputs

if __name__ == '__main__':
    unittest.main()
