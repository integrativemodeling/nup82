#!/usr/bin/env python

import unittest
import os
import sys
import subprocess

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), '..'))

class Tests(unittest.TestCase):

    def run_modeller_script(self, script_dir, script_name):
        """Run a Modeller script"""
        os.chdir(os.path.join(TOPDIR, 'data', 'MODELLER',
                              'MODELLER_scripts_5CWS', script_dir))
        # Run script
        p = subprocess.check_call(["python", script_name, "--test"])

    def check_modeller_model(self, model_name, resrng):
        """Test the model output from Modeller"""
        # Make sure PDB was produced with the requested residue range
        with open('%s.B99990001.pdb' % model_name) as fh:
            pdb_lines = [x for x in fh.readlines() if x.startswith('ATOM')]
        rng = (int(pdb_lines[0][22:26]), int(pdb_lines[-1][22:26]))
        self.assertEqual(rng, resrng)

    def test_simple(self):
        """Test model building"""
        os.chdir(os.path.join(TOPDIR, 'template'))
        p = subprocess.check_call(["python", "1_modeling_initial_random.py",
                     "-r", "1000", "-out", "output",
                     "-em2d", "../data/em2d/2.pgm", "-weight", "10000.0"])
        # todo: assert outputs

    def test_nsp1(self):
        """Test generation of comparative models for Nsp1"""
        self.run_modeller_script('Nsp1', 'Build-Nsp1.py')
        self.check_modeller_model('NS1', (637, 727))
        self.check_modeller_model('NS2', (742, 778))
        self.check_modeller_model('NS3', (788, 823))

    def test_nup159(self):
        """Test generation of comparative model for Nup159"""
        self.run_modeller_script('Nup159', 'Build-Nup159.py')
        self.check_modeller_model('Nup159', (1211, 1412))

    def test_nup82(self):
        """Test generation of comparative models for Nup182"""
        self.run_modeller_script('Nup82', 'Build-n82.py')
        self.check_modeller_model('Nup82', (522, 713))

    def test_cc_tr1(self):
        """Test generation of comparative models for coiled coil complex 1"""
        self.run_modeller_script('cc_tr1', 'cc_tr1.py')
        self.check_modeller_model('cc_tr1', (637, 612))

    def test_cc_tr2(self):
        """Test generation of comparative models for coiled coil complex 2"""
        self.run_modeller_script('cc_tr2', 'cc_tr2.py')
        self.check_modeller_model('cc_tr2', (742, 669))

    def test_cc_tr3(self):
        """Test generation of comparative models for coiled coil complex 3"""
        self.run_modeller_script('cc_tr3', 'cc_tr3.py')
        self.check_modeller_model('cc_tr3', (788, 713))


if __name__ == '__main__':
    unittest.main()
