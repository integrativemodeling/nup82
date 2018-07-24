#!/usr/bin/env python

import unittest
import os
import sys
import subprocess
import ihm.reader

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

    def test_mmcif(self):
        """Test generation of mmCIF output"""
        os.chdir(os.path.join(TOPDIR, 'template'))
        if os.path.exists("nup82.cif"):
            os.unlink("nup82.cif")
        # Potentially override methods that need network access
        env = os.environ.copy()
        env['PYTHONPATH'] = os.path.join(TOPDIR, 'test', 'mock') \
                            + ':' + env.get('PYTHONPATH', '')
        p = subprocess.check_call(["python", "1_modeling_initial_random.py",
                     "--dry-run", "-em2d", "../data/em2d/2.pgm",
                     "-weight", "10000.0", "--mmcif=nup82.cif"], env=env)
        # Check output file
        self._check_mmcif_file('nup82.cif')

    def _check_mmcif_file(self, fname):
        with open(fname) as fh:
            s, = ihm.reader.read(fh)
        self.assertEqual(len(s.citations), 1)
        self.assertEqual(s.citations[0].doi, '10.1016/j.cell.2016.10.028')
        self.assertEqual(len(s.software), 9)
        self.assertEqual(len(s.orphan_starting_models), 30)
        # Should be a single state, of one model
        self.assertEqual(len(s.state_groups), 1)
        self.assertEqual(len(s.state_groups[0]), 1)
        self.assertEqual(len(s.state_groups[0][0]), 1)
        model = s.state_groups[0][0][0][0]
        self.assertEqual(len(model._atoms), 0)
        self.assertEqual(len(model._spheres), 3274)
        # Should be 1 ensemble (cluster)
        self.assertEqual([e.num_models for e in s.ensembles], [370])
        # Check localization densities
        self.assertEqual([len(e.densities) for e in s.ensembles], [10])
        self.assertEqual([len(e.sequence) for e in s.entities],
                         [92, 713, 1460, 823, 1113])
        self.assertEqual([a.details for a in s.asym_units],
                         ['Dyn2.1', 'Dyn2.2', 'Nup82.1', 'Nup82.2', 'Nup159.1',
                          'Nup159.2', 'Nsp1.1', 'Nsp1.2', 'Nup116.1',
                          'Nup116.2'])
        # 27 restraints - 3 crosslinks, 21 EM2D images, 3 SAXS profiles
        self.assertEqual(len(s.restraints), 27)
        xl1, xl2, xl3 = s.restraints[:3]
        self.assertEqual(xl1.linker_type, 'DSS')
        self.assertEqual(len(xl1.experimental_cross_links), 240)
        self.assertEqual(len(xl1.cross_links), 942)
        self.assertEqual(xl1.dataset.location.path,
                'data/XL_wtNup82_DSS_standardized_no_FG_2copies_Ambiguity3.csv')
        # No final psi/sigma values available
        self.assertEqual(sum(len(x.fits) for x in xl1.cross_links), 0)

        self.assertEqual(xl2.linker_type, 'DSS')
        self.assertEqual(xl3.linker_type, 'EDC')

        em2d_rsrs = s.restraints[3:-3]
        for em2d in em2d_rsrs:
            # We only know the resolution
            self.assertAlmostEqual(em2d.image_resolution, 35.0, places=1)
            self.assertEqual(em2d.number_raw_micrographs, None)
            self.assertEqual(len(em2d.fits), 0)

        sas_rsrs = s.restraints[-3:]
        for sas in sas_rsrs:
            self.assertEqual(sas.segment, False)
            self.assertEqual(sas.fitting_method, "FoXS")
            self.assertEqual(sas.multi_state, False)
            self.assertEqual(sas.number_of_gaussians, None)

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
