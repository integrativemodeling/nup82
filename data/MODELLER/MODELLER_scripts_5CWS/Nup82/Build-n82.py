from modeller import *
from modeller.automodel import *
from modeller.scripts import complete_pdb
import sys

import fnmatch
import os

log.verbose()
env = environ()

#It begins from here
aln = alignment(env)
aln.append(file='Nup82_long.ali')
#aln.append(file='Nup82_short.ali')
aln.write(file='Nup82_f.ali', alignment_format='PIR')
aln.write(file='Nup82_f.pap', alignment_format='PAP')
aln.id_table(matrix_file='Nup82_f.mat')
aln.check()

class MyModel(automodel):
    def special_patches(self, aln):
        # Rename both chains and renumber the residues in each
        self.rename_segments(segment_ids=['E'],
                             renumber_residues=[522])
    def special_restraints(self, aln):
        rsr = self.restraints
        at = self.atoms

        rsr.add(secondary_structure.alpha(self.residue_range('522:E', '560:E')))


a = MyModel(env,
            alnfile='Nup82_f.ali',
            knowns=('5cws_E'),  
            sequence='Nup82',
            assess_methods=(assess.DOPE, assess.GA341))

a.starting_model = 1
a.ending_model = 20
if '--test' in sys.argv:
    a.ending_model = 1

a.make()
