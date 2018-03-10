from modeller import *
from modeller.automodel import *
from modeller.scripts import complete_pdb

import fnmatch
import os

log.verbose()
env1 = environ()
env2 = environ()
env3 = environ()

#It begins from here
aln1 = alignment(env1)
aln1.append(file='Nsp1a.ali')
aln1.write(file='Nsp1a_f.ali', alignment_format='PIR')
aln1.write(file='Nsp1a_f.pap', alignment_format='PAP')
aln1.id_table(matrix_file='Nsp1a_f.mat')
aln1.check()

aln2 = alignment(env2)
aln2.append(file='Nsp1b.ali')
aln2.write(file='Nsp1b_f.ali', alignment_format='PIR')
aln2.write(file='Nsp1b_f.pap', alignment_format='PAP')
aln2.id_table(matrix_file='Nsp1b_f.mat')
aln2.check()

aln3 = alignment(env3)
aln3.append(file='Nsp1c.ali')
aln3.write(file='Nsp1c_f.ali', alignment_format='PIR')
aln3.write(file='Nsp1c_f.pap', alignment_format='PAP')
aln3.id_table(matrix_file='Nsp1c_f.mat')
aln3.check()

class MyModel1(automodel):
    def special_patches(self, aln1):
        # Rename both chains and renumber the residues in each
        self.rename_segments(segment_ids=['C'],
                             renumber_residues=[637])
class MyModel2(automodel):
    def special_patches(self, aln2):
        self.rename_segments(segment_ids=['C'],
                             renumber_residues=[742])
class MyModel3(automodel):
    def special_patches(self, aln3):
        self.rename_segments(segment_ids=['C'],
                             renumber_residues=[788])


a = MyModel1(env1,
            alnfile='Nsp1a_f.ali',
            knowns=('5cws_C1'),  
            sequence='NS1',
            assess_methods=(assess.DOPE, assess.GA341))

a.starting_model = 1
a.ending_model = 1

a.make()

b = MyModel2(env2,
            alnfile='Nsp1b_f.ali',
            knowns=('5cws_C2'),
            sequence='NS2',
            assess_methods=(assess.DOPE, assess.GA341))

b.starting_model = 1
b.ending_model = 1

b.make()

c = MyModel3(env3,
            alnfile='Nsp1c_f.ali',
            knowns=('5cws_C3'),
            sequence='NS3',
            assess_methods=(assess.DOPE, assess.GA341))

c.starting_model = 1
c.ending_model = 1

c.make()
