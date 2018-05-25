from modeller import *
from modeller.automodel import *

log.verbose()
env = environ()
env.io.atom_files_directory = ['../../../']

class MyModel(automodel):
    def special_patches(self, aln):
        self.rename_segments(('C', 'D', 'E'), (637, 1211, 522))

a = MyModel(env,
            alnfile='cc_tr1.ali',
            knowns='5cws',
            sequence='cc_tr1')

a.starting_model = 1
a.ending_model = 1
a.make()

