from modeller import *
from modeller.automodel import *

log.verbose()
env = environ()
env.io.atom_files_directory = ['../../../']

class MyModel(automodel):
    def special_patches(self, aln):
        self.rename_segments(('C', 'D', 'E'), (788, 1382, 678))

a = MyModel(env,
            alnfile='cc_tr3.ali',
            knowns='5cws',
            sequence='cc_tr3')

a.starting_model = 1
a.ending_model = 1
a.make()

