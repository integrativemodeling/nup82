from modeller import *
from modeller.automodel import *

log.verbose()
env = environ()
env.io.atom_files_directory = ['../../../']

class MyModel(automodel):
    def special_patches(self, aln):
        self.rename_segments(('C', 'D', 'E'), (742, 1332, 625))

a = MyModel(env,
            alnfile='cc_tr2.ali',
            knowns='5cws',
            sequence='cc_tr2')

a.starting_model = 1
a.ending_model = 1
a.make()

