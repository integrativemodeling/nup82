from modeller import *
from modeller.automodel import *
from modeller.scripts import complete_pdb
import fnmatch
import os

log.verbose()
env = environ()


for files in os.listdir('.'):
    if fnmatch.fnmatch(files, 'Nup84.pdb'):
        print files
        mdl = model(env, file=files)
        mdl.chains[0].name = 'K'
        mdl.write(files)

