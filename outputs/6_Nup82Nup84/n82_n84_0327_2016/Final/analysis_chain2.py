from modeller import *
from modeller.automodel import *
from modeller.scripts import complete_pdb
import fnmatch
import os

log.verbose()
env = environ()


for files in os.listdir('.'):
    if fnmatch.fnmatch(files, '*_A.rebuilt.pdb'):
        print files
        mdl = model(env, file=files)
        mdl.chains[0].name = 'A'
        mdl.write(files)

for files in os.listdir('.'):
    if fnmatch.fnmatch(files, '*_B.rebuilt.pdb'):
        print files
        mdl = model(env, file=files)
        mdl.chains[0].name = 'B'
        mdl.write(files)

for files in os.listdir('.'):
    if fnmatch.fnmatch(files, '*_C.rebuilt.pdb'):
        print files
        mdl = model(env, file=files)
        mdl.chains[0].name = 'C'
        mdl.write(files)

for files in os.listdir('.'):
    if fnmatch.fnmatch(files, '*_D.rebuilt.pdb'):
        print files
        mdl = model(env, file=files)
        mdl.chains[0].name = 'D'
        mdl.write(files)

for files in os.listdir('.'):
    if fnmatch.fnmatch(files, '*_E.rebuilt.pdb'):
        print files
        mdl = model(env, file=files)
        mdl.chains[0].name = 'E'
        mdl.write(files)

for files in os.listdir('.'):
    if fnmatch.fnmatch(files, '*_F.rebuilt.pdb'):
        print files
        mdl = model(env, file=files)
        mdl.chains[0].name = 'F'
        mdl.write(files)

for files in os.listdir('.'):
    if fnmatch.fnmatch(files, '*_G.rebuilt.pdb'):
        print files
        mdl = model(env, file=files)
        mdl.chains[0].name = 'G'
        mdl.write(files)

for files in os.listdir('.'):
    if fnmatch.fnmatch(files, '*_H.rebuilt.pdb'):
        print files
        mdl = model(env, file=files)
        mdl.chains[0].name = 'H'
        mdl.write(files)

for files in os.listdir('.'):
    if fnmatch.fnmatch(files, '*_I.rebuilt.pdb'):
        print files
        mdl = model(env, file=files)
        mdl.chains[0].name = 'I'
        mdl.write(files)

for files in os.listdir('.'):
    if fnmatch.fnmatch(files, '*_J.rebuilt.pdb'):
        print files
        mdl = model(env, file=files)
        mdl.chains[0].name = 'J'
        mdl.write(files)

for files in os.listdir('.'):
    if fnmatch.fnmatch(files, '*_K.rebuilt.pdb'):
        print files
        mdl = model(env, file=files)
        mdl.chains[0].name = 'K'
        mdl.write(files)

for files in os.listdir('.'):
    if fnmatch.fnmatch(files, '*_L.rebuilt.pdb'):
        print files
        mdl = model(env, file=files)
        mdl.chains[0].name = 'L'
        mdl.write(files)

for files in os.listdir('.'):
    if fnmatch.fnmatch(files, '*_M.rebuilt.pdb'):
        print files
        mdl = model(env, file=files)
        mdl.chains[0].name = 'M'
        mdl.write(files)

for files in os.listdir('.'):
    if fnmatch.fnmatch(files, '*_N.rebuilt.pdb'):
        print files
        mdl = model(env, file=files)
        mdl.chains[0].name = 'N'
        mdl.write(files)

for files in os.listdir('.'):
    if fnmatch.fnmatch(files, '*_O.rebuilt.pdb'):
        print files
        mdl = model(env, file=files)
        mdl.chains[0].name = 'O'
        mdl.write(files)

for files in os.listdir('.'):
    if fnmatch.fnmatch(files, '*_P.rebuilt.pdb'):
        print files
        mdl = model(env, file=files)
        mdl.chains[0].name = 'P'
        mdl.write(files)

for files in os.listdir('.'):
    if fnmatch.fnmatch(files, '*_Q.rebuilt.pdb'):
        print files
        mdl = model(env, file=files)
        mdl.chains[0].name = 'Q'
        mdl.write(files)

