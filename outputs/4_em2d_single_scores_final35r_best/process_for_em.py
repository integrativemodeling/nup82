from math import sqrt
from Bio.PDB import *
#from Bio.PDB.PDBParser import PDBParser
import argparse
import os

parser = argparse.ArgumentParser(description='Extract PDB chains')
parser.add_argument('-pdb', action="store", dest="pdbfile", help="the input pdb file" )
parser.add_argument('-out', action="store", dest="outfile", help="the output pdb file" )
parser=parser.parse_args()

pdbparser = PDBParser()

structure = pdbparser.get_structure(parser.pdbfile,parser.pdbfile)
model = structure[0]
chainA = model['A']
#chainB = model['B']
#chainC = model['C']
#chainD = model['D']

# http://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ 

# Dyn2.1
class chainA_Select(Select):
    def accept_chain(self, chain):
        if chain.get_id()=='A':
            return 1
        else:
            return 0

# Dyn2.2
class chainB_Select(Select):
    def accept_chain(self, chain):
        if chain.get_id()=='B':
            return 1
        else:
            return 0

# Nup82.1
class chainC_Select(Select):
    def accept_chain(self, chain):
        if chain.get_id()=='C':
            return 1
        else:
            return 0

# Nup82.2
class chainD_Select(Select):
    def accept_chain(self, chain):
        if chain.get_id()=='D':
            return 1
        else:
            return 0

# Nup159.1
class chainE_Select(Select):
    def accept_chain(self, chain):
        if chain.get_id()=='E':
            return 1
        else:
            return 0
    
    def accept_residue(self, residue):
        #if residue.get_resname()=='THR':
        if residue.get_id()[1]  > 1116:
            return 1
        else:
            return 0

# Nup159.2
class chainF_Select(Select):
    def accept_chain(self, chain):
        if chain.get_id()=='F':
            return 1
        else:
            return 0

    def accept_residue(self, residue):
        #if residue.get_resname()=='THR':
        if residue.get_id()[1]  > 1116:
            return 1
        else:
            return 0

# Nsp1.1
class chainG_Select(Select):
    def accept_chain(self, chain):
        if chain.get_id()=='G':
            return 1
        else:
            return 0

    def accept_residue(self, residue):
        #if residue.get_resname()=='THR':
        if residue.get_id()[1]  > 636:
            return 1
        else:
            return 0

# Nsp1.2
class chainH_Select(Select):
    def accept_chain(self, chain):
        if chain.get_id()=='H':
            return 1
        else:
            return 0

    def accept_residue(self, residue):
        #if residue.get_resname()=='THR':
        if residue.get_id()[1]  > 636:
            return 1
        else:
            return 0

io = PDBIO()
io.set_structure(structure)
io.save('chainA.pdb', chainA_Select())
io.save('chainB.pdb', chainB_Select())
io.save('chainC.pdb', chainC_Select())
io.save('chainD.pdb', chainD_Select())
io.save('chainE.pdb', chainE_Select())
io.save('chainF.pdb', chainF_Select())
io.save('chainG.pdb', chainG_Select())
io.save('chainH.pdb', chainH_Select())

structureA = pdbparser.get_structure('chainA.pdb', 'chainA.pdb')
structureB = pdbparser.get_structure('chainB.pdb', 'chainB.pdb')
structureC = pdbparser.get_structure('chainC.pdb', 'chainC.pdb')
structureD = pdbparser.get_structure('chainD.pdb', 'chainD.pdb')
structureE = pdbparser.get_structure('chainE.pdb', 'chainE.pdb')
structureF = pdbparser.get_structure('chainF.pdb', 'chainF.pdb')
structureG = pdbparser.get_structure('chainG.pdb', 'chainG.pdb')
structureH = pdbparser.get_structure('chainH.pdb', 'chainH.pdb')

outfile=open(parser.outfile, "w")
io = PDBIO(1)
io.set_structure(structureA)
io.save(outfile, write_end=True)
io.set_structure(structureB)
io.save(outfile, write_end=True)
io.set_structure(structureC)
io.save(outfile, write_end=True)
io.set_structure(structureD)
io.save(outfile, write_end=True)
io.set_structure(structureE)
io.save(outfile, write_end=True)
io.set_structure(structureF)
io.save(outfile, write_end=True)
io.set_structure(structureG)
io.save(outfile, write_end=True)
io.set_structure(structureH)
io.save(outfile, write_end=True)
outfile.close()


