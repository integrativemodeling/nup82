from math import sqrt

from Bio.PDB import *
#from Bio.PDB.PDBParser import PDBParser
import argparse

parser = argparse.ArgumentParser(description='Extract PDB chains')
parser.add_argument('-pdb', action="store", dest="pdbfile", help="the pdb file" )

parser=parser.parse_args()

pdbparser = PDBParser()

structure = pdbparser.get_structure(parser.pdbfile,parser.pdbfile)
model = structure[0]
chainA = model['A']
chainB = model['B']
chainC = model['C']
chainD = model['D']
chainE = model['E']
chainF = model['F']
chainG = model['G']
chainH = model['H']
chainI = model['I']
chainJ = model['J']
chainK = model['K']
chainL = model['L']
chainM = model['M']
chainN = model['N']
chainO = model['O']
chainP = model['P']
chainQ = model['Q']

class chainA_Select(Select):
    def accept_chain(self, chain):
        if chain.get_id()=='A':
            return 1
        else:
            return 0

class chainB_Select(Select):
    def accept_chain(self, chain):
        if chain.get_id()=='B':
            return 1
        else:
            return 0

class chainC_Select(Select):
    def accept_chain(self, chain):
        if chain.get_id()=='C':
            return 1
        else:
            return 0

class chainD_Select(Select):
    def accept_chain(self, chain):
        if chain.get_id()=='D':
            return 1
        else:
            return 0

class chainE_Select(Select):
    def accept_chain(self, chain):
        if chain.get_id()=='E':
            return 1
        else:
            return 0

class chainF_Select(Select):
    def accept_chain(self, chain):
        if chain.get_id()=='F':
            return 1
        else:
            return 0

class chainG_Select(Select):
    def accept_chain(self, chain):
        if chain.get_id()=='G':
            return 1
        else:
            return 0

class chainH_Select(Select):
    def accept_chain(self, chain):
        if chain.get_id()=='H':
            return 1
        else:
            return 0

class chainI_Select(Select):
    def accept_chain(self, chain):
        if chain.get_id()=='I':
            return 1
        else:
            return 0

class chainJ_Select(Select):
    def accept_chain(self, chain):
        if chain.get_id()=='J':
            return 1
        else:
            return 0

class chainK_Select(Select):
    def accept_chain(self, chain):
        if chain.get_id()=='K':
            return 1
        else:
            return 0

class chainL_Select(Select):
    def accept_chain(self, chain):
        if chain.get_id()=='L':
            return 1
        else:
            return 0

class chainM_Select(Select):
    def accept_chain(self, chain):
        if chain.get_id()=='M':
            return 1
        else:
            return 0

class chainN_Select(Select):
    def accept_chain(self, chain):
        if chain.get_id()=='N':
            return 1
        else:
            return 0

class chainO_Select(Select):
    def accept_chain(self, chain):
        if chain.get_id()=='O':
            return 1
        else:
            return 0

class chainP_Select(Select):
    def accept_chain(self, chain):
        if chain.get_id()=='P':
            return 1
        else:
            return 0

class chainQ_Select(Select):
    def accept_chain(self, chain):
        if chain.get_id()=='Q':
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
io.save('chainI.pdb', chainI_Select())
io.save('chainJ.pdb', chainJ_Select())
io.save('chainK.pdb', chainK_Select())
io.save('chainL.pdb', chainL_Select())
io.save('chainM.pdb', chainM_Select())
io.save('chainN.pdb', chainN_Select())
io.save('chainO.pdb', chainO_Select())
io.save('chainP.pdb', chainP_Select())
io.save('chainQ.pdb', chainQ_Select())

