#!/usr/bin/env python
#####################################################
# Last Update: September 10th, 2015
# by Seung Joong Kim and Riccardo Pellarin
# at Andrej Sali group, University of California San Francisco (UCSF)
#####################################################
import IMP
import IMP.core
#import IMP.base
import IMP.algebra
import IMP.atom
import IMP.container

#import crosslinking_nup82
import IMP.pmi.restraints.crosslinking
import IMP.pmi.restraints.stereochemistry
#import IMP.pmi.restraints.em
import em2d_nup82
#import IMP.pmi.restraints.em2d
import IMP.pmi.restraints.basic
import IMP.pmi.restraints.proteomics
#import representation_nup82
import IMP.pmi.representation
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.tools
import IMP.pmi.output
import IMP.pmi.samplers
import random

import os


#####################################################
# Parsing parameter inputs
#####################################################
import argparse

parser = argparse.ArgumentParser(description='Performing the INITIAL/REFINEMENT Monte Carlo job, with crosslinks and selected/ALL domain mapping data. Example of usage: setup_environment.sh python ./sj_SEA_XLDM.py -f models_1877.rmf -n 0')
parser.add_argument('-copy', action="store", dest="ncopy", help="copy numbers (stoichiometry) for SEA4 and Seh1" )
parser.add_argument('-sym', action="store", dest="symmetry", help="symmetry option for SEA4 and Seh1" )
parser.add_argument('-rmf', action="store", dest="rmf_input", help="rmf file name to continue" )
parser.add_argument('-rmf_n', action="store", dest="rmf_frame_number", help="rmf frame number to continue" )
parser.add_argument('-em2d', action="store", dest="em2d_input", help="em2d image file name to read" )
parser.add_argument('-r', action="store", dest="nrepeats", help="number of Monte Carlo cycles" )
parser.add_argument('-x', action="store", dest="XL_input", help="Cross-links file name to read" )
parser.add_argument('-out', action="store", dest="folder_output", help="folder name for output" )
parser.add_argument('-o', action="store", dest="rmf_output", help="rmf file name for output" )
parser.add_argument('-s', action="store", dest="stat_output", help="stat file name for output" )
parser.add_argument('-REFINE', action="store", dest="refinement", help="refinement True or False" )
parser.add_argument('-weight', action="store", dest="weight", help="weight for the EM 2D restraint" )
parser.add_argument('-res_cry', action="store", dest="res_cry", help="resolution of the crystal structures" )
parser.add_argument('-res_hom', action="store", dest="res_hom", help="resolution of the comparative (homology) models" )
parser.add_argument('-res_ev', action="store", dest="res_ev", help="resolution of the excluded volume restraints" )
parser.add_argument('-res_compo', action="store", dest="res_compo", help="resolution of the composite restraints" )
parser.add_argument('-draw_hierarchy', action="store", dest="draw_hierarchy", help="draw hierarchy" )
inputs = parser.parse_args()

#----------------------------------------------------
# Setting up the input parameters
#----------------------------------------------------
if (inputs.ncopy is None) :
    inputs.ncopy = "2"
if (inputs.symmetry == "True") or (inputs.symmetry == "true") or (inputs.symmetry == "Yes") or (inputs.symmetry == "yes") :
    inputs.symmetry = True
else:
    inputs.symmetry = False
    
if (inputs.rmf_input is not None) :
    f=open(inputs.rmf_input,"r")
    f.close()
if (inputs.rmf_frame_number is None) :
    inputs.rmf_frame_number = 0
if (inputs.em2d_input is not None) :
    f=open(inputs.em2d_input,"r")
    f.close()
    
if (inputs.XL_input is None) :
    inputs.XL_input = "./p85_xl1.txt"
else :
    f=open(inputs.XL_input,"r")
    f.close()
if (inputs.nrepeats is None) :
    inputs.nrepeats = 1000
if (inputs.folder_output is None) :
    inputs.folder_output = "output"
if (inputs.rmf_output is None) :
    inputs.rmf_output = "models.rmf"
if (inputs.stat_output is None) :
    inputs.stat_output = "stat.dat"
if (inputs.refinement == "True") or (inputs.refinement == "true") or (inputs.refinement == "Yes") or (inputs.refinement == "yes") :
    inputs.refinement = True
else:
    inputs.refinement = False
if (inputs.weight is None) :
    inputs.weight = 10000.0
    
if (inputs.res_cry is None) :
    inputs.res_cry = 1.0
if (inputs.res_hom is None) :
    inputs.res_hom = 5.0
if (inputs.res_ev is None) :
    inputs.res_ev = 10.0
if (inputs.res_compo is None) :
    inputs.res_compo = 100.0
if (inputs.draw_hierarchy == "True") or (inputs.draw_hierarchy == "true") or (inputs.draw_hierarchy == "Yes") or (inputs.draw_hierarchy == "yes") :
    inputs.draw_hierarchy = True
else:
    inputs.draw_hierarchy = False
print inputs


#####################################################
# setting up topology
#####################################################
m = IMP.Model()
simo = IMP.pmi.representation.Representation(m,upperharmonic=True,disorderedlength=False)
#simo = representation_nup82.Representation(m,upperharmonic=True,disorderedlength=False)


#####################################################
# setting up parameters
#####################################################
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
print "rank = ", rank

# rigid body movement params
rbmaxtrans = 4.00
rbmaxrot = 0.04

# flexible bead movement
fbmaxtrans = 4.00
outputobjects = []
sampleobjects = []

res_cry = int(inputs.res_cry)
res_ev = int(inputs.res_ev)
res_conn = int(inputs.res_cry)
beadsize = 10
beadsize20 = 20
beadsize25 = 25
beadsize100 = 100
em2d_weight = float(inputs.weight)

datadirectory = "../data/"
fasta_files = "../data/protein_fasta."

n84_fastafile  ='../data_nup84/protein_fasta.Nup84.txt'
n85_fastafile  ='../data_nup84/protein_fasta.Nup85.txt'
n120_fastafile ='../data_nup84/protein_fasta.Nup120.txt'
n133_fastafile ='../data_nup84/protein_fasta.Nup133.txt'
n145c_fastafile='../data_nup84/protein_fasta.Nup145c.txt'
seh1_fastafile ='../data_nup84/protein_fasta.Seh1.txt'
sec13_fastafile='../data_nup84/protein_fasta.Sec13.txt'

# After removal of disordered regions in PDB files.
n84_pdbfile    ='../data_nup84/ScNup84_7-488_506-726_new2.pdb'
n85_pdbfile    ='../data_nup84/ScNup85_44-744_new2.pdb'
n120_pdbfile   ='../data_nup84/ScNup120_1-1037_new.pdb'
n133n_pdbfile  ='../data_nup84/ScNup133N_56-480_new.pdb'
n133c_pdbfile  ='../data_nup84/ScNup133C_490_1157_new.pdb'
n145c_pdbfile  ='../data_nup84/ScNup145C_126-553_new2.pdb'
seh1_pdbfile   ='../data_nup84/ScSeh1_1-346_new2.pdb'
sec13_pdbfile  ='../data_nup84/ScSec13_2-296_new2.pdb'


#####################################################
# REPRESENTATION
#####################################################
# compname  hier_name      color   fastafile                  fastaid    pdbname                       chain res_range       read_em_files bead_size  rigid_body super_rigid_body em_num_components em_txt_file_name em_mrc_file_name chain_of_super_rb
domains = \
[("Dyn2.1",  "Dyn2.1",     0.48,   fasta_files+"Dyn2.txt",    "Dyn2",    datadirectory+"4DS1.pdb",     "A",  (  1, 92,0),    None,         beadsize,       1,      [1],        4,               None,            None, [1]),
 ("Dyn2.2",  "Dyn2.2",     0.65,   fasta_files+"Dyn2.txt",    "Dyn2",    datadirectory+"4DS1.pdb",     "C",  (  1, 92,0),    None,         beadsize,       1,      [1],        4,               None,            None, [1]),

 
 ("Nup82.1", "Nup82.1_1",   0.0,   fasta_files+"Nup82.txt",   "Nup82",   datadirectory+"3PBP1.pdb",    "G",  (  1,452,0),    None,         beadsize,       1,      [1],        18,              None,            None, [1]),
 ("Nup82.1", "Nup82.1_11",  0.0,   fasta_files+"Nup82.txt",   "Nup82",   "BEADS",                      " ",  (453,521,0),    None,         beadsize,       1,      [1],        0,               None,            None, [1]),
 ("Nup82.1", "Nup82.1_2",   0.0,   fasta_files+"Nup82.txt",   "Nup82",   datadirectory+"cc_tr1_1.pdb", "E",  (522,612,0),    None,         beadsize,       1,      [1],        4,               None,            None, [1]),
 ("Nup82.1", "Nup82.1_21",  0.0,   fasta_files+"Nup82.txt",   "Nup82",   "BEADS",                      " ",  (613,624,0),    None,         beadsize,       1,      [1],        0,               None,            None, [1]),
 ("Nup82.1", "Nup82.1_3",   0.0,   fasta_files+"Nup82.txt",   "Nup82",   datadirectory+"cc_tr2_1.pdb", "E",  (625,669,0),    None,         beadsize,       1,      [1],        2,               None,            None, [1]),
 ("Nup82.1", "Nup82.1_31",  0.0,   fasta_files+"Nup82.txt",   "Nup82",   "BEADS",                      " ",  (670,677,0),    None,         beadsize,       1,      [1],        0,               None,            None, [1]),
 ("Nup82.1", "Nup82.1_4",   0.0,   fasta_files+"Nup82.txt",   "Nup82",   datadirectory+"cc_tr3_1.pdb", "E",  (678,713,0),    None,         beadsize,       1,      [1],        2,               None,            None, [1]),

 ("Nup82.2", "Nup82.2_1",   0.15,  fasta_files+"Nup82.txt",   "Nup82",   datadirectory+"3PBP2.pdb",    "G",  (  1,452,0),    None,         beadsize,       1,      [1],        18,              None,            None, [1]),
 ("Nup82.2", "Nup82.2_11",  0.15,  fasta_files+"Nup82.txt",   "Nup82",   "BEADS",                      " ",  (453,521,0),    None,         beadsize,       1,      [1],        0,               None,            None, [1]),
 ("Nup82.2", "Nup82.2_2",   0.15,  fasta_files+"Nup82.txt",   "Nup82",   datadirectory+"cc_tr1_2.pdb", "E",  (522,612,0),    None,         beadsize,       1,      [1],        4,               None,            None, [1]),
 ("Nup82.2", "Nup82.2_21",  0.15,  fasta_files+"Nup82.txt",   "Nup82",   "BEADS",                      " ",  (613,624,0),    None,         beadsize,       1,      [1],        0,               None,            None, [1]),
 ("Nup82.2", "Nup82.2_3",   0.15,  fasta_files+"Nup82.txt",   "Nup82",   datadirectory+"cc_tr2_2.pdb", "E",  (625,669,0),    None,         beadsize,       1,      [1],        2,               None,            None, [1]),
 ("Nup82.2", "Nup82.2_31",  0.15,  fasta_files+"Nup82.txt",   "Nup82",   "BEADS",                      " ",  (670,677,0),    None,         beadsize,       1,      [1],        0,               None,            None, [1]),
 ("Nup82.2", "Nup82.2_4",   0.15,  fasta_files+"Nup82.txt",   "Nup82",   datadirectory+"cc_tr3_2.pdb", "E",  (678,713,0),    None,         beadsize,       1,      [1],        2,               None,            None, [1]),


 ("Nup159.1", "Nup159.1_1", 1.0,   fasta_files+"Nup159.txt",  "Nup159",  datadirectory+"1XIP1.pdb",    "A",  (   1, 381,0),  None,         beadsize,       2,      [2],        0,               None,            None, [1]),
 ("Nup159.1", "Nup159.1_11",1.0,   fasta_files+"Nup159.txt",  "Nup159",  "BEADS",                      " ",  ( 382,1116,0),  None,         beadsize100,    2,      [2],        0,               None,            None, [1]),
 ("Nup159.1", "Nup159.1_2", 1.0,   fasta_files+"Nup159.txt",  "Nup159",  datadirectory+"4DS1.pdb",     "B",  (1117,1126,0),  None,         beadsize,       1,      [1],        1,               None,            None, [1]),
 ("Nup159.1", "Nup159.1_21",1.0,   fasta_files+"Nup159.txt",  "Nup159",  "BEADS",                      " ",  (1127,1210,0),  None,         beadsize,       1,      [1],        0,               None,            None, [1]),
 ("Nup159.1", "Nup159.1_3", 1.0,   fasta_files+"Nup159.txt",  "Nup159",  datadirectory+"cc_tr1_1.pdb", "D",  (1211,1321,0),  None,         beadsize,       1,      [1],        4,               None,            None, [1]),
 ("Nup159.1", "Nup159.1_31",1.0,   fasta_files+"Nup159.txt",  "Nup159",  "BEADS",                      " ",  (1322,1331,0),  None,         beadsize,       1,      [1],        0,               None,            None, [1]),
 ("Nup159.1", "Nup159.1_4", 1.0,   fasta_files+"Nup159.txt",  "Nup159",  datadirectory+"cc_tr2_1.pdb", "D",  (1332,1372,0),  None,         beadsize,       1,      [1],        2,               None,            None, [1]),
 ("Nup159.1", "Nup159.1_41",1.0,   fasta_files+"Nup159.txt",  "Nup159",  "BEADS",                      " ",  (1373,1381,0),  None,         beadsize,       1,      [1],        0,               None,            None, [1]),
 ("Nup159.1", "Nup159.1_5", 1.0,   fasta_files+"Nup159.txt",  "Nup159",  datadirectory+"cc_tr3_1.pdb", "D",  (1382,1412,0),  None,         beadsize,       1,      [1],        2,               None,            None, [1]),
 ("Nup159.1", "Nup159.1_51",1.0,   fasta_files+"Nup159.txt",  "Nup159",  "BEADS",                      " ",  (1413,1428,0),  None,         beadsize,       1,      [1],        0,               None,            None, [1]),
 ("Nup159.1", "Nup159.1_6", 1.0,   fasta_files+"Nup159.txt",  "Nup159",  datadirectory+"3PBP1.pdb",    "I",  (1429,1460,0),  None,         beadsize,       1,      [1],        2,               None,            None, [1]),

 ("Nup159.2", "Nup159.2_1", 0.9,   fasta_files+"Nup159.txt",  "Nup159",  datadirectory+"1XIP2.pdb",    "A",  (   1, 381,0),  None,         beadsize,       3,      [3],        0,               None,            None, [1]),
 ("Nup159.2", "Nup159.2_11",0.9,   fasta_files+"Nup159.txt",  "Nup159",  "BEADS",                      " ",  ( 382,1116,0),  None,         beadsize100,    3,      [3],        0,               None,            None, [1]),
 ("Nup159.2", "Nup159.2_2", 0.9,   fasta_files+"Nup159.txt",  "Nup159",  datadirectory+"4DS1.pdb",     "D",  (1117,1126,0),  None,         beadsize,       1,      [1],        1,               None,            None, [1]),
 ("Nup159.2", "Nup159.2_21",0.9,   fasta_files+"Nup159.txt",  "Nup159",  "BEADS",                      " ",  (1127,1210,0),  None,         beadsize,       1,      [1],        0,               None,            None, [1]),
 ("Nup159.2", "Nup159.2_3", 0.9,   fasta_files+"Nup159.txt",  "Nup159",  datadirectory+"cc_tr1_2.pdb", "D",  (1211,1321,0),  None,         beadsize,       1,      [1],        4,               None,            None, [1]),
 ("Nup159.2", "Nup159.2_31",0.9,   fasta_files+"Nup159.txt",  "Nup159",  "BEADS",                      " ",  (1322,1331,0),  None,         beadsize,       1,      [1],        0,               None,            None, [1]),
 ("Nup159.2", "Nup159.2_4", 0.9,   fasta_files+"Nup159.txt",  "Nup159",  datadirectory+"cc_tr2_2.pdb", "D",  (1332,1372,0),  None,         beadsize,       1,      [1],        2,               None,            None, [1]),
 ("Nup159.2", "Nup159.2_41",0.9,   fasta_files+"Nup159.txt",  "Nup159",  "BEADS",                      " ",  (1373,1381,0),  None,         beadsize,       1,      [1],        0,               None,            None, [1]),
 ("Nup159.2", "Nup159.2_5", 0.9,   fasta_files+"Nup159.txt",  "Nup159",  datadirectory+"cc_tr3_2.pdb", "D",  (1382,1412,0),  None,         beadsize,       1,      [1],        2,               None,            None, [1]),
 ("Nup159.2", "Nup159.2_51",0.9,   fasta_files+"Nup159.txt",  "Nup159",  "BEADS",                      " ",  (1413,1428,0),  None,         beadsize,       1,      [1],        0,               None,            None, [1]),
 ("Nup159.2", "Nup159.2_6", 0.9,   fasta_files+"Nup159.txt",  "Nup159",  datadirectory+"3PBP2.pdb",    "I",  (1429,1460,0),  None,         beadsize,       1,      [1],        2,               None,            None, [1]),


 ("Nsp1.1",  "Nsp1.1_10",   0.3,   fasta_files+"Nsp1.txt",    "Nsp1",    "BEADS",                      " ",  (  1,636,0),    None,         beadsize100,    4,      [4],        0,               None,            None, [1]),
 ("Nsp1.1",  "Nsp1.1_1",    0.3,   fasta_files+"Nsp1.txt",    "Nsp1",    datadirectory+"cc_tr1_1.pdb", "C",  (637,727,0),    None,         beadsize,       1,      [1],        4,               None,            None, [1]),
 ("Nsp1.1",  "Nsp1.1_11",   0.3,   fasta_files+"Nsp1.txt",    "Nsp1",    "BEADS",                      " ",  (728,741,0),    None,         beadsize,       1,      [1],        0,               None,            None, [1]),
 ("Nsp1.1",  "Nsp1.1_2",    0.3,   fasta_files+"Nsp1.txt",    "Nsp1",    datadirectory+"cc_tr2_1.pdb", "C",  (742,778,0),    None,         beadsize,       1,      [1],        2,               None,            None, [1]),
 ("Nsp1.1",  "Nsp1.1_21",   0.3,   fasta_files+"Nsp1.txt",    "Nsp1",    "BEADS",                      " ",  (779,787,0),    None,         beadsize,       1,      [1],        0,               None,            None, [1]),
 ("Nsp1.1",  "Nsp1.1_3",    0.3,   fasta_files+"Nsp1.txt",    "Nsp1",    datadirectory+"cc_tr3_1.pdb", "C",  (788,823,0),    None,         beadsize,       1,      [1],        2,               None,            None, [1]),
 
 ("Nsp1.2",  "Nsp1.2_10",   0.38,  fasta_files+"Nsp1.txt",    "Nsp1",    "BEADS",                      " ",  (  1,636,0),    None,         beadsize100,    5,      [5],        0,               None,            None, [1]),
 ("Nsp1.2",  "Nsp1.2_1",    0.38,  fasta_files+"Nsp1.txt",    "Nsp1",    datadirectory+"cc_tr1_2.pdb", "C",  (637,727,0),    None,         beadsize,       1,      [1],        4,               None,            None, [1]),
 ("Nsp1.2",  "Nsp1.2_11",   0.38,  fasta_files+"Nsp1.txt",    "Nsp1",    "BEADS",                      " ",  (728,741,0),    None,         beadsize,       1,      [1],        0,               None,            None, [1]),
 ("Nsp1.2",  "Nsp1.2_2",    0.38,  fasta_files+"Nsp1.txt",    "Nsp1",    datadirectory+"cc_tr2_2.pdb", "C",  (742,778,0),    None,         beadsize,       1,      [1],        2,               None,            None, [1]),
 ("Nsp1.2",  "Nsp1.2_21",   0.38,  fasta_files+"Nsp1.txt",    "Nsp1",    "BEADS",                      " ",  (779,787,0),    None,         beadsize,       1,      [1],        0,               None,            None, [1]),
 ("Nsp1.2",  "Nsp1.2_3",    0.38,  fasta_files+"Nsp1.txt",    "Nsp1",    datadirectory+"cc_tr3_2.pdb", "C",  (788,823,0),    None,         beadsize,       1,      [1],        2,               None,            None, [1]),


 ("Nup116.1","Nup116.1_10", 0.75,  fasta_files+"Nup116.txt",  "Nup116",  "BEADS",                      " ",  (  1, 750,0),   None,         beadsize100,    6,      [6],        0,               None,            None, [1]),
 ("Nup116.1","Nup116.1_1",  0.75,  fasta_files+"Nup116.txt",  "Nup116",  datadirectory+"3PBP1.pdb",    "H",  (751,1113,0),   None,         beadsize25,     1,      [1],        0,               None,            None, [1]),

 ("Nup116.2","Nup116.2_10", 0.8,   fasta_files+"Nup116.txt",  "Nup116",  "BEADS",                      " ",  (  1, 750,0),   None,         beadsize100,    7,      [7],        0,               None,            None, [1]),
 ("Nup116.2","Nup116.2_1",  0.8,   fasta_files+"Nup116.txt",  "Nup116",  datadirectory+"3PBP2.pdb",    "H",  (751,1113,0),   None,         beadsize25,     1,      [1],        0,               None,            None, [1]),

 
 ("Nup84",    "Nup84",      0.0,   n84_fastafile,             "Nup84",   n84_pdbfile,                  "A",  (  1, 726,0),   True,         beadsize20,     0,     [0],            7,               None,            None, [0]),
 ("Nup85",    "Nup85_1",    0.2,   n85_fastafile,             "Nup85",   n85_pdbfile,                  "B",  (  1, 529,0),   True,         beadsize20,     0,     [0],            8,               None,            None, [0]),
 ("Nup85",    "Nup85_2",    0.2,   n85_fastafile,             "Nup85",   n85_pdbfile,                  "B",  (530, 744,0),   True,         beadsize20,     0,     [0],            3,               None,            None, [0]),
 ("Nup120",   "Nup120_1",   0.35,  n120_fastafile,            "Nup120",  n120_pdbfile,                 "C",  (  1, 710,0),   True,         beadsize20,     0,     [0],            7,               None,            None, [0]),
 ("Nup120",   "Nup120_2",   0.35,  n120_fastafile,            "Nup120",  n120_pdbfile,                 "C",  (711,1037,0),   True,         beadsize20,     0,     [0],            3,               None,            None, [0]),
 ("Nup133",   "Nup133_1",   0.5,   n133_fastafile,            "Nup133",  n133n_pdbfile,                "D",  (  1, 480,0),   True,         beadsize20,     0,     [0],            5,               None,            None, [0]),
 ("Nup133",   "Nup133_2",   0.5,   n133_fastafile,            "Nup133",  n133c_pdbfile,                "D",  (481,1157,0),   True,         beadsize20,     0,     [0],            7,               None,            None, [0]),
 ("Nup145c",  "Nup145c",    0.65,  n145c_fastafile,           "Nup145c", n145c_pdbfile,                "E",  (  1, 712,0),   True,         beadsize20,     0,     [0],            7,               None,            None, [0]),
 ("Seh1",     "Seh1",       0.8,   seh1_fastafile,            "Seh1",    seh1_pdbfile,                 "F",  (  1, 349,0),   True,         beadsize20,     0,     [0],            4,               None,            None, [0]),
 ("Sec13",    "Sec13",      0.95,  sec13_fastafile,           "Sec13",   sec13_pdbfile,                "G",  (  1, 297,0),   True,         beadsize20,     0,     [0],            3,               None,            None, [0])
]

bm1 = IMP.pmi.macros.BuildModel1(simo)
#bm1.set_gmm_models_directory(datadirectory + "em_gmm_model/")

if (inputs.rmf_input is not None) :
    n82=['Dyn2.1', 'Dyn2.2', 'Nup82.1', 'Nup82.2', 'Nup159.1', 'Nup159.2', 'Nsp1.1', 'Nsp1.2', 'Nup116.1', 'Nup116.2']
    n84=['Nup84', 'Nup85', 'Nup120', 'Nup133', 'Nup145c', 'Seh1', 'Sec13']
    
    #n82=set([s[0] for s in domains])
    for d in list(n82):
        bm1.set_rmf_file(d, inputs.rmf_input, int(inputs.rmf_frame_number))
        
    for d in list(n84):
        bm1.set_rmf_file(d, "../data_nup84/rmfs_nup84/nup84_16.0.rmf3", 0)

bm1.build_model(data_structure = domains, sequence_connectivity_scale=2.0, sequence_connectivity_resolution=1.0)
#bm1.scale_bead_radii(40,0.8)
#resdensities = bm1.get_density_hierarchies([t[1] for t in domainxl_cliques_psi = 0.25s])
#print resdensities; exit()

#model_ps = []
#for h in self.densities:
#    model_ps += IMP.atom.get_leaves(h)


#####################################################
# randomize the initial configuration
#####################################################
if (inputs.rmf_input is None) :
    simo.shuffle_configuration(50)
    #simo.shuffle_configuration(100)


#####################################################
# defines the movers
#####################################################
simo.set_rigid_bodies_max_rot(rbmaxrot)
simo.set_floppy_bodies_max_trans(fbmaxtrans)
simo.set_rigid_bodies_max_trans(rbmaxtrans)

### These two below are already executed in BuildModel1
#simo.set_floppy_bodies()
#simo.setup_bonds()

#prot = simo.prot
outputobjects.append(simo)
sampleobjects.append(simo)


#####################################################
# Restraints setup
# Excluded Volume restraint
#####################################################
ev = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(simo, resolution = res_ev)
ev.add_to_model()
outputobjects.append(ev)
print(ev.get_output())
print "ExcludedVolumeSphere !!\n"


#####################################################
# Restraints setup
# External Barrier restraint
#####################################################
eb = IMP.pmi.restraints.basic.ExternalBarrier(simo, radius = 500)
eb.add_to_model()
outputobjects.append(eb)
print(eb.get_output())
print "ExternalBarrier !!\n"


#####################################################
# Restraints setup
# Cross-link restraints for the Nup84 complex
#####################################################
columnmap={}
columnmap["Protein1"]=0
columnmap["Protein2"]=2
columnmap["Residue1"]=1
columnmap["Residue2"]=3
columnmap["IDScore"]=4
columnmap["XLUniqueID"]=5

ids_map=IMP.pmi.tools.map()
ids_map.set_map_element(1.0,1.0)

if (True):
    xl1 = IMP.pmi.restraints.crosslinking.ISDCrossLinkMS(simo,
                                       '../data_nup84/yeast_Nup84_DSS.new.dat',
                                       length=21.0,
                                       slope=0.01,
                                       columnmapping=columnmap,
                                       ids_map=ids_map,resolution=1.0,
                                       filelabel="DSS",
                                       label="DSS")
    xl1.add_to_model()
    sampleobjects.append(xl1)
    outputobjects.append(xl1)
    xl1.set_psi_is_sampled(False)
    psi=xl1.get_psi(1.0)[0]
    psi.set_scale(0.05)


    xl2 = IMP.pmi.restraints.crosslinking.ISDCrossLinkMS(simo,
                                       '../data_nup84/EDC_XL_122013.new.dat',
                                       length=16.0,
                                       slope=0.01,
                                       columnmapping=columnmap,
                                       ids_map=ids_map,resolution=1.0,
                                       filelabel="EDC",
                                       label="EDC")
    xl2.add_to_model()
    sampleobjects.append(xl2)
    outputobjects.append(xl2)
    xl2.set_psi_is_sampled(False)
    psi=xl2.get_psi(1.0)[0]
    psi.set_scale(0.05)


#####################################################
# Restraints setup
# Cross-link restraints
#####################################################
columnmap = {}
columnmap["Protein1"] = "Protein 1"
columnmap["Protein2"] = "Protein 2"
columnmap["Residue1"] = "Residue 1"
columnmap["Residue2"] = "Residue 2"
columnmap["IDScore"] = "p value"
columnmap["XLUniqueID"] = "XLUniqueID"

ids_map = IMP.pmi.tools.map()
ids_map.set_map_element(1.0, 1.0)

if (False):
    #----------------------------------------------------
    # wild type ScNup82 complex DSS XL data
    #----------------------------------------------------
    xl1 = IMP.pmi.restraints.crosslinking.ISDCrossLinkMS(simo,
                                                        '../data/XL_wtNup82_DSS_standardized_no_FG_2copies_Ambiguity3.csv',
                                                        length = 21.0,
                                                        slope = 0.00,
                                                        columnmapping = columnmap,
                                                        ids_map = ids_map,
                                                        resolution = 1.0,
                                                        inner_slope = 0.01,                                                    
                                                        filelabel = "scDSS",
                                                        label = "scDSS",
                                                        attributes_for_label = ["XLUniqueID"],
                                                        csvfile = True)
    xl1.add_to_model()
    sampleobjects.append(xl1)
    outputobjects.append(xl1)
    xl1.set_psi_is_sampled(False)
    psi1 = xl1.get_psi(1.0)[0]
    psi1.set_scale(0.05)

    #----------------------------------------------------
    # wild type SkNup82 complex DSS XL data
    #----------------------------------------------------
    xl2 = IMP.pmi.restraints.crosslinking.ISDCrossLinkMS(simo,
                                                        '../data/XL_skNup82_DSS_standardized_equiv_no_FG_2copies_Ambiguity3.csv',
                                                        length = 21.0,
                                                        slope = 0.00,
                                                        columnmapping = columnmap,
                                                        ids_map = ids_map,
                                                        resolution = 1.0,
                                                        inner_slope = 0.01,                                                    
                                                        filelabel = "skDSS",
                                                        label = "skDSS",
                                                        attributes_for_label = ["XLUniqueID"],
                                                        csvfile = True)
    xl2.add_to_model()
    sampleobjects.append(xl2)
    outputobjects.append(xl2)
    xl2.set_psi_is_sampled(False)
    psi2 = xl2.get_psi(1.0)[0]
    psi2.set_scale(0.05)
    
    #----------------------------------------------------
    # wild type ScNup82 complex EDC XL data
    #----------------------------------------------------
    xl3 = IMP.pmi.restraints.crosslinking.ISDCrossLinkMS(simo,
                                                        '../data/XL_wtNup82_EDC_standardized_no_FG_2copies_Ambiguity3.csv',
                                                        length = 16.0,
                                                        slope = 0.00,
                                                        columnmapping = columnmap,
                                                        ids_map = ids_map,
                                                        resolution = 1.0,
                                                        inner_slope = 0.01,                                                    
                                                        filelabel = "scEDC",
                                                        label = "scEDC",
                                                        attributes_for_label = ["XLUniqueID"],
                                                        csvfile = True)
    xl3.add_to_model()
    sampleobjects.append(xl3)
    outputobjects.append(xl3)
    xl3.set_psi_is_sampled(False)
    psi3 = xl3.get_psi(1.0)[0]
    psi3.set_scale(0.05)


if (False):
    """
    #----------------------------------------------------
    # XL Cliques (combined skDSS / wtDSS / wtEDC)
    #----------------------------------------------------
    xl_cliques_psi = 1.0
    xl_cliques_sigma = 5.0
    
    xl4 = IMP.pmi.restraints.crosslinking.ISDCrossLinkMS(simo,
                                                        '../data/XL_cliques_2copies.csv',
                                                        length = 10.0,
                                                        slope = 0.00,
                                                        columnmapping = columnmap,
                                                        ids_map = ids_map,
                                                        resolution = 1.0,
                                                        inner_slope = 0.01,                                                    
                                                        filelabel = "cliques",
                                                        label = "cliques",
                                                        attributes_for_label = ["XLUniqueID"],
                                                        csvfile = True)
    xl4.add_to_model()
    #xl4.set_weight(xl_cliques_weight)
    sampleobjects.append(xl4)
    outputobjects.append(xl4)
    xl4.set_psi_is_sampled(False)
    psi4 = xl4.get_psi(xl_cliques_psi)[0]
    psi4.set_scale(0.05)

    #xl4.set_sigma_is_sampled(False)
    #sigma4 = xl4.get_sigma(xl_cliques_sigma)[0]
    #sigma4.set_scale(1.0)
    """


#####################################################
# Restraints setup
# Distance restraints for homo-dimers
#####################################################
if (True):
    dist_min = 3.0
    dist_max = 30.0
    dr_weight = 100.0
    
    dr1 = IMP.pmi.restraints.basic.DistanceRestraint(simo,(1417,1417,"Nup159.1"), (1417,1417,"Nup159.2"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
    dr1.add_to_model()
    dr1.set_label("Nup159_1417-1417")
    dr1.set_weight(dr_weight)
    outputobjects.append(dr1)
    print(dr1.get_output())

    dr2 = IMP.pmi.restraints.basic.DistanceRestraint(simo,(1432,1432,"Nup159.1"), (1432,1432,"Nup159.2"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
    dr2.add_to_model()
    dr2.set_label("Nup159_1432-1432")
    dr2.set_weight(dr_weight)
    outputobjects.append(dr2)
    print(dr2.get_output())
        
    dr3 = IMP.pmi.restraints.basic.DistanceRestraint(simo,(1384,1384,"Nup159.1"), (1384,1384,"Nup159.2"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
    dr3.add_to_model()
    dr3.set_label("Nup159_1384-1384")
    dr3.set_weight(dr_weight)
    outputobjects.append(dr3)
    print(dr3.get_output())

    dr4 = IMP.pmi.restraints.basic.DistanceRestraint(simo,(1414,1414,"Nup159.1"), (1414,1414,"Nup159.2"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
    dr4.add_to_model()
    dr4.set_label("Nup159_1414-1414")
    dr4.set_weight(dr_weight)
    outputobjects.append(dr4)
    print(dr4.get_output())

    dr5 = IMP.pmi.restraints.basic.DistanceRestraint(simo,(1387,1387,"Nup159.1"), (1387,1387,"Nup159.2"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
    dr5.add_to_model()
    dr5.set_label("Nup159_1387-1387")
    dr5.set_weight(dr_weight)
    outputobjects.append(dr5)
    print(dr5.get_output())
    
    dr6 = IMP.pmi.restraints.basic.DistanceRestraint(simo,(517,517,"Nup82.1"), (517,517,"Nup82.2"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
    dr6.add_to_model()
    dr6.set_label("Nup82_517-517")
    dr6.set_weight(dr_weight)
    outputobjects.append(dr6)
    print(dr6.get_output())
    
    # by Ed Hurt
    if (False):
        """
        dr21 = IMP.pmi.restraints.basic.DistanceRestraint(simo,(1343,1343,"Nup159.1"), (1343,1343,"Nup159.2"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
        dr21.add_to_model()
        dr21.set_label("Nup159_1343-1343")
        dr21.set_weight(dr_weight)
        outputobjects.append(dr21)
        print(dr21.get_output())
        
        dr22 = IMP.pmi.restraints.basic.DistanceRestraint(simo,(541,541,"Nup82.1"), (541,541,"Nup82.2"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
        dr22.add_to_model()
        dr22.set_label("Nup82_541-541")
        dr22.set_weight(dr_weight)
        outputobjects.append(dr22)
        print(dr22.get_output())
        
        d23 = IMP.pmi.restraints.basic.DistanceRestraint(simo,(580,580,"Nup82.1"), (580,580,"Nup82.2"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
        d23.add_to_model()
        d23.set_label("Nup82_580-580")
        d23.set_weight(dr_weight)
        outputobjects.append(d23)
        print(d23.get_output())
        """
    
    print "\nDistance Restraints applied for homo-dimers !!"
    print "weight = ", dr_weight, "dist_min = ", dist_min, "dist_max = ", dist_max, "\n"

    # for skNup82 homo-dimers - DISABLED
    if (False):
        """
        dist_min = 3.0
        dist_max = 35.0
        dr_weight = 10.0
        
        dr7 = IMP.pmi.restraints.basic.DistanceRestraint(simo,(649,649,"Nup82.1"), (692,692,"Nup82.2"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
        dr7.add_to_model()
        dr7.set_label("Nup82_649-692a")
        dr7.set_weight(dr_weight)
        outputobjects.append(dr7)
        print(dr7.get_output())

        dr8 = IMP.pmi.restraints.basic.DistanceRestraint(simo,(649,649,"Nup82.2"), (692,692,"Nup82.1"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
        dr8.add_to_model()
        dr8.set_label("Nup82_649-692b")
        dr8.set_weight(dr_weight)
        outputobjects.append(dr8)
        print(dr8.get_output())

        dr9 = IMP.pmi.restraints.basic.DistanceRestraint(simo,(541,541,"Nup82.1"), (569,569,"Nup82.2"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
        dr9.add_to_model()
        dr9.set_label("Nup82_541-569a")
        dr9.set_weight(dr_weight)
        outputobjects.append(dr9)
        print(dr9.get_output())

        dr10 = IMP.pmi.restraints.basic.DistanceRestraint(simo,(541,541,"Nup82.2"), (569,569,"Nup82.1"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
        dr10.add_to_model()
        dr10.set_label("Nup82_541-569b")
        dr10.set_weight(dr_weight)
        outputobjects.append(dr10)
        print(dr10.get_output())

        print "\nDistance Restraints applied for skNup82 homo-dimers !!"
        print "weight = ", dr_weight, "dist_min = ", dist_min, "dist_max = ", dist_max, "\n"
        """



#####################################################
# Restraints setup
# Distance restraints for XL cliques
#####################################################
if (False):
    protein1 = columnmap["Protein1"]
    protein2 = columnmap["Protein2"]
    residue1 = columnmap["Residue1"]
    residue2 = columnmap["Residue2"]
    idscore = columnmap["IDScore"]
    xluniqueid = columnmap["XLUniqueID"]

    db = IMP.pmi.tools.get_db_from_csv('../data/XL_cliques_2copies.csv')

    dist_min = 3.0
    dist_max = 35.0
    dr_weight = 10.0
    
    for nxl, entry in enumerate(db):
        #print nxl, entry

        mol1 = entry[protein1]
        res1 = int(entry[residue1])
        mol2 = entry[protein2]
        res2 = int(entry[residue2])
        id_score = float(entry[idscore])
        xlunique_id = int(entry[xluniqueid])

        dr = IMP.pmi.restraints.basic.DistanceRestraint(simo,(res1,res1,mol1), (res2,res2,mol2), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
        dr.add_to_model()
        temp_label = mol1 + "_" + str(res1) + "-" + mol2 + "_" + str(res2)
        dr.set_label(temp_label)
        dr.set_weight(dr_weight)
        outputobjects.append(dr)
        print(dr.get_output())

    print "\nDistance Restraints applied for XL cliques !!"
    print "weight = ", dr_weight, "dist_min = ", dist_min, "dist_max = ", dist_max, "\n"



#####################################################
# Restraints setup
# EM 3D restraint using GMM
#####################################################
if (False):
    """
    # tail module em density
    mass = sum((IMP.atom.Mass(p).get_mass() for h in resdensities for p in IMP.atom.get_leaves(h)))
    gem = IMP.pmi.restraints.em.GaussianEMRestraint(resdensities,
                                                    '../data/em_density/spidervol_overlap.mrc.gmm.100.txt',
                                                    #'../data/em_density/emanvol.mrc.gmm.100.txt',
                                                    target_mass_scale=mass,
                                                    slope=0.000001,
                                                    target_radii_scale=3.0)
    gem.add_to_model()
    gem.set_weight(100.0)
    #gem.set_weight(200.0)
    #gem.center_model_on_target_density(simo)
    outputobjects.append(gem)
    
    
    # tail module em density
    mass2 = sum((IMP.atom.Mass(p).get_mass() for h in resdensities for p in IMP.atom.get_leaves(h)))
    gem2 = IMP.pmi.restraints.em.GaussianEMRestraint(resdensities,
                                                    #'../data/em_density/spidervol_overlap.mrc.gmm.100.txt',
                                                    '../data/em_density/emanvol.mrc.gmm.100.txt',
                                                    target_mass_scale=mass2,
                                                    slope=0.000001,
                                                    target_radii_scale=3.0)
    gem2.add_to_model()
    gem2.set_weight(100.0)
    #gem2.set_weight(200.0)
    #gem2.center_model_on_target_density(simo)
    outputobjects.append(gem2)
    """


#####################################################
# Metropolis Monte Carlo sampling with Replica Exchange
#####################################################
sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))
print "\nEVAL 1 : ", sf.evaluate(False), " (initial) - ", rank

if (True):
    simo.optimize_floppy_bodies(150)
    print "\nEVAL 2 : ", sf.evaluate(False), " (after calling optimize_floppy_bodies(150)) - ", rank

    initial_nframes = 200
    mc1 = IMP.pmi.macros.ReplicaExchange0(m,
                                        simo,
                                        monte_carlo_sample_objects = sampleobjects,
                                        output_objects = outputobjects,
                                        crosslink_restraints = [xl1, xl2],
                                        #crosslink_restraints = [xl4],
                                        monte_carlo_temperature = 1.0,
                                        replica_exchange_minimum_temperature = 1.0,
                                        replica_exchange_maximum_temperature = 2.5,
                                        number_of_best_scoring_models = 10,
                                        monte_carlo_steps = 10,
                                        number_of_frames = initial_nframes,
                                        write_initial_rmf = True,
                                        initial_rmf_name_suffix = "initial",
                                        stat_file_name_suffix = "stat",
                                        best_pdb_name_suffix = "model",
                                        do_clean_first = True,
                                        do_create_directories = True,
                                        global_output_directory = "pre-EM2D_output",
                                        rmf_dir = "rmfs/",
                                        best_pdb_dir = "pdbs/",
                                        replica_stat_file_suffix = "stat_replica")
    mc1.execute_macro()
    rex1 = mc1.get_replica_exchange_object()
    print "\nEVAL 3 : ", sf.evaluate(False), " (after performing the pre-sampling) - ", rank
else:
    rex1 = None
    print "\n>> NO pre-sampling"


#####################################################
# Restraints setup
# XL Restraints between Nup82 - Nup84
#####################################################
if (True):
    #----------------------------------------------------
    # Nup82 - Nup84 DSS XL data
    #----------------------------------------------------
    xl5 = IMP.pmi.restraints.crosslinking.ISDCrossLinkMS(simo,
                                                        '../data_nup84/XL_n82n84_DSS.csv',
                                                        length = 21.0,
                                                        slope = 0.00,
                                                        columnmapping = columnmap,
                                                        ids_map = ids_map,
                                                        resolution = 1.0,
                                                        inner_slope = 0.05,                                                    
                                                        filelabel = "wtDSS",
                                                        label = "wtDSS",
                                                        attributes_for_label = ["XLUniqueID"],
                                                        csvfile = True)
    xl5.add_to_model()
    sampleobjects.append(xl5)
    outputobjects.append(xl5)
    xl5.set_psi_is_sampled(False)
    psi2 = xl5.get_psi(1.0)[0]
    psi2.set_scale(0.05)


#####################################################
# Restraints setup+
# EM 2D restraint for each class
#####################################################
if (inputs.em2d_input is not None):
    #images = [inputs.em2d_input]
    images = []
    for class_num in range(0, 11):  #for 0-10
        pgm = "../data/em2d/" + str(class_num) + ".pgm"
        images.append(pgm)
    for class_num in range(12, 19): #for 12-18
        pgm = "../data/em2d/" + str(class_num) + ".pgm"
        images.append(pgm)
    for class_num in range(20, 23): #for 20-22
        pgm = "../data/em2d/" + str(class_num) + ".pgm"
        images.append(pgm)
    print images
    
    em2d = em2d_nup82.ElectronMicroscopy2D (simo,
                                            images,
                                            pixel_size = 3.23,
                                            image_resolution = 35.0,
                                            projection_number = 100,
                                            #projection_number = 400,
                                            resolution = 1.0,
                                            n_components = 1)
    em2d.add_to_model()
    em2d.set_weight(em2d_weight)
    outputobjects.append(em2d)

    sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))
    print "\nEVAL 4 : ", sf.evaluate(False), " (after applying the EM 2D restraint) - ", rank


#####################################################
# Metropolis Monte Carlo sampling with Replica Exchange
#####################################################
mc2 = IMP.pmi.macros.ReplicaExchange0(m,
                                    simo,
                                    monte_carlo_sample_objects = sampleobjects,
                                    output_objects = outputobjects,
                                    crosslink_restraints = [xl1, xl2, xl5],
                                    #crosslink_restraints = [xl5],
                                    monte_carlo_temperature = 1.0,
                                    replica_exchange_minimum_temperature = 1.0,
                                    replica_exchange_maximum_temperature = 2.5,
                                    number_of_best_scoring_models = 100,
                                    #number_of_best_scoring_models = int(inputs.nrepeats)-200,
                                    monte_carlo_steps = 10,
                                    number_of_frames = int(inputs.nrepeats),
                                    write_initial_rmf = True,
                                    initial_rmf_name_suffix = "initial",
                                    stat_file_name_suffix = "stat",
                                    best_pdb_name_suffix = "model",
                                    do_clean_first = True,
                                    do_create_directories = True,
                                    global_output_directory = inputs.folder_output,
                                    rmf_dir = "rmfs/",
                                    best_pdb_dir = "pdbs/",
                                    replica_stat_file_suffix = "stat_replica",
                                    replica_exchange_object = rex1)
mc2.execute_macro()
print "\nEVAL 5 : ", sf.evaluate(False), " (final evaluation) - ", rank

