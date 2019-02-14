#!/usr/bin/env python
#####################################################
# Last Update: April 19th, 2015
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
import IMP.pmi1.restraints.crosslinking
import IMP.pmi1.restraints.stereochemistry
#import IMP.pmi1.restraints.em
import em2d_nup82
#import IMP.pmi1.restraints.em2d
import IMP.pmi1.restraints.basic
import IMP.pmi1.restraints.proteomics
#import representation_nup82
import IMP.pmi1.representation
import IMP.pmi1.macros
import IMP.pmi1.restraints
import IMP.pmi1.tools
import IMP.pmi1.output
import IMP.pmi1.samplers
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
simo = IMP.pmi1.representation.Representation(m,upperharmonic=True,disorderedlength=False)


#####################################################
# setting up parameters
#####################################################
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
print "rank = ", rank

rbmaxtrans = 2.00
fbmaxtrans = 3.00
rbmaxrot = 0.04
outputobjects = []
sampleobjects = []

res_cry = int(inputs.res_cry)
res_ev = int(inputs.res_ev)
res_conn = int(inputs.res_cry)
beadsize = 20
beadsize50 = 50
beadsize100 = 100
em2d_weight = float(inputs.weight)

datadirectory = "../data/"
fasta_files = "../data/protein_fasta."


#####################################################
# REPRESENTATION
#####################################################
# compname  hier_name      color   fastafile                  fastaid    pdbname                       chain res_range       read_em_files bead_size  rigid_body super_rigid_body em_num_components em_txt_file_name em_mrc_file_name chain_of_super_rb
domains = \
[("Dyn2.1", "Dyn2.1",      0.48,   fasta_files+"Dyn2.txt",    "Dyn2",    datadirectory+"4DS1.pdb",     "A",  (1,92,0),       None,         beadsize,       19,     [100    ,19],   4,               None,            None, [3]),
 ("Dyn2.2", "Dyn2.2",      0.65,   fasta_files+"Dyn2.txt",    "Dyn2",    datadirectory+"4DS1.pdb",     "C",  (1,92,0),       None,         beadsize,       19,     [100    ,19],   4,               None,            None, [3]),
 
 ("Nup82.1", "Nup82.1_1",   0.0,   fasta_files+"Nup82.txt",   "Nup82",   datadirectory+"3PBP.pdb",     "G",  (1,521,0),      None,         beadsize,       0,      [100,101,0],    18,              None,            None, [0]),
 ("Nup82.1", "Nup82.1_2",   0.0,   fasta_files+"Nup82.txt",   "Nup82",   "IDEAL_HELIX",                " ",  (522,571,0),    None,         beadsize,       1,      [100,101,1],    4,               None,            None, [0]),
 ("Nup82.1", "Nup82.1_3",   0.0,   fasta_files+"Nup82.txt",   "Nup82",   "BEADS",                      " ",  (572,588,0),    None,         beadsize,       2,      [100,101,2],    0,               None,            None, [0]),
 ("Nup82.1", "Nup82.1_4",   0.0,   fasta_files+"Nup82.txt",   "Nup82",   "IDEAL_HELIX",                " ",  (589,613,0),    None,         beadsize,       3,      [100,101,3],    2,               None,            None, [0]),
 ("Nup82.1", "Nup82.1_5",   0.0,   fasta_files+"Nup82.txt",   "Nup82",   "BEADS",                      " ",  (614,624,0),    None,         beadsize,       4,      [100,101,4],    0,               None,            None, [0]),
 ("Nup82.1", "Nup82.1_6_1", 0.0,   fasta_files+"Nup82.txt",   "Nup82",   "IDEAL_HELIX",                " ",  (625,641,0),    None,         beadsize,       5,      [100,101,5],    1,               None,            None, [0]),
 ("Nup82.1", "Nup82.1_6_2", 0.0,   fasta_files+"Nup82.txt",   "Nup82",   "BEADS",                      " ",  (642,647,0),    None,         beadsize,       51,     [100,101,51],   0,               None,            None, [0]),
 ("Nup82.1", "Nup82.1_6_3", 0.0,   fasta_files+"Nup82.txt",   "Nup82",   "IDEAL_HELIX",                " ",  (648,671,0),    None,         beadsize,       52,     [100,101,52],   2,               None,            None, [0]),
 ("Nup82.1", "Nup82.1_7",   0.0,   fasta_files+"Nup82.txt",   "Nup82",   "BEADS",                      " ",  (672,676,0),    None,         beadsize,       6,      [100,101,6],    0,               None,            None, [0]),
 ("Nup82.1", "Nup82.1_8",   0.0,   fasta_files+"Nup82.txt",   "Nup82",   "IDEAL_HELIX",                " ",  (677,709,0),    None,         beadsize,       7,      [100,101,7],    2,               None,            None, [0]),
 ("Nup82.1", "Nup82.1_9",   0.0,   fasta_files+"Nup82.txt",   "Nup82",   "BEADS",                      " ",  (710,713,0),    None,         beadsize,       8,      [100,101,8],    0,               None,            None, [0]),
 ("Nup82.2", "Nup82.2_1",   0.15,  fasta_files+"Nup82.txt",   "Nup82",   datadirectory+"3PBP.pdb",     "G",  (1,521,0),      None,         beadsize,       40,     [100,102,40],   18,              None,            None, [4]),
 ("Nup82.2", "Nup82.2_2",   0.15,  fasta_files+"Nup82.txt",   "Nup82",   "IDEAL_HELIX",                " ",  (522,571,0),    None,         beadsize,       41,     [100,102,41],   4,               None,            None, [4]),
 ("Nup82.2", "Nup82.2_3",   0.15,  fasta_files+"Nup82.txt",   "Nup82",   "BEADS",                      " ",  (572,588,0),    None,         beadsize,       42,     [100,102,42],   0,               None,            None, [4]),
 ("Nup82.2", "Nup82.2_4",   0.15,  fasta_files+"Nup82.txt",   "Nup82",   "IDEAL_HELIX",                " ",  (589,613,0),    None,         beadsize,       43,     [100,102,43],   2,               None,            None, [4]),
 ("Nup82.2", "Nup82.2_5",   0.15,  fasta_files+"Nup82.txt",   "Nup82",   "BEADS",                      " ",  (614,624,0),    None,         beadsize,       44,     [100,102,44],   0,               None,            None, [4]),
 ("Nup82.2", "Nup82.2_6_1", 0.15,  fasta_files+"Nup82.txt",   "Nup82",   "IDEAL_HELIX",                " ",  (625,641,0),    None,         beadsize,       45,     [100,102,45],   1,               None,            None, [4]),
 ("Nup82.2", "Nup82.2_6_2", 0.15,  fasta_files+"Nup82.txt",   "Nup82",   "BEADS",                      " ",  (642,647,0),    None,         beadsize,       451,    [100,102,451],  0,               None,            None, [4]),
 ("Nup82.2", "Nup82.2_6_3", 0.15,  fasta_files+"Nup82.txt",   "Nup82",   "IDEAL_HELIX",                " ",  (648,671,0),    None,         beadsize,       452,    [100,102,452],  2,               None,            None, [4]),
 ("Nup82.2", "Nup82.2_7",   0.15,  fasta_files+"Nup82.txt",   "Nup82",   "BEADS",                      " ",  (672,676,0),    None,         beadsize,       46,     [100,102,46],   0,               None,            None, [4]),
 ("Nup82.2", "Nup82.2_8",   0.15,  fasta_files+"Nup82.txt",   "Nup82",   "IDEAL_HELIX",                " ",  (677,709,0),    None,         beadsize,       47,     [100,102,47],   2,               None,            None, [4]),
 ("Nup82.2", "Nup82.2_9",   0.15,  fasta_files+"Nup82.txt",   "Nup82",   "BEADS",                      " ",  (710,713,0),    None,         beadsize,       48,     [100,102,48],   0,               None,            None, [4]),

 ("Nup159.1", "Nup159.1_1", 1.0,   fasta_files+"Nup159.txt",  "Nup159",  datadirectory+"1XIP.pdb",     "A",  (1,381,0),      None,         beadsize,       17,     [100,104,17],   0,               None,            None, [2]),
 ("Nup159.1", "Nup159.1_2", 1.0,   fasta_files+"Nup159.txt",  "Nup159",  "BEADS",                      " ",  (382,1116,0),   None,         beadsize100,    18,     [100,104,18],   0,               None,            None, [2]),
 ("Nup159.1", "Nup159.1_3", 1.0,   fasta_files+"Nup159.txt",  "Nup159",  datadirectory+"4DS1.pdb",     "B",  (1117,1126,0),  None,         beadsize,       19,     [100,104,19],   1,               None,            None, [2]),
 ("Nup159.1", "Nup159.1_4", 1.0,   fasta_files+"Nup159.txt",  "Nup159",  "BEADS",                      " ",  (1127,1208,0),  None,         beadsize,       20,     [100,104,20],   0,               None,            None, [2]),
 ("Nup159.1", "Nup159.1_5", 1.0,   fasta_files+"Nup159.txt",  "Nup159",  "IDEAL_HELIX",                " ",  (1209,1241,0),  None,         beadsize,       21,     [100,104,21],   2,               None,            None, [2]),
 ("Nup159.1", "Nup159.1_6", 1.0,   fasta_files+"Nup159.txt",  "Nup159",  "BEADS",                      " ",  (1242,1265,0),  None,         beadsize,       22,     [100,104,22],   0,               None,            None, [2]),
 ("Nup159.1", "Nup159.1_7", 1.0,   fasta_files+"Nup159.txt",  "Nup159",  "IDEAL_HELIX",                " ",  (1266,1320,0),  None,         beadsize,       23,     [100,104,23],   2,               None,            None, [2]),
 ("Nup159.1", "Nup159.1_8", 1.0,   fasta_files+"Nup159.txt",  "Nup159",  "BEADS",                      " ",  (1321,1334,0),  None,         beadsize,       24,     [100,104,24],   0,               None,            None, [2]),
 ("Nup159.1", "Nup159.1_9", 1.0,   fasta_files+"Nup159.txt",  "Nup159",  "IDEAL_HELIX",                " ",  (1335,1369,0),  None,         beadsize,       25,     [100,104,25],   2,               None,            None, [2]),
 ("Nup159.1", "Nup159.1_10",1.0,   fasta_files+"Nup159.txt",  "Nup159",  "BEADS",                      " ",  (1370,1379,0),  None,         beadsize,       26,     [100,104,26],   0,               None,            None, [2]),
 ("Nup159.1", "Nup159.1_11",1.0,   fasta_files+"Nup159.txt",  "Nup159",  "IDEAL_HELIX",                " ",  (1380,1410,0),  None,         beadsize,       27,     [100,104,27],   2,               None,            None, [2]),
 ("Nup159.1", "Nup159.1_12",1.0,   fasta_files+"Nup159.txt",  "Nup159",  "BEADS",                      " ",  (1411,1428,0),  None,         beadsize,       28,     [100,104,28],   0,               None,            None, [2]),
 ("Nup159.1", "Nup159.1_13",1.0,   fasta_files+"Nup159.txt",  "Nup159",  datadirectory+"3PBP.pdb",     "I",  (1429,1460,0),  None,         beadsize,       0,      [100,104,0],    2,               None,            None, [2]),
 ("Nup159.2", "Nup159.2_1", 0.9,   fasta_files+"Nup159.txt",  "Nup159",  datadirectory+"1XIP.pdb",     "A",  (1,381,0),      None,         beadsize,       57,     [100,108,57],   0,               None,            None, [8]),
 ("Nup159.2", "Nup159.2_2", 0.9,   fasta_files+"Nup159.txt",  "Nup159",  "BEADS",                      " ",  (382,1116,0),   None,         beadsize100,    58,     [100,108,58],   0,               None,            None, [8]),
 ("Nup159.2", "Nup159.2_3", 0.9,   fasta_files+"Nup159.txt",  "Nup159",  datadirectory+"4DS1.pdb",     "D",  (1117,1126,0),  None,         beadsize,       19,     [100,108,19],   1,               None,            None, [8]),
 ("Nup159.2", "Nup159.2_4", 0.9,   fasta_files+"Nup159.txt",  "Nup159",  "BEADS",                      " ",  (1127,1208,0),  None,         beadsize,       60,     [100,108,60],   0,               None,            None, [8]),
 ("Nup159.2", "Nup159.2_5", 0.9,   fasta_files+"Nup159.txt",  "Nup159",  "IDEAL_HELIX",                " ",  (1209,1241,0),  None,         beadsize,       61,     [100,108,61],   2,               None,            None, [8]),
 ("Nup159.2", "Nup159.2_6", 0.9,   fasta_files+"Nup159.txt",  "Nup159",  "BEADS",                      " ",  (1242,1265,0),  None,         beadsize,       62,     [100,108,62],   0,               None,            None, [8]),
 ("Nup159.2", "Nup159.2_7", 0.9,   fasta_files+"Nup159.txt",  "Nup159",  "IDEAL_HELIX",                " ",  (1266,1320,0),  None,         beadsize,       63,     [100,108,63],   2,               None,            None, [8]),
 ("Nup159.2", "Nup159.2_8", 0.9,   fasta_files+"Nup159.txt",  "Nup159",  "BEADS",                      " ",  (1321,1334,0),  None,         beadsize,       64,     [100,108,64],   0,               None,            None, [8]),
 ("Nup159.2", "Nup159.2_9", 0.9,   fasta_files+"Nup159.txt",  "Nup159",  "IDEAL_HELIX",                " ",  (1335,1369,0),  None,         beadsize,       65,     [100,108,65],   2,               None,            None, [8]),
 ("Nup159.2", "Nup159.2_10",0.9,   fasta_files+"Nup159.txt",  "Nup159",  "BEADS",                      " ",  (1370,1379,0),  None,         beadsize,       66,     [100,108,66],   0,               None,            None, [8]),
 ("Nup159.2", "Nup159.2_11",0.9,   fasta_files+"Nup159.txt",  "Nup159",  "IDEAL_HELIX",                " ",  (1380,1410,0),  None,         beadsize,       67,     [100,108,67],   2,               None,            None, [8]),
 ("Nup159.2", "Nup159.2_12",0.9,   fasta_files+"Nup159.txt",  "Nup159",  "BEADS",                      " ",  (1411,1428,0),  None,         beadsize,       68,     [100,108,68],   0,               None,            None, [8]),
 ("Nup159.2", "Nup159.2_13",0.9,   fasta_files+"Nup159.txt",  "Nup159",  datadirectory+"3PBP.pdb",     "I",  (1429,1460,0),  None,         beadsize,       40,     [100,108,40],   2,               None,            None, [8]),
 
 ("Nsp1.1",  "Nsp1.1_1",   0.3,    fasta_files+"Nsp1.txt",    "Nsp1",    "BEADS",                      " ",  (1,649,0),      None,         beadsize100,    9,      [100,103,9],    0,               None,            None, [1]),
 ("Nsp1.1",  "Nsp1.1_2",   0.3,    fasta_files+"Nsp1.txt",    "Nsp1",    "IDEAL_HELIX",                " ",  (650,661,0),    None,         beadsize,       10,     [100,103,10],   1,               None,            None, [1]),
 ("Nsp1.1",  "Nsp1.1_3",   0.3,    fasta_files+"Nsp1.txt",    "Nsp1",    "BEADS",                      " ",  (662,662,0),    None,         beadsize,       11,     [100,103,11],   0,               None,            None, [1]),
 ("Nsp1.1",  "Nsp1.1_4",   0.3,    fasta_files+"Nsp1.txt",    "Nsp1",    "IDEAL_HELIX",                " ",  (663,727,0),    None,         beadsize,       12,     [100,103,12],   4,               None,            None, [1]),
 ("Nsp1.1",  "Nsp1.1_5",   0.3,    fasta_files+"Nsp1.txt",    "Nsp1",    "BEADS",                      " ",  (728,737,0),    None,         beadsize,       13,     [100,103,13],   0,               None,            None, [1]),
 ("Nsp1.1",  "Nsp1.1_6",   0.3,    fasta_files+"Nsp1.txt",    "Nsp1",    "IDEAL_HELIX",                " ",  (738,771,0),    None,         beadsize,       14,     [100,103,14],   2,               None,            None, [1]),
 ("Nsp1.1",  "Nsp1.1_7",   0.3,    fasta_files+"Nsp1.txt",    "Nsp1",    "BEADS",                      " ",  (772,788,0),    None,         beadsize,       15,     [100,103,15],   0,               None,            None, [1]),
 ("Nsp1.1",  "Nsp1.1_8",   0.3,    fasta_files+"Nsp1.txt",    "Nsp1",    "IDEAL_HELIX",                " ",  (789,823,0),    None,         beadsize,       16,     [100,103,16],   2,               None,            None, [1]),
 ("Nsp1.2",  "Nsp1.2_1",   0.38,   fasta_files+"Nsp1.txt",    "Nsp1",    "BEADS",                      " ",  (1,649,0),      None,         beadsize100,    79,     [100,107,79],   0,               None,            None, [7]),
 ("Nsp1.2",  "Nsp1.2_2",   0.38,   fasta_files+"Nsp1.txt",    "Nsp1",    "IDEAL_HELIX",                " ",  (650,661,0),    None,         beadsize,       80,     [100,107,80],   1,               None,            None, [7]),
 ("Nsp1.2",  "Nsp1.2_3",   0.38,   fasta_files+"Nsp1.txt",    "Nsp1",    "BEADS",                      " ",  (662,662,0),    None,         beadsize,       81,     [100,107,81],   0,               None,            None, [7]),
 ("Nsp1.2",  "Nsp1.2_4",   0.38,   fasta_files+"Nsp1.txt",    "Nsp1",    "IDEAL_HELIX",                " ",  (663,727,0),    None,         beadsize,       82,     [100,107,82],   4,               None,            None, [7]),
 ("Nsp1.2",  "Nsp1.2_5",   0.38,   fasta_files+"Nsp1.txt",    "Nsp1",    "BEADS",                      " ",  (728,737,0),    None,         beadsize,       83,     [100,107,83],   0,               None,            None, [7]),
 ("Nsp1.2",  "Nsp1.2_6",   0.38,   fasta_files+"Nsp1.txt",    "Nsp1",    "IDEAL_HELIX",                " ",  (738,771,0),    None,         beadsize,       84,     [100,107,84],   2,               None,            None, [7]),
 ("Nsp1.2",  "Nsp1.2_7",   0.38,   fasta_files+"Nsp1.txt",    "Nsp1",    "BEADS",                      " ",  (772,788,0),    None,         beadsize,       85,     [100,107,85],   0,               None,            None, [7]),
 ("Nsp1.2",  "Nsp1.2_8",   0.38,   fasta_files+"Nsp1.txt",    "Nsp1",    "IDEAL_HELIX",                " ",  (789,823,0),    None,         beadsize,       86,     [100,107,86],   2,               None,            None, [7]),
 
 ("Nup116.1", "Nup116.1_1", 0.75,  fasta_files+"Nup116.txt",  "Nup116",  "BEADS",                      " ",  (1,750,0),      None,         beadsize100,    91,     [100,105,91],   0,               None,            None, [5]),
 ("Nup116.1", "Nup116.1_2", 0.75,  fasta_files+"Nup116.txt",  "Nup116",  datadirectory+"3PBP.pdb",     "H",  (751,1113,0),   None,         beadsize,       0,      [100,105,0],    0,               None,            None, [5]),
 ("Nup116.2", "Nup116.2_1", 0.8,   fasta_files+"Nup116.txt",  "Nup116",  "BEADS",                      " ",  (1,750,0),      None,         beadsize100,    92,     [100,106,92],   0,               None,            None, [6]),
 ("Nup116.2", "Nup116.2_2", 0.8,   fasta_files+"Nup116.txt",  "Nup116",  datadirectory+"3PBP.pdb",     "H",  (751,1113,0),   None,         beadsize,       40,     [100,106,40],   0,               None,            None, [6])]

bm1 = IMP.pmi1.macros.BuildModel1(simo)
#bm1.set_gmm_models_directory(datadirectory + "em_gmm_model/")

if (inputs.rmf_input is not None) :
    n82=set([s[0] for s in domains])
    for d in list(n82):
        bm1.set_rmf_file(d, inputs.rmf_input, int(inputs.rmf_frame_number))

bm1.build_model(data_structure = domains)
#bm1.scale_bead_radii(40,0.8)
#resdensities = bm1.get_density_hierarchies([t[1] for t in domains])
#print resdensities; exit()

#model_ps = []
#for h in self.densities:
#    model_ps += IMP.atom.get_leaves(h)


#####################################################
# randomize the initial configuration
#####################################################
if (inputs.rmf_input is None) :
    simo.shuffle_configuration(100)


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
ev = IMP.pmi1.restraints.stereochemistry.ExcludedVolumeSphere(simo, resolution = res_ev)
ev.add_to_model()
outputobjects.append(ev)
print(ev.get_output())
print "ExcludedVolumeSphere !!\n"


#####################################################
# Restraints setup
# External Barrier restraint
#####################################################
eb = IMP.pmi1.restraints.basic.ExternalBarrier(simo, radius = 300)
eb.add_to_model()
outputobjects.append(eb)
print(eb.get_output())
print "ExternalBarrier !!\n"


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

ids_map = IMP.pmi1.tools.map()
ids_map.set_map_element(1.0, 1.0)

if (True):
    #----------------------------------------------------
    # wild type ScNup82 complex DSS XL data
    #----------------------------------------------------
    xl1 = IMP.pmi1.restraints.crosslinking.ISDCrossLinkMS(simo,
                                                        '../data/XL_wtNup82_DSS_standardized_2copies.csv',
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
    xl2 = IMP.pmi1.restraints.crosslinking.ISDCrossLinkMS(simo,
                                                        '../data/XL_skNup82_DSS_standardized_equiv_2copies.csv',
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
    xl3 = IMP.pmi1.restraints.crosslinking.ISDCrossLinkMS(simo,
                                                        '../data/XL_wtNup82_EDC_standardized_2copies.csv',
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


#####################################################
# Restraints setup
# Distance restraints for homo-dimers
#####################################################
if (True):
    dist_min = 3.0
    dist_max = 25.0
    dr_weight = 100.0
    
    dr1 = IMP.pmi1.restraints.basic.DistanceRestraint(simo,(1417,1417,"Nup159.1"), (1417,1417,"Nup159.2"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
    dr1.add_to_model()
    dr1.set_label("Nup159_1417-1417")
    dr1.set_weight(dr_weight)
    outputobjects.append(dr1)
    print(dr1.get_output())

    dr2 = IMP.pmi1.restraints.basic.DistanceRestraint(simo,(1432,1432,"Nup159.1"), (1432,1432,"Nup159.2"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
    dr2.add_to_model()
    dr2.set_label("Nup159_1432-1432")
    dr2.set_weight(dr_weight)
    outputobjects.append(dr2)
    print(dr2.get_output())
        
    dr3 = IMP.pmi1.restraints.basic.DistanceRestraint(simo,(1384,1384,"Nup159.1"), (1384,1384,"Nup159.2"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
    dr3.add_to_model()
    dr3.set_label("Nup159_1384-1384")
    dr3.set_weight(dr_weight)
    outputobjects.append(dr3)
    print(dr3.get_output())

    dr4 = IMP.pmi1.restraints.basic.DistanceRestraint(simo,(1414,1414,"Nup159.1"), (1414,1414,"Nup159.2"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
    dr4.add_to_model()
    dr4.set_label("Nup159_1414-1414")
    dr4.set_weight(dr_weight)
    outputobjects.append(dr4)
    print(dr4.get_output())

    dr5 = IMP.pmi1.restraints.basic.DistanceRestraint(simo,(1387,1387,"Nup159.1"), (1387,1387,"Nup159.2"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
    dr5.add_to_model()
    dr5.set_label("Nup159_1387-1387")
    dr5.set_weight(dr_weight)
    outputobjects.append(dr5)
    print(dr5.get_output())
    
    dr6 = IMP.pmi1.restraints.basic.DistanceRestraint(simo,(517,517,"Nup82.1"), (517,517,"Nup82.2"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
    dr6.add_to_model()
    dr6.set_label("Nup82_517-517")
    dr6.set_weight(dr_weight)
    outputobjects.append(dr6)
    print(dr6.get_output())

    print "\nDistance Restraints applied for homo-dimers !!"
    print "weight = ", dr_weight, "dist_min = ", dist_min, "dist_max = ", dist_max, "\n"


#####################################################
# Restraints setup
# EM 3D restraint using GMM
#####################################################
if (False):
    """
    # tail module em density
    mass = sum((IMP.atom.Mass(p).get_mass() for h in resdensities for p in IMP.atom.get_leaves(h)))
    gem = IMP.pmi1.restraints.em.GaussianEMRestraint(resdensities,
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
    gem2 = IMP.pmi1.restraints.em.GaussianEMRestraint(resdensities,
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
sf = IMP.core.RestraintsScoringFunction(IMP.pmi1.tools.get_restraint_set(m))
print "\nEVAL 1 : ", sf.evaluate(False), " (initial) - ", rank

if (False):
    simo.optimize_floppy_bodies(150)
    print "\nEVAL 2 : ", sf.evaluate(False), " (after calling optimize_floppy_bodies(150)) - ", rank

    initial_nframes = 5
    mc1 = IMP.pmi1.macros.ReplicaExchange0(m,
                                        simo,
                                        monte_carlo_sample_objects = sampleobjects,
                                        output_objects = outputobjects,
                                        crosslink_restraints = [xl1, xl2, xl3],
                                        #crosslink_restraints = [xl1],
                                        monte_carlo_temperature = 1.0,
                                        replica_exchange_minimum_temperature = 1.0,
                                        replica_exchange_maximum_temperature = 2.5,
                                        number_of_best_scoring_models = initial_nframes,
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
# EM 2D restraint for each class
#####################################################
if (inputs.em2d_input is not None):
    images = [inputs.em2d_input]
    
    em2d = em2d_nup82.ElectronMicroscopy2D (simo,
                                            images,
                                            resolution = 1.0,
                                            pixel_size = 3.23,
                                            image_resolution = 15.0,
                                            #projection_number = 100)
                                            projection_number = 400)
    em2d.add_to_model()
    em2d.set_weight(em2d_weight)
    outputobjects.append(em2d)

    sf = IMP.core.RestraintsScoringFunction(IMP.pmi1.tools.get_restraint_set(m))
    print "\nEVAL 4 : ", sf.evaluate(False), " (after applying the EM 2D restraint) - ", rank


#####################################################
# Metropolis Monte Carlo sampling with Replica Exchange
#####################################################
mc2 = IMP.pmi1.macros.ReplicaExchange0(m,
                                    simo,
                                    monte_carlo_sample_objects = sampleobjects,
                                    output_objects = outputobjects,
                                    crosslink_restraints = [xl1, xl2, xl3],
                                    #crosslink_restraints = [xl1],
                                    monte_carlo_temperature = 1.0,
                                    replica_exchange_minimum_temperature = 1.0,
                                    replica_exchange_maximum_temperature = 2.5,
                                    number_of_best_scoring_models = 500,
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

