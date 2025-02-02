#!/usr/bin/env python
#####################################################
# Last Update: September 10th, 2015
# by Seung Joong Kim and Riccardo Pellarin
# at Andrej Sali group, University of California San Francisco (UCSF)
#####################################################
from __future__ import print_function
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
import IMP.pmi1.mmcif
import ihm
import ihm.analysis
import ihm.dataset
import ihm.cross_linkers
try:
    import ihm.reference
except ImportError:
    pass
try:
    import ihm.citations
except ImportError:
    pass
#import representation_nup82
import IMP.pmi1.representation
import IMP.pmi1.macros
import IMP.pmi1.restraints
import IMP.pmi1.tools
import IMP.pmi1.output
import IMP.pmi1.samplers
import saxs
import random

import os
import sys

sys.path.append('../util/')
import make_archive


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
parser.add_argument('--dry-run', action='store_true',
                    help="Dry run (do not do any sampling)")
parser.add_argument('--mmcif', action='store', type=str, default=None,
                    help="Record modeling protocol in a named mmCIF file")

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
print(inputs)


#####################################################
# setting up topology
#####################################################
m = IMP.Model()
simo = IMP.pmi1.representation.Representation(m,upperharmonic=True,disorderedlength=False)
simo.dry_run = inputs.dry_run
#simo = representation_nup82.Representation(m,upperharmonic=True,disorderedlength=False)


#####################################################
# setting up parameters
#####################################################
try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
except ImportError:
    rank = 0

print("rank = ", rank)

if inputs.mmcif:
    # Record the modeling protocol to an mmCIF file
    po = IMP.pmi1.mmcif.ProtocolOutput(open(inputs.mmcif, 'w'))
    simo.add_protocol_output(po)
    po.system.title = ('Structure of the S. cerevisiae nuclear pore complex '
                       'cytoplasmic mRNA export platform, Nup82')
    # Add publication
    po.system.citations.append(ihm.Citation.from_pubmed_id(27839866))

    # Point to repositories where files are deposited
    zenodo_id = '1256259'
    doi = '10.5281/zenodo.' + zenodo_id
    url_top = 'https://zenodo.org/record/%s/files' % zenodo_id
    for subdir, zipname in make_archive.ARCHIVES.items():
        simo.add_metadata(ihm.location.Repository(
              doi=doi, root="../%s" % subdir,
              url="%s/%s.zip" % (url_top, zipname),
              top_directory=os.path.basename(subdir)))
        simo.add_metadata(ihm.location.Repository(
              doi=doi, root="..", url="%s/nup82-master.zip" % url_top,
              top_directory="nup82-master"))

# We used HHpred to detect remote homologs for some input subunits
s = ihm.Software(
          name='HHpred', classification='protein homology detection',
          description='Protein homology detection by HMM-HMM comparison',
          version='2.0.16',
          location='https://toolkit.tuebingen.mpg.de/hhpred')
if hasattr(ihm, 'citations'):
    s.citation = ihm.citations.hhpred
simo.add_metadata(s)

# We used PSIPRED to predict secondary structure for subunits
s = ihm.Software(
          name='PSIPRED', classification='secondary structure prediction',
          description='Protein secondary structure prediction based on '
                      'position-specific scoring matrices',
          version='4.0',
          location='http://bioinf.cs.ucl.ac.uk/psipred/')
if hasattr(ihm, 'citations'):
    s.citation = ihm.citations.psipred
simo.add_metadata(s)

# We used DISOPRED to predict (and remove) disordered regions in the subunits
s = ihm.Software(
          name='DISOPRED', classification='disorder prediction',
          description='prediction of protein disorder', version=3,
          location='http://bioinf.cs.ucl.ac.uk/psipred/?disopred=1')
if hasattr(ihm, 'citations'):
    s.citation = ihm.citations.disopred
simo.add_metadata(s)

# We used DomPred to predict domains
s = ihm.Software(
          name='DomPred', classification='domain prediction',
          description='prediction of protein domains',
          location='http://bioinf.cs.ucl.ac.uk/psipred/?dompred=1')
if hasattr(ihm, 'citations'):
    s.citation = ihm.Citation(
        pmid='17430199',
        title='Computer-assisted protein domain boundary prediction '
              'using the DomPred server.',
        journal='Curr Protein Pept Sci', volume=8, page_range=('181', '188'),
        year=2007, authors=['Bryson K', 'Cozzetto D', 'Jones DT'],
        doi='10.2174/138920307780363415')
simo.add_metadata(s)

# We used COILS/PCOILS and Multicoil2 to predict coiled-coil regions
s = ihm.Software(
          name='COILS/PCOILS', classification='coiled-coil prediction',
          description='prediction of protein coiled-coil regions',
          location='https://toolkit.tuebingen.mpg.de/pcoils')
if hasattr(ihm, 'citations'):
    s.citation = ihm.Citation(
        pmid='2031185',
        title='Predicting coiled coils from protein sequences.',
        journal='Science', volume=252, page_range=('1162', '1164'),
        year=1991, authors=['Lupas A', 'Van Dyke M', 'Stock J'],
        doi='10.1126/science.252.5009.1162')
simo.add_metadata(s)

s = ihm.Software(
          name='Multicoil2', classification='coiled-coil prediction',
          description='prediction of protein coiled-coil regions',
          location='http://groups.csail.mit.edu/cb/multicoil2/cgi-bin/multicoil2.cgi')
if hasattr(ihm, 'citations'):
    s.citation = ihm.Citation(
        pmid='21901122',
        title='Multicoil2: predicting coiled coils and their oligomerization '
              'states from sequence in the twilight zone.',
        journal='PLoS One', volume=6, page_range='e23519', year=2011,
        authors=['Trigg J', 'Gutwin K', 'Keating AE', 'Berger B'],
        doi='10.1371/journal.pone.0023519')
simo.add_metadata(s)

rbmaxtrans = 2.00
fbmaxtrans = 3.00
rbmaxrot = 0.04
outputobjects = []
sampleobjects = []

res_cry = int(inputs.res_cry)
res_ev = int(inputs.res_ev)
res_conn = int(inputs.res_cry)
beadsize = 10
beadsize25 = 25
beadsize100 = 100
em2d_weight = float(inputs.weight)

datadirectory = "../data/"
fasta_files = "../data/protein_fasta."


#####################################################
# REPRESENTATION
#####################################################
# compname  hier_name      color   fastafile                  fastaid    pdbname                       chain res_range       read_em_files bead_size  rigid_body super_rigid_body em_num_components em_txt_file_name em_mrc_file_name chain_of_super_rb
domains = \
[("Dyn2.1",  "Dyn2.1",     0.48,   fasta_files+"Dyn2.txt",    "Dyn2",    datadirectory+"4DS1.pdb",     "A",  (  1, 92,0),    None,         beadsize,       0,      [0],        4,               None,            None, [3]),
 ("Dyn2.2",  "Dyn2.2",     0.65,   fasta_files+"Dyn2.txt",    "Dyn2",    datadirectory+"4DS1.pdb",     "C",  (  1, 92,0),    None,         beadsize,       0,      [0],        4,               None,            None, [3]),


 ("Nup82.1", "Nup82.1_1",   0.0,   fasta_files+"Nup82.txt",   "Nup82",   datadirectory+"3PBP1.pdb",    "G",  (  1,452,0),    None,         beadsize,       1,      [11,1],     18,              None,            None, [0]),
 ("Nup82.1", "Nup82.1_11",  0.0,   fasta_files+"Nup82.txt",   "Nup82",   "BEADS",                      " ",  (453,521,0),    None,         beadsize,       101,    [11,101],   0,               None,            None, [0]),
 ("Nup82.1", "Nup82.1_2",   0.0,   fasta_files+"Nup82.txt",   "Nup82",   datadirectory+"cc_tr1_1.pdb", "E",  (522,612,0),    None,         beadsize,       2,      [11,2],     4,               None,            None, [0]),
 ("Nup82.1", "Nup82.1_21",  0.0,   fasta_files+"Nup82.txt",   "Nup82",   "BEADS",                      " ",  (613,624,0),    None,         beadsize,       102,    [11,102],   0,               None,            None, [0]),
 ("Nup82.1", "Nup82.1_3",   0.0,   fasta_files+"Nup82.txt",   "Nup82",   datadirectory+"cc_tr2_1.pdb", "E",  (625,669,0),    None,         beadsize,       3,      [11,3],     2,               None,            None, [0]),
 ("Nup82.1", "Nup82.1_31",  0.0,   fasta_files+"Nup82.txt",   "Nup82",   "BEADS",                      " ",  (670,677,0),    None,         beadsize,       103,    [11,103],   0,               None,            None, [0]),
 ("Nup82.1", "Nup82.1_4",   0.0,   fasta_files+"Nup82.txt",   "Nup82",   datadirectory+"cc_tr3_1.pdb", "E",  (678,713,0),    None,         beadsize,       4,      [11,4],     2,               None,            None, [0]),

 ("Nup82.2", "Nup82.2_1",   0.15,  fasta_files+"Nup82.txt",   "Nup82",   datadirectory+"3PBP2.pdb",    "G",  (  1,452,0),    None,         beadsize,       5,      [15,5],     18,              None,            None, [4]),
 ("Nup82.2", "Nup82.2_11",  0.15,  fasta_files+"Nup82.txt",   "Nup82",   "BEADS",                      " ",  (453,521,0),    None,         beadsize,       105,    [15,105],   0,               None,            None, [4]),
 ("Nup82.2", "Nup82.2_2",   0.15,  fasta_files+"Nup82.txt",   "Nup82",   datadirectory+"cc_tr1_2.pdb", "E",  (522,612,0),    None,         beadsize,       6,      [15,6],     4,               None,            None, [4]),
 ("Nup82.2", "Nup82.2_21",  0.15,  fasta_files+"Nup82.txt",   "Nup82",   "BEADS",                      " ",  (613,624,0),    None,         beadsize,       106,    [15,106],   0,               None,            None, [4]),
 ("Nup82.2", "Nup82.2_3",   0.15,  fasta_files+"Nup82.txt",   "Nup82",   datadirectory+"cc_tr2_2.pdb", "E",  (625,669,0),    None,         beadsize,       7,      [15,7],     2,               None,            None, [4]),
 ("Nup82.2", "Nup82.2_31",  0.15,  fasta_files+"Nup82.txt",   "Nup82",   "BEADS",                      " ",  (670,677,0),    None,         beadsize,       107,    [15,107],   0,               None,            None, [4]),
 ("Nup82.2", "Nup82.2_4",   0.15,  fasta_files+"Nup82.txt",   "Nup82",   datadirectory+"cc_tr3_2.pdb", "E",  (678,713,0),    None,         beadsize,       8,      [15,8],     2,               None,            None, [4]),


 ("Nup159.1", "Nup159.1_1", 1.0,   fasta_files+"Nup159.txt",  "Nup159",  datadirectory+"1XIP1.pdb",    "A",  (   1, 381,0),  None,         beadsize,       98,     [12,98],    0,               None,            None, [2]),
 ("Nup159.1", "Nup159.1_11",1.0,   fasta_files+"Nup159.txt",  "Nup159",  "BEADS",                      " ",  ( 382,1116,0),  None,         beadsize100,    198,    [12,198],   0,               None,            None, [2]),
 ("Nup159.1", "Nup159.1_2", 1.0,   fasta_files+"Nup159.txt",  "Nup159",  datadirectory+"4DS1.pdb",     "B",  (1117,1126,0),  None,         beadsize,       0,      [12,0],     1,               None,            None, [2]),
 ("Nup159.1", "Nup159.1_21",1.0,   fasta_files+"Nup159.txt",  "Nup159",  "BEADS",                      " ",  (1127,1210,0),  None,         beadsize,       110,    [12,110],   0,               None,            None, [2]),
 ("Nup159.1", "Nup159.1_3", 1.0,   fasta_files+"Nup159.txt",  "Nup159",  datadirectory+"cc_tr1_1.pdb", "D",  (1211,1321,0),  None,         beadsize,       2,      [12,2],     4,               None,            None, [2]),
 ("Nup159.1", "Nup159.1_31",1.0,   fasta_files+"Nup159.txt",  "Nup159",  "BEADS",                      " ",  (1322,1331,0),  None,         beadsize,       112,    [12,112],   0,               None,            None, [2]),
 ("Nup159.1", "Nup159.1_4", 1.0,   fasta_files+"Nup159.txt",  "Nup159",  datadirectory+"cc_tr2_1.pdb", "D",  (1332,1372,0),  None,         beadsize,       3,      [12,3],     2,               None,            None, [2]),
 ("Nup159.1", "Nup159.1_41",1.0,   fasta_files+"Nup159.txt",  "Nup159",  "BEADS",                      " ",  (1373,1381,0),  None,         beadsize,       113,    [12,113],   0,               None,            None, [2]),
 ("Nup159.1", "Nup159.1_5", 1.0,   fasta_files+"Nup159.txt",  "Nup159",  datadirectory+"cc_tr3_1.pdb", "D",  (1382,1412,0),  None,         beadsize,       4,      [12,4],     2,               None,            None, [2]),
 ("Nup159.1", "Nup159.1_51",1.0,   fasta_files+"Nup159.txt",  "Nup159",  "BEADS",                      " ",  (1413,1428,0),  None,         beadsize,       114,    [12,114],   0,               None,            None, [2]),
 ("Nup159.1", "Nup159.1_6", 1.0,   fasta_files+"Nup159.txt",  "Nup159",  datadirectory+"3PBP1.pdb",    "I",  (1429,1460,0),  None,         beadsize,       1,      [12,1],     2,               None,            None, [2]),

 ("Nup159.2", "Nup159.2_1", 0.9,   fasta_files+"Nup159.txt",  "Nup159",  datadirectory+"1XIP2.pdb",    "A",  (   1, 381,0),  None,         beadsize,       99,     [16,99],    0,               None,            None, [8]),
 ("Nup159.2", "Nup159.2_11",0.9,   fasta_files+"Nup159.txt",  "Nup159",  "BEADS",                      " ",  ( 382,1116,0),  None,         beadsize100,    199,    [16,199],   0,               None,            None, [8]),
 ("Nup159.2", "Nup159.2_2", 0.9,   fasta_files+"Nup159.txt",  "Nup159",  datadirectory+"4DS1.pdb",     "D",  (1117,1126,0),  None,         beadsize,       0,      [16,0],     1,               None,            None, [8]),
 ("Nup159.2", "Nup159.2_21",0.9,   fasta_files+"Nup159.txt",  "Nup159",  "BEADS",                      " ",  (1127,1210,0),  None,         beadsize,       115,    [16,115],   0,               None,            None, [8]),
 ("Nup159.2", "Nup159.2_3", 0.9,   fasta_files+"Nup159.txt",  "Nup159",  datadirectory+"cc_tr1_2.pdb", "D",  (1211,1321,0),  None,         beadsize,       6,      [16,6],     4,               None,            None, [8]),
 ("Nup159.2", "Nup159.2_31",0.9,   fasta_files+"Nup159.txt",  "Nup159",  "BEADS",                      " ",  (1322,1331,0),  None,         beadsize,       116,    [16,116],   0,               None,            None, [8]),
 ("Nup159.2", "Nup159.2_4", 0.9,   fasta_files+"Nup159.txt",  "Nup159",  datadirectory+"cc_tr2_2.pdb", "D",  (1332,1372,0),  None,         beadsize,       7,      [16,7],     2,               None,            None, [8]),
 ("Nup159.2", "Nup159.2_41",0.9,   fasta_files+"Nup159.txt",  "Nup159",  "BEADS",                      " ",  (1373,1381,0),  None,         beadsize,       117,    [16,117],   0,               None,            None, [8]),
 ("Nup159.2", "Nup159.2_5", 0.9,   fasta_files+"Nup159.txt",  "Nup159",  datadirectory+"cc_tr3_2.pdb", "D",  (1382,1412,0),  None,         beadsize,       8,      [16,8],     2,               None,            None, [8]),
 ("Nup159.2", "Nup159.2_51",0.9,   fasta_files+"Nup159.txt",  "Nup159",  "BEADS",                      " ",  (1413,1428,0),  None,         beadsize,       118,    [16,118],   0,               None,            None, [8]),
 ("Nup159.2", "Nup159.2_6", 0.9,   fasta_files+"Nup159.txt",  "Nup159",  datadirectory+"3PBP2.pdb",    "I",  (1429,1460,0),  None,         beadsize,       5,      [16,5],     2,               None,            None, [8]),


 ("Nsp1.1",  "Nsp1.1_10",   0.3,   fasta_files+"Nsp1.txt",    "Nsp1",    "BEADS",                      " ",  (  1,636,0),    None,         beadsize100,    161,    [13,161],   0,               None,            None, [1]),
 ("Nsp1.1",  "Nsp1.1_1",    0.3,   fasta_files+"Nsp1.txt",    "Nsp1",    datadirectory+"cc_tr1_1.pdb", "C",  (637,727,0),    None,         beadsize,       2,      [13,2],     4,               None,            None, [1]),
 ("Nsp1.1",  "Nsp1.1_11",   0.3,   fasta_files+"Nsp1.txt",    "Nsp1",    "BEADS",                      " ",  (728,741,0),    None,         beadsize,       162,    [13,162],   0,               None,            None, [1]),
 ("Nsp1.1",  "Nsp1.1_2",    0.3,   fasta_files+"Nsp1.txt",    "Nsp1",    datadirectory+"cc_tr2_1.pdb", "C",  (742,778,0),    None,         beadsize,       3,      [13,3],     2,               None,            None, [1]),
 ("Nsp1.1",  "Nsp1.1_21",   0.3,   fasta_files+"Nsp1.txt",    "Nsp1",    "BEADS",                      " ",  (779,787,0),    None,         beadsize,       163,    [13,163],   0,               None,            None, [1]),
 ("Nsp1.1",  "Nsp1.1_3",    0.3,   fasta_files+"Nsp1.txt",    "Nsp1",    datadirectory+"cc_tr3_1.pdb", "C",  (788,823,0),    None,         beadsize,       4,      [13,4],     2,               None,            None, [1]),

 ("Nsp1.2",  "Nsp1.2_10",   0.38,  fasta_files+"Nsp1.txt",    "Nsp1",    "BEADS",                      " ",  (  1,636,0),    None,         beadsize100,    165,    [17,165],   0,               None,            None, [7]),
 ("Nsp1.2",  "Nsp1.2_1",    0.38,  fasta_files+"Nsp1.txt",    "Nsp1",    datadirectory+"cc_tr1_2.pdb", "C",  (637,727,0),    None,         beadsize,       6,      [17,6],     4,               None,            None, [7]),
 ("Nsp1.2",  "Nsp1.2_11",   0.38,  fasta_files+"Nsp1.txt",    "Nsp1",    "BEADS",                      " ",  (728,741,0),    None,         beadsize,       166,    [17,166],   0,               None,            None, [7]),
 ("Nsp1.2",  "Nsp1.2_2",    0.38,  fasta_files+"Nsp1.txt",    "Nsp1",    datadirectory+"cc_tr2_2.pdb", "C",  (742,778,0),    None,         beadsize,       7,      [17,7],     2,               None,            None, [7]),
 ("Nsp1.2",  "Nsp1.2_21",   0.38,  fasta_files+"Nsp1.txt",    "Nsp1",    "BEADS",                      " ",  (779,787,0),    None,         beadsize,       167,    [17,167],   0,               None,            None, [7]),
 ("Nsp1.2",  "Nsp1.2_3",    0.38,  fasta_files+"Nsp1.txt",    "Nsp1",    datadirectory+"cc_tr3_2.pdb", "C",  (788,823,0),    None,         beadsize,       8,      [17,8],     2,               None,            None, [7]),


 ("Nup116.1","Nup116.1_10", 0.75,  fasta_files+"Nup116.txt",  "Nup116",  "BEADS",                      " ",  (  1, 750,0),   None,         beadsize100,    191,    [14,191],   0,               None,            None, [5]),
 ("Nup116.1","Nup116.1_1",  0.75,  fasta_files+"Nup116.txt",  "Nup116",  datadirectory+"3PBP1.pdb",    "H",  (751,1113,0),   None,         beadsize25,     1,      [14,1],     0,               None,            None, [5]),

 ("Nup116.2","Nup116.2_10", 0.8,   fasta_files+"Nup116.txt",  "Nup116",  "BEADS",                      " ",  (  1, 750,0),   None,         beadsize100,    195,    [18,195],   0,               None,            None, [6]),
 ("Nup116.2","Nup116.2_1",  0.8,   fasta_files+"Nup116.txt",  "Nup116",  datadirectory+"3PBP2.pdb",    "H",  (751,1113,0),   None,         beadsize25,     5,      [18,5],     0,               None,            None, [6])
]

bm1 = IMP.pmi1.macros.BuildModel1(simo)
#bm1.set_gmm_models_directory(datadirectory + "em_gmm_model/")

if (inputs.rmf_input is not None) :
    n82=set([s[0] for s in domains])
    for d in list(n82):
        bm1.set_rmf_file(d, inputs.rmf_input, int(inputs.rmf_frame_number))

bm1.build_model(data_structure = domains, sequence_connectivity_scale=0.5, sequence_connectivity_resolution=1.0)
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
    #simo.shuffle_configuration(50)
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
print("ExcludedVolumeSphere !!\n")


#####################################################
# Restraints setup
# External Barrier restraint
#####################################################
eb = IMP.pmi1.restraints.basic.ExternalBarrier(simo, radius = 300)
eb.add_to_model()
outputobjects.append(eb)
print(eb.get_output())
print("ExternalBarrier !!\n")


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
    xl2 = IMP.pmi1.restraints.crosslinking.ISDCrossLinkMS(simo,
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
    xl3 = IMP.pmi1.restraints.crosslinking.ISDCrossLinkMS(simo,
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

    xl4 = IMP.pmi1.restraints.crosslinking.ISDCrossLinkMS(simo,
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

    dr1 = IMP.pmi1.restraints.basic.DistanceRestraint(simo,(1417,1417,"Nup159.1"), (1417,1417,"Nup159.2"), distancemin=dist_min, distancemax=dist_max, resolution=1.0, label="Nup159_1417-1417")
    dr1.add_to_model()
    dr1.set_weight(dr_weight)
    outputobjects.append(dr1)
    print(dr1.get_output())

    dr2 = IMP.pmi1.restraints.basic.DistanceRestraint(simo,(1432,1432,"Nup159.1"), (1432,1432,"Nup159.2"), distancemin=dist_min, distancemax=dist_max, resolution=1.0, label="Nup159_1432-1432")
    dr2.add_to_model()
    dr2.set_weight(dr_weight)
    outputobjects.append(dr2)
    print(dr2.get_output())

    dr3 = IMP.pmi1.restraints.basic.DistanceRestraint(simo,(1384,1384,"Nup159.1"), (1384,1384,"Nup159.2"), distancemin=dist_min, distancemax=dist_max, resolution=1.0, label="Nup159_1384-1384")
    dr3.add_to_model()
    dr3.set_weight(dr_weight)
    outputobjects.append(dr3)
    print(dr3.get_output())

    dr4 = IMP.pmi1.restraints.basic.DistanceRestraint(simo,(1414,1414,"Nup159.1"), (1414,1414,"Nup159.2"), distancemin=dist_min, distancemax=dist_max, resolution=1.0, label="Nup159_1414-1414")
    dr4.add_to_model()
    dr4.set_weight(dr_weight)
    outputobjects.append(dr4)
    print(dr4.get_output())

    dr5 = IMP.pmi1.restraints.basic.DistanceRestraint(simo,(1387,1387,"Nup159.1"), (1387,1387,"Nup159.2"), distancemin=dist_min, distancemax=dist_max, resolution=1.0, label="Nup159_1387-1387")
    dr5.add_to_model()
    dr5.set_weight(dr_weight)
    outputobjects.append(dr5)
    print(dr5.get_output())

    dr6 = IMP.pmi1.restraints.basic.DistanceRestraint(simo,(517,517,"Nup82.1"), (517,517,"Nup82.2"), distancemin=dist_min, distancemax=dist_max, resolution=1.0, label="Nup82_517-517")
    dr6.add_to_model()
    dr6.set_weight(dr_weight)
    outputobjects.append(dr6)
    print(dr6.get_output())

    # by Ed Hurt
    if (False):
        """
        dr21 = IMP.pmi1.restraints.basic.DistanceRestraint(simo,(1343,1343,"Nup159.1"), (1343,1343,"Nup159.2"), distancemin=dist_min, distancemax=dist_max, resolution=1.0, label="Nup159_1343-1343")
        dr21.add_to_model()
        dr21.set_weight(dr_weight)
        outputobjects.append(dr21)
        print(dr21.get_output())

        dr22 = IMP.pmi1.restraints.basic.DistanceRestraint(simo,(541,541,"Nup82.1"), (541,541,"Nup82.2"), distancemin=dist_min, distancemax=dist_max, resolution=1.0, label="Nup82_541-541")
        dr22.add_to_model()
        dr22.set_weight(dr_weight)
        outputobjects.append(dr22)
        print(dr22.get_output())

        d23 = IMP.pmi1.restraints.basic.DistanceRestraint(simo,(580,580,"Nup82.1"), (580,580,"Nup82.2"), distancemin=dist_min, distancemax=dist_max, resolution=1.0, label="Nup82_580-580")
        d23.add_to_model()
        d23.set_weight(dr_weight)
        outputobjects.append(d23)
        print(d23.get_output())
        """

    print("\nDistance Restraints applied for homo-dimers !!")
    print("weight = ", dr_weight, "dist_min = ", dist_min, "dist_max = ", dist_max, "\n")

    # for skNup82 homo-dimers - DISABLED
    if (False):
        """
        dist_min = 3.0
        dist_max = 35.0
        dr_weight = 10.0

        dr7 = IMP.pmi1.restraints.basic.DistanceRestraint(simo,(649,649,"Nup82.1"), (692,692,"Nup82.2"), distancemin=dist_min, distancemax=dist_max, resolution=1.0, label="Nup82_649-692a")
        dr7.add_to_model()
        dr7.set_weight(dr_weight)
        outputobjects.append(dr7)
        print(dr7.get_output())

        dr8 = IMP.pmi1.restraints.basic.DistanceRestraint(simo,(649,649,"Nup82.2"), (692,692,"Nup82.1"), distancemin=dist_min, distancemax=dist_max, resolution=1.0, label="Nup82_649-692b")
        dr8.add_to_model()
        dr8.set_weight(dr_weight)
        outputobjects.append(dr8)
        print(dr8.get_output())

        dr9 = IMP.pmi1.restraints.basic.DistanceRestraint(simo,(541,541,"Nup82.1"), (569,569,"Nup82.2"), distancemin=dist_min, distancemax=dist_max, resolution=1.0, label="Nup82_649-692b")
        dr9.add_to_model()
        dr9.set_weight(dr_weight)
        outputobjects.append(dr9)
        print(dr9.get_output())

        dr10 = IMP.pmi1.restraints.basic.DistanceRestraint(simo,(541,541,"Nup82.2"), (569,569,"Nup82.1"), distancemin=dist_min, distancemax=dist_max, resolution=1.0, label="Nup82_541-569b")
        dr10.add_to_model()
        dr10.set_weight(dr_weight)
        outputobjects.append(dr10)
        print(dr10.get_output())

        print("\nDistance Restraints applied for skNup82 homo-dimers !!")
        print("weight = ", dr_weight, "dist_min = ", dist_min, "dist_max = ", dist_max, "\n")
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

    db = IMP.pmi1.tools.get_db_from_csv('../data/XL_cliques_2copies.csv')

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

        temp_label = mol1 + "_" + str(res1) + "-" + mol2 + "_" + str(res2)
        dr = IMP.pmi1.restraints.basic.DistanceRestraint(simo,(res1,res1,mol1), (res2,res2,mol2), distancemin=dist_min, distancemax=dist_max, resolution=1.0, label=temp_label)
        dr.add_to_model()
        dr.set_weight(dr_weight)
        outputobjects.append(dr)
        print(dr.get_output())

    print("\nDistance Restraints applied for XL cliques !!")
    print("weight = ", dr_weight, "dist_min = ", dist_min, "dist_max = ", dist_max, "\n")



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
print("\nEVAL 1 : ", sf.evaluate(False), " (initial) - ", rank)

if (True):
    if not inputs.dry_run:
        simo.optimize_floppy_bodies(150)
        print("\nEVAL 2 : ", sf.evaluate(False), " (after calling optimize_floppy_bodies(150)) - ", rank)

    initial_nframes = 1000
    mc1 = IMP.pmi1.macros.ReplicaExchange0(m,
                                        simo,
                                        monte_carlo_sample_objects = sampleobjects,
                                        output_objects = outputobjects,
                                        crosslink_restraints = [xl1, xl2, xl3],
                                        #crosslink_restraints = [xl4],
                                        monte_carlo_temperature = 1.0,
                                        replica_exchange_minimum_temperature = 1.0,
                                        replica_exchange_maximum_temperature = 2.5,
                                        number_of_best_scoring_models = 100,
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
                                        test_mode=simo.dry_run,
                                        replica_stat_file_suffix = "stat_replica")
    mc1.execute_macro()
    rex1 = mc1.get_replica_exchange_object()
    print("\nEVAL 3 : ", sf.evaluate(False), " (after performing the pre-sampling) - ", rank)
else:
    rex1 = None
    print("\n>> NO pre-sampling")


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
                                            image_resolution = 35.0,
                                            projection_number = 100,
                                            #projection_number = 400,
                                            n_components = 2)
    em2d.add_to_model()
    em2d.set_weight(em2d_weight)
    outputobjects.append(em2d)

    sf = IMP.core.RestraintsScoringFunction(IMP.pmi1.tools.get_restraint_set(m))
    print("\nEVAL 4 : ", sf.evaluate(False), " (after applying the EM 2D restraint) - ", rank)


#####################################################
# Metropolis Monte Carlo sampling with Replica Exchange
#####################################################
mc2 = IMP.pmi1.macros.ReplicaExchange0(m,
                                    simo,
                                    monte_carlo_sample_objects = sampleobjects,
                                    output_objects = outputobjects,
                                    crosslink_restraints = [xl1, xl2, xl3],
                                    #crosslink_restraints = [xl4],
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
                                    test_mode=simo.dry_run,
                                    replica_exchange_object = rex1)
mc2.execute_macro()
print("\nEVAL 5 : ", sf.evaluate(False), " (final evaluation) - ", rank)

if inputs.mmcif:
    # Link entities to UniProt
    if hasattr(ihm, 'reference'):
        for subunit, accession in (
                ('Dyn2.1', 'Q02647'), ('Nup82.1', 'P40368'),
                ('Nup159.1', 'P40477'), ('Nsp1.1', 'P14907'),
                ('Nup116.1', 'Q02630')):
            ref = ihm.reference.UniProtSequence.from_accession(accession)
            e = po.asym_units[subunit].entity.references.append(ref)

    # Add refinement step with all 2D EM data
    # (see 2_modeling_allEM_except11_19.py)
    # note that we also exclude image #2 because we've already added that above
    images = ["../data/em2d/%d.pgm" % i for i in range(23)
                                        if i not in (2,11,19)]
    em2d = em2d_nup82.ElectronMicroscopy2D(simo, images, resolution = 1.0,
                                           pixel_size = 3.23,
                                           image_resolution = 35.0,
                                           projection_number = 100,
                                           n_components = 2)
    em2d.add_to_model()
    em2d.set_weight(em2d_weight)
    outputobjects.append(em2d)
    sf = IMP.core.RestraintsScoringFunction(IMP.pmi1.tools.get_restraint_set(m))
    mc3 = IMP.pmi1.macros.ReplicaExchange0(m, simo,
                                    monte_carlo_sample_objects = sampleobjects,
                                    output_objects = outputobjects,
                                    crosslink_restraints = [xl1, xl2, xl3],
                                    monte_carlo_temperature = 1.0,
                                    replica_exchange_minimum_temperature = 1.0,
                                    replica_exchange_maximum_temperature = 2.5,
                                    number_of_best_scoring_models = 500,
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
                                    test_mode=simo.dry_run,
                                    replica_exchange_object = rex1)
    mc3.execute_macro()

    # Read in final model
    for c in simo.get_component_names():
        simo.set_coordinates_from_rmf(c,
                '../outputs/4_em2d_single_scores_final35r_best/225-17.rmf3', 0,
                force_rigid_update=True)
    pp = po._add_simple_postprocessing(num_models_begin=10000,
                                       num_models_end=463)

    r = ihm.location.Repository(doi=doi, url="%s/cluster0.dcd" % url_top)
    f = ihm.location.OutputFileLocation(path='.', repo=r,
                details="All ensemble structures for cluster 0")
    # precision of 9.0 is from Table 1
    e = po._add_simple_ensemble(pp, name="Cluster 0", num_models=370,
                                drmsd=9.0, num_models_deposited=1,
                                localization_densities={}, ensemble_file=f)

    # Add localization densities
    for fname, domains in [
             ('Dyn2', [po.asym_units['Dyn2.1'], po.asym_units['Dyn2.2']]),
             ('Nup82.1_NTD', [po.asym_units['Nup82.1'](1,521)]),
             ('Nup82.1_CTD', [po.asym_units['Nup82.1'](522,713)]),
             ('Nup82.2_NTD', [po.asym_units['Nup82.2'](1,521)]),
             ('Nup82.2_CTD', [po.asym_units['Nup82.2'](522,713)]),
             ('Nup159.1_CTD', [po.asym_units['Nup159.1'](1117,1460)]),
             ('Nup159.2_CTD', [po.asym_units['Nup159.2'](1117,1460)]),
             ('Nsp1.1_CTD', [po.asym_units['Nsp1.1'](637,823)]),
             ('Nsp1.2_CTD', [po.asym_units['Nsp1.1'](637,823)]) ]:
        loc = ihm.location.OutputFileLocation(
                '../outputs/3_analysis_allEM2D/kmeans_10000_2/'
                'cluster.0_370/%s.mrc' % fname)
        for domain in domains:
            den = ihm.model.LocalizationDensity(file=loc, asym_unit=domain)
            e.densities.append(den)

    model = po.add_model(e.model_group)
    model._is_restrained = False # We have no restraint data for this model

    # Correct number of output models to account for multiple runs
    # todo: step 1 isn't quite right because it shows modeling with a single
    # em2D image; in fact we ran multiple "step 1"s, each with a different image
    model.protocol.steps[1].num_models_end = 1350000
    # todo: add em2D filter between steps 1 and 2
    model.protocol.steps[2].num_models_begin = 1350000
    model.protocol.steps[2].num_models_end = 10000

    # Add SAXS datasets (for validation)
    f = saxs.SAXSFits(po)
    validation_datasets = list(f.add_from_csv(model))

    # Add validation step to the protocol
    analysis = ihm.analysis.Analysis()
    analysis.steps.append(ihm.analysis.ValidationStep(
                   assembly=po.system.complete_assembly,
                   dataset_group=ihm.dataset.DatasetGroup(validation_datasets),
                   feature='other', num_models_begin=1, num_models_end=1))
    model.protocol.analyses.append(analysis)

    # Correct crosslinker type from skDSS to DSS
    for r in po.system.restraints:
        if hasattr(r, 'linker') and r.linker.auth_name == 'skDSS':
            r.linker = ihm.cross_linkers.dss

    po.flush()
