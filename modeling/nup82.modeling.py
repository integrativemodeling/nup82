import IMP
import IMP.core
import IMP.base
import IMP.algebra
import IMP.atom
import IMP.container

import IMP.pmi.restraints.crosslinking
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.em
import IMP.pmi.restraints.basic
import IMP.pmi.representation
import IMP.pmi.tools
import IMP.pmi.samplers
import IMP.pmi.output
import IMP.pmi.macros

import os

# setting up parameters

rbmaxtrans = 2.00
fbmaxtrans = 3.00
rbmaxrot=0.04
outputobjects = []
sampleobjects = []
beadsize=20

# setting up topology

m = IMP.Model()
simo = IMP.pmi.representation.Representation(m,upperharmonic=True,disorderedlength=False)

datadirectory="../data/"


# compname         hier_name     color    fastafile                      fastaid  pdbname                   chain resrange    read                     beadsize rigid_body super_rigid_body emnum_components, emtxtfilename    emmrcfilename                    

domains=[("Nup82",  "Nup82_1",     0.0,     "../data/protein_fasta.Nup82.txt",   "Nup82",   datadirectory+"/3PBP.pdb",          "A",  (1,521,0),      None,                   beadsize,       0,         [0],       0,                None,            None),
         ("Nup82",  "Nup82_2",     0.0,     "../data/protein_fasta.Nup82.txt",   "Nup82",   "IDEAL_HELIX",                      " ",  (522,571,0),    None,                   beadsize,       1,         [0,1],     0,                None,            None),
         ("Nup82",  "Nup82_3",     0.0,     "../data/protein_fasta.Nup82.txt",   "Nup82",   "BEADS",                            " ",  (572,588,0),    None,                   beadsize,       2,         [0,2],     0,                None,            None),
         ("Nup82",  "Nup82_4",     0.0,     "../data/protein_fasta.Nup82.txt",   "Nup82",   "IDEAL_HELIX",                      " ",  (589,613,0),    None,                   beadsize,       3,         [0,3],     0,                None,            None),
         ("Nup82",  "Nup82_5",     0.0,     "../data/protein_fasta.Nup82.txt",   "Nup82",   "BEADS",                            " ",  (614,624,0),    None,                   beadsize,       4,         [0,4],     0,                None,            None),
         ("Nup82",  "Nup82_6_1",   0.0,     "../data/protein_fasta.Nup82.txt",   "Nup82",   "IDEAL_HELIX",                      " ",  (625,641,0),    None,                   beadsize,       5,         [0,5],     0,                None,            None),
         ("Nup82",  "Nup82_6_2",   0.0,     "../data/protein_fasta.Nup82.txt",   "Nup82",   "BEADS",                            " ",  (642,647,0),    None,                   beadsize,       51,        [0,51],     0,                None,            None),
         ("Nup82",  "Nup82_6_3",   0.0,     "../data/protein_fasta.Nup82.txt",   "Nup82",   "IDEAL_HELIX",                      " ",  (648,671,0),    None,                   beadsize,       52,        [0,52],     0,                None,            None),
         ("Nup82",  "Nup82_7",     0.0,     "../data/protein_fasta.Nup82.txt",   "Nup82",   "BEADS",                            " ",  (672,676,0),    None,                   beadsize,       6,         [0,6],     0,                None,            None),
         ("Nup82",  "Nup82_8",     0.0,     "../data/protein_fasta.Nup82.txt",   "Nup82",   "IDEAL_HELIX",                      " ",  (677,709,0),    None,                   beadsize,       7,         [0,7],     0,                None,            None),
         ("Nup82",  "Nup82_9",     0.0,     "../data/protein_fasta.Nup82.txt",   "Nup82",   "BEADS",                            " ",  (710,713,0),    None,                   beadsize,       8,         [0,8],     0,                None,            None),
         ("Nsp1",   "Nsp1_1",      0.55,     "../data/protein_fasta.Nsp1.txt",    "Nsp1",    "BEADS",                            " ",  (1,649,0),      None,                   beadsize,       9,         [9,10],     0,                None,            None),
         ("Nsp1",   "Nsp1_2",      0.55,     "../data/protein_fasta.Nsp1.txt",    "Nsp1",    "IDEAL_HELIX",                      " ",  (650,671,0),    None,                   beadsize,       10,        [9,11],     0,                None,            None),
         ("Nsp1",   "Nsp1_3",      0.55,     "../data/protein_fasta.Nsp1.txt",    "Nsp1",    "BEADS",                            " ",  (672,672,0),    None,                   beadsize,       11,        [9,12],     0,                None,            None),
         ("Nsp1",   "Nsp1_4",      0.55,     "../data/protein_fasta.Nsp1.txt",    "Nsp1",    "IDEAL_HELIX",                      " ",  (673,725,0),    None,                   beadsize,       12,        [9,13],     0,                None,            None),
         ("Nsp1",   "Nsp1_5",      0.55,     "../data/protein_fasta.Nsp1.txt",    "Nsp1",    "BEADS",                            " ",  (726,737,0),    None,                   beadsize,       13,        [9,14],     0,                None,            None),
         ("Nsp1",   "Nsp1_6",      0.55,     "../data/protein_fasta.Nsp1.txt",    "Nsp1",    "IDEAL_HELIX",                      " ",  (738,771,0),    None,                   beadsize,       14,        [9,15],     0,                None,            None),
         ("Nsp1",   "Nsp1_7",      0.55,     "../data/protein_fasta.Nsp1.txt",    "Nsp1",    "BEADS",                            " ",  (772,788,0),    None,                   beadsize,       15,        [9,16],     0,                None,            None),
         ("Nsp1",   "Nsp1_8",      0.55,     "../data/protein_fasta.Nsp1.txt",    "Nsp1",    "IDEAL_HELIX",                      " ",  (789,823,0),    None,                   beadsize,       16,        [9,17],     0,                None,            None),
         ("Nup159", "Nup159_1",    1.0,     "../data/protein_fasta.Nup159.txt",  "Nup159",  datadirectory+"1XIP.pdb",           "A",  (1,521,0),      None,                   beadsize,       17,        [18,19],     0,                None,            None),
         ("Nup159", "Nup159_2",    1.0,     "../data/protein_fasta.Nup159.txt",  "Nup159",  "BEADS",                            " ",  (522,1099,0),   None,                   beadsize,       18,        [18,20],     0,                None,            None),
         ("Nup159", "Nup159_3",    1.0,     "../data/protein_fasta.Nup159.txt",  "Nup159",  datadirectory+"4DS1.pdb",           "B",  (1100,1146,0),  None,                   beadsize,       19,        [18,21],     0,                None,            None),
         ("Nup159", "Nup159_4",    1.0,     "../data/protein_fasta.Nup159.txt",  "Nup159",  "BEADS",                            " ",  (1147,1208,0),  None,                   beadsize,       20,        [18,22],     0,                None,            None),
         ("Nup159", "Nup159_5",    1.0,     "../data/protein_fasta.Nup159.txt",  "Nup159",  "IDEAL_HELIX",                      " ",  (1209,1241,0),  None,                   beadsize,       21,        [18,23],     0,                None,            None),
         ("Nup159", "Nup159_6",    1.0,     "../data/protein_fasta.Nup159.txt",  "Nup159",  "BEADS",                            " ",  (1242,1265,0),  None,                   beadsize,       22,        [18,24],     0,                None,            None),
         ("Nup159", "Nup159_7",    1.0,     "../data/protein_fasta.Nup159.txt",  "Nup159",  "IDEAL_HELIX",                      " ",  (1266,1320,0),  None,                   beadsize,       23,        [18,25],     0,                None,            None),
         ("Nup159", "Nup159_8",    1.0,     "../data/protein_fasta.Nup159.txt",  "Nup159",  "BEADS",                            " ",  (1321,1334,0),  None,                   beadsize,       24,        [18,26],     0,                None,            None),
         ("Nup159", "Nup159_9",    1.0,     "../data/protein_fasta.Nup159.txt",  "Nup159",  "IDEAL_HELIX",                      " ",  (1335,1369,0),  None,                   beadsize,       25,        [18,27],     0,                None,            None),
         ("Nup159", "Nup159_10",   1.0,     "../data/protein_fasta.Nup159.txt",  "Nup159",  "BEADS",                            " ",  (1370,1379,0),  None,                   beadsize,       26,        [18,28],     0,                None,            None),
         ("Nup159", "Nup159_11",   1.0,     "../data/protein_fasta.Nup159.txt",  "Nup159",  "IDEAL_HELIX",                      " ",  (1380,1410,0),  None,                   beadsize,       27,        [18,29],     0,                None,            None),
         ("Nup159", "Nup159_12",   1.0,     "../data/protein_fasta.Nup159.txt",  "Nup159",  "BEADS",                            " ",  (1411,1432,0),  None,                   beadsize,       28,        [18,30],     0,                None,            None),
         ("Nup159", "Nup159_13",   1.0,     "../data/protein_fasta.Nup159.txt",  "Nup159",  datadirectory+"/3PBP.pdb",          "C",  (1433,1460,0),  None,                   beadsize,       0,         [18,31],     0,                None,            None),
         ("Nup116", "Nup116",      0.75,    "../data/protein_fasta.Nup116.txt",  "YMR047C", datadirectory+"/3PBP.pdb",          "B",  (1,-1,0),       None,                   beadsize,       0,         [0],         0,                    None,            None),
         ("Dyn2",   "Dyn2",        0.25,    "../data/protein_fasta.Dyn2.txt",    "Dyn2",    datadirectory+'/4DS1.pdb',          "A",  (1,-1,0),       None,                   beadsize,       19,        [18,32],     0,                None,            None)]


bm1=IMP.pmi.macros.BuildModel1(simo)
bm1.build_model(domains)
resdensities=bm1.get_density_hierarchies()


# randomize the initial configuration

simo.shuffle_configuration(100)

# defines the movers

simo.set_rigid_bodies_max_rot(rbmaxrot)
simo.set_floppy_bodies_max_trans(fbmaxtrans)
simo.set_rigid_bodies_max_trans(rbmaxtrans)

outputobjects.append(simo)
sampleobjects.append(simo)

# scoring function

ev = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(simo,resolution=10)
ev.add_to_model()
outputobjects.append(ev)



columnmap={}
columnmap["Protein1"]="Protein 1"
columnmap["Protein2"]="Protein 2"
columnmap["Residue1"]="Residue 1"
columnmap["Residue2"]="Residue 2"
columnmap["IDScore"]="p value"

ids_map=IMP.pmi.tools.map()
ids_map.set_map_element(1.0,1.0)

xl1 = IMP.pmi.restraints.crosslinking.ISDCrossLinkMS(simo,
                                   '../data/Publication_Nup82XL_list.csv',
                                   length=21.0,
                                   slope=0.01,
                                   columnmapping=columnmap,
                                   ids_map=ids_map,
                                   resolution=1.0,
                                   label="DSS",
                                   csvfile=True)
xl1.add_to_model()
sampleobjects.append(xl1)
outputobjects.append(xl1)
xl1.set_psi_is_sampled(False)
psi=xl1.get_psi(1.0)[0]
psi.set_scale(0.05)



# sampling



simo.optimize_floppy_bodies(100)

mc1=IMP.pmi.macros.ReplicaExchange0(m,
                                    simo,
                                    monte_carlo_sample_objects=sampleobjects,
                                    output_objects=outputobjects,
                                    crosslink_restraints=[xl1],
                                    monte_carlo_temperature=1.0,
                                    replica_exchange_minimum_temperature=1.0,
                                    replica_exchange_maximum_temperature=2.5,
                                    number_of_best_scoring_models=500,
                                    monte_carlo_steps=10,
                                    number_of_frames=50000,
                                    write_initial_rmf=True,
                                    initial_rmf_name_suffix="initial",
                                    stat_file_name_suffix="stat",
                                    best_pdb_name_suffix="model",
                                    do_clean_first=True,
                                    do_create_directories=True,
                                    global_output_directory="output",
                                    rmf_dir="rmfs/",
                                    best_pdb_dir="pdbs/",
                                    replica_stat_file_suffix="stat_replica")
mc1.execute_macro()



