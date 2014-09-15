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


#         compname hier_name     color    fastafile                      fastaid  pdbname                                       chain   resrange    read                     beadsize       rigid_body super_rigid_body emnum_components, emtxtfilename    emmrcfilename        chain of super rigid bodies            

domains=[("Nup82",  "Nup82_1",     0.0,     "../data/protein_fasta.Nup82.txt",   "Nup82",   datadirectory+"/3PBP.pdb",          "A",  (1,521,0),      True,                   beadsize,       0,         [100,0],       10,               None,           None, [0]),
         ("Nup82",  "Nup82_2",     0.0,     "../data/protein_fasta.Nup82.txt",   "Nup82",   "IDEAL_HELIX",                      " ",  (522,571,0),    True,                   beadsize,       1,         [100,0,1],     1,               None,            None, [0]),
         ("Nup82",  "Nup82_3",     0.0,     "../data/protein_fasta.Nup82.txt",   "Nup82",   "BEADS",                            " ",  (572,588,0),    True,                   beadsize,       2,         [100,0,2],     0,               None,            None, [0]),
         ("Nup82",  "Nup82_4",     0.0,     "../data/protein_fasta.Nup82.txt",   "Nup82",   "IDEAL_HELIX",                      " ",  (589,613,0),    True,                   beadsize,       3,         [100,0,3],     1,               None,            None, [0]),
         ("Nup82",  "Nup82_5",     0.0,     "../data/protein_fasta.Nup82.txt",   "Nup82",   "BEADS",                            " ",  (614,624,0),    True,                   beadsize,       4,         [100,0,4],     0,               None,            None, [0]),
         ("Nup82",  "Nup82_6_1",   0.0,     "../data/protein_fasta.Nup82.txt",   "Nup82",   "IDEAL_HELIX",                      " ",  (625,641,0),    True,                   beadsize,       5,         [100,0,5],     1,               None,            None, [0]),
         ("Nup82",  "Nup82_6_2",   0.0,     "../data/protein_fasta.Nup82.txt",   "Nup82",   "BEADS",                            " ",  (642,647,0),    True,                   beadsize,       51,        [100,0,51],     0,               None,            None, [0]),
         ("Nup82",  "Nup82_6_3",   0.0,     "../data/protein_fasta.Nup82.txt",   "Nup82",   "IDEAL_HELIX",                      " ",  (648,671,0),    True,                   beadsize,       52,        [100,0,52],     1,               None,            None, [0]),
         ("Nup82",  "Nup82_7",     0.0,     "../data/protein_fasta.Nup82.txt",   "Nup82",   "BEADS",                            " ",  (672,676,0),    True,                   beadsize,       6,         [100,0,6],     0,               None,            None, [0]),
         ("Nup82",  "Nup82_8",     0.0,     "../data/protein_fasta.Nup82.txt",   "Nup82",   "IDEAL_HELIX",                      " ",  (677,709,0),    True,                   beadsize,       7,         [100,0,7],     1,               None,            None, [0]),
         ("Nup82",  "Nup82_9",     0.0,     "../data/protein_fasta.Nup82.txt",   "Nup82",   "BEADS",                            " ",  (710,713,0),    True,                   beadsize,       8,         [100,0,8],     0,               None,            None, [0]),
         ("Nsp1",   "Nsp1_1",      0.55,     "../data/protein_fasta.Nsp1.txt",    "Nsp1",    "BEADS",                            " ",  (1,649,0),      None,                   beadsize,       9,        [100,9,10],     0,               None,            None, [1]),
         ("Nsp1",   "Nsp1_2",      0.55,     "../data/protein_fasta.Nsp1.txt",    "Nsp1",    "IDEAL_HELIX",                      " ",  (650,671,0),    True,                   beadsize,       10,       [100,9,11],     1,               None,            None, [1]),
         ("Nsp1",   "Nsp1_3",      0.55,     "../data/protein_fasta.Nsp1.txt",    "Nsp1",    "BEADS",                            " ",  (672,672,0),    True,                   beadsize,       11,       [100,9,12],     0,               None,            None, [1]),
         ("Nsp1",   "Nsp1_4",      0.55,     "../data/protein_fasta.Nsp1.txt",    "Nsp1",    "IDEAL_HELIX",                      " ",  (673,725,0),    True,                   beadsize,       12,       [100,9,13],     1,               None,            None, [1]),
         ("Nsp1",   "Nsp1_5",      0.55,     "../data/protein_fasta.Nsp1.txt",    "Nsp1",    "BEADS",                            " ",  (726,737,0),    True,                   beadsize,       13,       [100,9,14],     0,               None,            None, [1]),
         ("Nsp1",   "Nsp1_6",      0.55,     "../data/protein_fasta.Nsp1.txt",    "Nsp1",    "IDEAL_HELIX",                      " ",  (738,771,0),    True,                   beadsize,       14,       [100,9,15],     1,               None,            None, [1]),
         ("Nsp1",   "Nsp1_7",      0.55,     "../data/protein_fasta.Nsp1.txt",    "Nsp1",    "BEADS",                            " ",  (772,788,0),    True,                   beadsize,       15,       [100,9,16],     0,               None,            None, [1]),
         ("Nsp1",   "Nsp1_8",      0.55,     "../data/protein_fasta.Nsp1.txt",    "Nsp1",    "IDEAL_HELIX",                      " ",  (789,823,0),    True,                   beadsize,       16,       [100,9,17],     1,               None,            None, [1]),
         ("Nup159", "Nup159_1",    1.0,     "../data/protein_fasta.Nup159.txt",  "Nup159",  datadirectory+"1XIP.pdb",           "A",  (1,521,0),      None,                   beadsize,       17,        [100,18,19],     0,               None,            None, [2]),
         ("Nup159", "Nup159_2",    1.0,     "../data/protein_fasta.Nup159.txt",  "Nup159",  "BEADS",                            " ",  (522,1099,0),   None,                   beadsize,       18,        [100,18,20],     0,               None,            None, [2]),
         ("Nup159", "Nup159_3",    1.0,     "../data/protein_fasta.Nup159.txt",  "Nup159",  datadirectory+"4DS1.pdb",           "B",  (1100,1146,0),  None,                   beadsize,       19,        [100,18,21],     0,               None,            None, [2]),
         ("Nup159", "Nup159_4",    1.0,     "../data/protein_fasta.Nup159.txt",  "Nup159",  "BEADS",                            " ",  (1147,1208,0),  True,                   beadsize,       20,        [100,18,22],     0,               None,            None, [2]),
         ("Nup159", "Nup159_5",    1.0,     "../data/protein_fasta.Nup159.txt",  "Nup159",  "IDEAL_HELIX",                      " ",  (1209,1241,0),  True,                   beadsize,       21,        [100,18,23],     1,               None,            None, [2]),
         ("Nup159", "Nup159_6",    1.0,     "../data/protein_fasta.Nup159.txt",  "Nup159",  "BEADS",                            " ",  (1242,1265,0),  True,                   beadsize,       22,        [100,18,24],     0,               None,            None, [2]),
         ("Nup159", "Nup159_7",    1.0,     "../data/protein_fasta.Nup159.txt",  "Nup159",  "IDEAL_HELIX",                      " ",  (1266,1320,0),  True,                   beadsize,       23,        [100,18,25],     1,               None,            None, [2]),
         ("Nup159", "Nup159_8",    1.0,     "../data/protein_fasta.Nup159.txt",  "Nup159",  "BEADS",                            " ",  (1321,1334,0),  True,                   beadsize,       24,        [100,18,26],     0,               None,            None, [2]),
         ("Nup159", "Nup159_9",    1.0,     "../data/protein_fasta.Nup159.txt",  "Nup159",  "IDEAL_HELIX",                      " ",  (1335,1369,0),  True,                   beadsize,       25,        [100,18,27],     1,               None,            None, [2]),
         ("Nup159", "Nup159_10",   1.0,     "../data/protein_fasta.Nup159.txt",  "Nup159",  "BEADS",                            " ",  (1370,1379,0),  True,                   beadsize,       26,        [100,18,28],     0,               None,            None, [2]),
         ("Nup159", "Nup159_11",   1.0,     "../data/protein_fasta.Nup159.txt",  "Nup159",  "IDEAL_HELIX",                      " ",  (1380,1410,0),  True,                   beadsize,       27,        [100,18,29],     1,               None,            None, [2]),
         ("Nup159", "Nup159_12",   1.0,     "../data/protein_fasta.Nup159.txt",  "Nup159",  "BEADS",                            " ",  (1411,1432,0),  True,                   beadsize,       28,        [100,18,30],     0,               None,            None, [2]),
         ("Nup159", "Nup159_13",   1.0,     "../data/protein_fasta.Nup159.txt",  "Nup159",  datadirectory+"/3PBP.pdb",          "C",  (1433,1460,0),  True,                   beadsize,       0,         [100,18,31,0],   1,               None,            None, [2]),
         ("Dyn2",   "Dyn2",        0.25,    "../data/protein_fasta.Dyn2.txt",    "Dyn2",    datadirectory+'/4DS1.pdb',          "A",  (1,-1,0),       None,                   beadsize,       19,        [100,18,32],     0,               None,            None, [3])]


bm1=IMP.pmi.macros.BuildModel1(simo)
bm1.set_gmm_models_directory("../data/em_gmm_model/")
bm1.build_model(domains)
resdensities=bm1.get_density_hierarchies([t[1] for t in domains])


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


# tail module em density

mass=sum((IMP.atom.Mass(p).get_mass() for h in resdensities for p in IMP.atom.get_leaves(h)))
gem = IMP.pmi.restraints.em.GaussianEMRestraint(resdensities,'../data/em_density/nup82_lp_090514.mrc.gmm.20.txt',
                                               target_mass_scale=mass,
                                                slope=0.000001,
                                                target_radii_scale=3.0)
gem.add_to_model()
gem.set_weight(100.0)
#gem.center_model_on_target_density(simo)
outputobjects.append(gem)


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



