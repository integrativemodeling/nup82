import IMP
import IMP.pmi
#import IMP.pmi.macros
import sys

#####################################################
# Parsing parameter inputs
#####################################################
import argparse

parser = argparse.ArgumentParser(description='generate clusters of the RMF solutions')
parser.add_argument('-mpi', action="store", dest="is_mpi", help="mpi enabled")
parser.add_argument('-preload', action="store", dest="load_distance_matrix_file", help="skip the matrix calcuklation and read the precalculated matrix")
parser.add_argument('-nmods', action="store", dest="nbestscoringmodels", help="number of models to be clustered")
parser.add_argument('-nclusters', action="store", dest="nclusters", help="number of clusters to be used by kmeans algorithm")
parser.add_argument('-prefilter', action="store", dest="prefiltervalue", help="prefilter the models by score")
inputs = parser.parse_args()

# mpi enabled
if (inputs.is_mpi=="True") or (inputs.is_mpi=="true") or (inputs.is_mpi=="Yes") or (inputs.is_mpi=="yes") :
    inputs.is_mpi = True
else:
    inputs.is_mpi = False

# skip the matrix calcuklation and read the precalculated matrix
if (inputs.load_distance_matrix_file=="True") or (inputs.load_distance_matrix_file=="true") or (inputs.load_distance_matrix_file=="Yes") or (inputs.load_distance_matrix_file=="yes") :
    inputs.load_distance_matrix_file = True
else:
    inputs.load_distance_matrix_file = False

# number of models to be clustered
if inputs.nbestscoringmodels==None:
    inputs.nbestscoringmodels = 500

# number of clusters to be used by kmeans algorithm
if inputs.nclusters==None:
    inputs.nclusters = 1

# prefilter the models by score
if inputs.prefiltervalue==None:
    inputs.prefiltervalue = 565.0

print inputs

is_mpi = inputs.is_mpi                                          # mpi enabled
load_distance_matrix_file = inputs.load_distance_matrix_file    # skip the matrix calcuklation and read the precalculated matrix
nbestscoringmodels = int(inputs.nbestscoringmodels)             # number of models to be clustered
nclusters = int(inputs.nclusters)                               # number of clusters to be used by kmeans algorithm
prefiltervalue = float(inputs.prefiltervalue)                   # prefilter the models by score


#####################################################
# initialize the macro
#####################################################
import macros_n82

model=IMP.Model()

mc=macros_n82.AnalysisReplicaExchange0(model,
                                        stat_file_name_suffix="stat",     # don't change
                                        merge_directories=["../modeling1_failed",
                                                           "../modeling2_failed",
                                                           "../modeling3_failed",
                                                           "../modeling3_failed1",
                                                           "../modeling3_failed2",
                                                           "../modeling4_failed",
                                                           "../modeling5_failed",
                                                           "../modeling6_failed",
                                                           "../modeling7_failed",
                                                           "../modeling8_failed"],
                                        #merge_directories=["../modeling11"], # change this list splitting the runs or adding new runs
                                        global_output_directory="output/")

# fields that have to be extracted for the stat file

feature_list=["ElectronMicroscopy2D_None",
              "ISDCrossLinkMS_Distance_intrarb",
              "ISDCrossLinkMS_Distance_interrb",
              "ISDCrossLinkMS_Data_Score",
              #"GaussianEMRestraint_None",
              "SimplifiedModel_Linker_Score_None",
              "ExcludedVolumeSphere_None",
              "ISDCrossLinkMS_Psi",
              "ISDCrossLinkMS_Sigma"]

# Dictionary of densities to be calculated
# the key is the name of the file and the value if the selection
# example:
#              {"med17-CTD":[(200,300,"med17")],"med17-CTD.med14":[(200,300,"med17"),"med14"]   }


reduced_density_dict={"Nup82.1_NTD":[(1,521,"Nup82.1")],
                      "Nup82.1_CTD":[(522,713,"Nup82.1")],
                      "Nup82.2_NTD":[(1,521,"Nup82.2")],
                      "Nup82.2_CTD":[(522,713,"Nup82.2")],
                      "Nup82_NTD":[(1,521,"Nup82.1"),(1,521,"Nup82.2")],
                      "Dyn2":["Dyn2.1","Dyn2.2"],
                      "Nsp1.1_CTD":[(637,823,"Nsp1.1")],
                      "Nsp1.2_CTD":[(637,823,"Nsp1.2")],
                      "Nup159.1_CTD":[(1117,1460,"Nup159.1")],
                      "Nup159.2_CTD":[(1117,1460,"Nup159.2")],
                      "Subunit1":["Nup82.1", "Dyn2.1", (637,823,"Nsp1.1"), (1117,1460,"Nup159.1")],
                      "Subunit2":["Nup82.2", "Dyn2.2", (637,823,"Nsp1.2"), (1117,1460,"Nup159.2")],
                      "Whole":["Nup82.1","Nup82.2", "Dyn2.1","Dyn2.2", (637,823,"Nsp1.1"),(637,823,"Nsp1.2"), (1117,1460,"Nup159.1"), (1117,1460,"Nup159.2")]}

# list of component names needed to calculate the RMSD for the clustering

components_names={"Nup82.1":"Nup82.1",
                  "Nup82.2":"Nup82.2",
                  "Dyn2.1":"Dyn2.1",
                  "Dyn2.2":"Dyn2.2",
                  "Nsp1.1_CTD":(637,823,"Nsp1.1"),
                  "Nsp1.2_CTD":(637,823,"Nsp1.2"),
                  "Nup159.1_CTD":(1117,1460,"Nup159.1"),
                  "Nup159.2_CTD":(1117,1460,"Nup159.2")}

    
mc.clustering("SimplifiedModel_Total_Score_None",  # don't change, field where to find the score
              "rmf_file",                          # don't change, field where to find the path for the rmf_file
              "rmf_frame_index",                   # don't change, field for the frame index
              prefiltervalue=prefiltervalue,               # prefilter the models by score
              number_of_best_scoring_models=nbestscoringmodels,   # number of models to be clustered
              #alignment_components=None,           # don't change, (list of proteins you want to use for structural alignment
              alignment_components=components_names,         # don't change, (list of proteins you want to use for structural alignment
              rmsd_calculation_components=components_names,  # list of proteins used to calculated the rmsd
              distance_matrix_file="distance.rawmatrix.pkl", # save the distance matrix
              outputdir="kmeans_"+str(nbestscoringmodels)+"_"+str(nclusters)+"/",  # directory name for the clustering
              feature_keys=feature_list,                     # extract these fields from the stat file
              load_distance_matrix_file=load_distance_matrix_file,                # skip the matrix calcuklation and read the precalculated matrix
              skip_clustering=False,                         # skip clustering
              display_plot=True,                            # display the heat map plot of the distance matrix
              exit_after_display=False,                      # exit after having displayed the distance matrix plot
              get_every=1,                                   # skip structures for faster computation
              #is_mpi=is_mpi,                                 # mpi enabled
              number_of_clusters=nclusters,                  # number of clusters to be used by kmeans algorithm
              voxel_size=5.0,                                # voxel size of the mrc files
              density_custom_ranges=reduced_density_dict)    # setup the list of densities to be calculated

