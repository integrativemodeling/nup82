These scripts demonstrate the use of [IMP](https://integrativemodeling.org), [MODELLER](https://salilab.org/modeller), and [PMI](https://github.com/salilab/pmi) in the modeling of the Nup82 complex using DSS/EDC chemical cross-links and electron microscopy (EM) 2D class averages.

First, [MODELLER](https://salilab.org/modeller) is used to generate
initial structures for the individual components in the Nup82 complex. Then, IMP
is used to model these components using DSS/EDC crosslinks and the electron microscopy 2D class averages for the entire Nup82 complex.

The modeling protocol will work with a default build of IMP, but for most effective sampling, IMP should be built with [MPI](https://integrativemodeling.org/2.5.0/doc/ref/namespaceIMP_1_1mpi.html) so that replica exchange can be used.

## List of files and descriptions:
1_run_initial_random_EMclass2.sh: script that runs the PMI script 1_modeling_initial_random.py. 

EM2D-Filter.py: Script to compute a score for a given model against a set of class averages. 

Get-Distribution-Statistics.sh: Script to analyze the results of the filtering; requires GNUPlot and ImageMagick. 

Satisfied-Classes-Count.sh: Quick script to format the results of the analysis into a summary text file. 

Models_Selection_Threshold.py: Script to select models given a threshold (either a user input or a adaptive threshold).

2_run_refinement.sh: script that runs the PMI script 2_modeling_allEM_except11_19.py. 

## (1) Running the IMP/PMI scripts for the initial stage of modeling
1) python 1_modeling_initial_random.py -r repeats -out outputdir -em2d class_average -weight em2d_weight

The script will generate independent and uncorrelated trajectories from different random initial configuration restrained by all data and information available for the system (e.g., chemical cross-linking data, excluded volume, structures) and restrained by one EM 2D class average.

## (2) Running the IMP scripts for the filtering stage of modeling:
First, we use IMP to register models against EM class averages and then select models that satisfy those class averages the best. Here, the filtering is applied to the Nup82 Complex. The filtering protocol will work with a default build of IMP. The scripts require GNUplot and ImageMagick.

1) python EM2D-Filter.py input_rmf_file list_of_class_averages angstrom_per_pixel number_of_projections model_resolution image_resolution frame_of_rmf_to_read 

Inputs: input_rmf_file: structure to be projected and registered against class averages 
list_of_class_averages: text file listing the list of class averages.
Output: A list of score for the structure against each image (1-cross correlation).

2) bash Get-Distribution-Satistics.sh 
Inputs: List of files containing the scores. This is hard coded in the script here.
Outputs: Histogram of the score for each images given a set of models Automatically made plots for the distributions of the scores Satistics (average, standard deviation, min, max) for the score of a set of models given a class average.

## (3) Running the IMP/PMI scripts for the refinement stage of modeling
1) python 2_modeling_allEM_except11_19.py -r repeats -out outputdir -rmf starting_rmf -rmf_n rmf_frame -em2d class_average -weight em2d_weight

The script will generate independent and uncorrelated trajectories, refined from the best-scoring configurations after the EM filter.  All data and information available for the system, including chemical cross-linking data, excluded volume, atomic structures, and 21 good EM 2D class averages (except classes 11 and 19) are used.

## Information
_Author(s)_: Seung Joong Kim, Ilan E. Chemmama, Riccardo Pellarin 

_Date_: October 6th, 2016

_License_: [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)
This work is licensed under the Creative Commons Attribution-ShareAlike 4.0
International License.

_Last known good IMP version_: [![build info](https://integrativemodeling.org/systems/?sysstat=23&branch=master)](https://integrativemodeling.org/systems/) [![build info](https://integrativemodeling.org/systems/?sysstat=23&branch=develop)](https://integrativemodeling.org/systems/)

_Testable_: Yes.

_Parallelizeable_: Yes

_Publications_:
 
 \*These authors contributed equally to this work as co-first authors.
