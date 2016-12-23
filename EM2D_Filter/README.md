These scripts demonstrate how to use EM2D as a filter. 
First, we use [IMP](http://integrativemodeling.org) to register models against EM class averages and then select models that satisfy those class averages the best.
Here, the filtering is applied to the Nup82 Complex.
The filtering protocol will work with a default build of IMP. The scripts require GNUplot and ImageMagick.
## List of files:
EM2D-Filter.py: Script to compute a score for a given model against a set of class averages.
Get-Distribution-Statistics.sh: Script to analyze the results of the filtering; requires GNUPlot and ImageMagick.
Satisfied-Classes-Count.sh: Quick script to format the results of the analysis into a summary text file.
Models_Selection_Threshold.py: Script to select models given a threshold (either a user input or a adaptive threshold).

## Running the IMP/PMI scripts for filtering stage of modeling:
1) python EM2D-Filter.py input_rmf_file list_of_class_averages angstrom_per_pixel number_of_projections model_resolution image_resolution frame_of_rmf_to_read
Inputs:   input_rmf_file: structure to be projected and registered against an class averages
          list_of_class_averages: text file listing the list of class averages.

Output:   A list of score for the structure against each image (1-cross correlation). 

2) bash Get-Distribution-Satistics.sh
Inputs:   List of files containing the scores. This is hard coded in the script here. 

Outputs:   Histogram of the score for each images given a set of models 
           Automatically made plots for the distributions of the scores
           Satistics (average, standard deviation, min, max) for the score of a set of models given a class average.

## Information

_Author(s)_: Ilan E. Chemmama and Seung Joong Kim

_Date_: December 21, 2016

_License_: [GPL](http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html).
This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

_Last known good IMP version_: [![build info](https://integrativemodeling.org/systems/?sysstat=6&branch=master)](http://integrativemodeling.org/systems/) [![build info](https://integrativemodeling.org/systems/?sysstat=6&branch=develop)](http://integrativemodeling.org/systems/)

_Testable_: Yes.

_Parallelizeable_: Yes

_Publications_:
