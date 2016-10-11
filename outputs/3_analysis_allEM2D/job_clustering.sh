#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -r n
#$ -j y
#$ -l mem_free=4G
#$ -l arch=linux-x64
##$ -l netapp=5G,scratch=5G
##$ -l netappsali=5G
##$ -l scrapp=500G
##$ -l scrapp2=500G
#$ -l h_rt=335:59:59
#$ -R y
#$ -V
##$ -q lab.q
##$ -l hostname="ad*"                #-- anything based on Intel cores
#$ -l hostname="!opt*"			    #-- anything but opt*
##$ -m e                            #-- uncomment to get email when the job finishes
##$ -pe ompi 64
##$ -t 1
#$ -t 1-4                           #-- specify the number of tasks
#$ -N cluster1-4
#########################################
  
# load MPI modules
#module load openmpi-1.6-nodlopen
#module load sali-libraries
#mpirun -V

export IMP=setup_environment.sh
#export PYTHONPATH=$PYTHONPATH:/netapp/sali/kimsj/bin/scikit-learn/lib64/python2.6/site-packages/:/netapp/sali/etjioe/IMP/local_modules/lib64/python/:/netapp/sali/etjioe/IMP/local_modules/lib/python/

echo "NSLOTS = $NSLOTS"
echo "JOB_ID = $JOB_ID"
echo "SGE_TASK_ID = $SGE_TASK_ID"

if [ -z $1 ]; then
    PREFILTER="815.0"
else
    PREFILTER="$1"
fi
echo "PREFILTER = $PREFILTER"

if [ -z $2 ]; then
    NMODS="500"
else
    NMODS="$2"
fi
echo "NMODS = $NMODS"


# write hostname and starting time 
hostname
date

#run
echo "$IMP python ./clustering.py -mpi False -preload True -nmods $NMODS -nclusters $SGE_TASK_ID -prefilter $PREFILTER"
$IMP python ./clustering.py -mpi False -preload True -nmods $NMODS -nclusters $SGE_TASK_ID -prefilter $PREFILTER

# done
hostname
date

