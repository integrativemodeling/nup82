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
#$ -l hostname="i*"                 #-- anything based on Intel cores
##$ -l hostname="!opt*"			    #-- anything but opt*
##$ -m e                            #-- uncomment to get email when the job finishes
#$ -pe ompi 32
##$ -t 1
#$ -t 1-10                        #-- specify the number of tasks
#$ -N n82_1-10
#########################################

#: '#lyre usage : nohup ./job_test.sh 20000 output > job_test.log &
NSLOTS=4    ## Should be an "EVEN number" or 1
SGE_TASK_ID=2
#'
# load MPI modules
#module load openmpi-1.6-nodlopen
#module load sali-libraries
#mpirun -V

export IMP=setup_environment.sh
MODELING_SCRIPT=1_modeling_initial_random.py
XL_FILE=XL.csv
EM2D_FILE=../data/em2d/2.pgm
EM2D_WEIGHT=10000.0


# Parameters
if [ -z $1 ]; then
    REPEAT="5000"
else
    REPEAT="$1"
fi
echo "number of REPEATs = $REPEAT"

if [ -z $2 ]; then
    OUTPUT="output"
else
    OUTPUT="$2"
fi
echo "OUTPUT foler = $OUTPUT"

echo "SGE_TASK_ID = $SGE_TASK_ID"
echo "JOB_ID = $JOB_ID"
echo "NSLOTS = $NSLOTS"

# write hostname and starting time 
hostname
date

let "SLEEP_TIME=$SGE_TASK_ID*2"
#sleep $SLEEP_TIME

PWD_PARENT=$(pwd)

i=$(expr $SGE_TASK_ID)
DIR=modeling$i

rm -rf $DIR
if [ ! -d $DIR ]; then
    mkdir $DIR
    cp -pr template/$MODELING_SCRIPT $DIR
    cp -pr template/em2d_nup82.py $DIR
fi
cd $DIR

PWD=$(pwd)
echo $PWD_PARENT : $PWD

if [ $PWD_PARENT != $PWD ]; then
    # run the job
    #mpirun -np $NSLOTS $IMP python ./$MODELING_SCRIPT -r $REPEAT -out $OUTPUT -em2d $EM2D_FILE -weight $EM2D_WEIGHT
    $IMP python ./$MODELING_SCRIPT -r $REPEAT -out $OUTPUT -em2d $EM2D_FILE -weight $EM2D_WEIGHT
    cd ..
fi

# done
hostname
date
