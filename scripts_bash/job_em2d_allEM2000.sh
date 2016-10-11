#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -r n
#$ -j y
#$ -l mem_free=2G
#$ -l arch=linux-x64
##$ -l netapp=5G,scratch=5G
##$ -l netappsali=5G
##$ -l scrapp=500G
##$ -l scrapp2=500G
#$ -l h_rt=335:59:59
#$ -R y
#$ -V
##$ -q lab.q
##$ -l hostname="i*"                 #-- anything based on Intel cores
#$ -l hostname="!opt*"			    #-- anything but opt*
##$ -m e                            #-- uncomment to get email when the job finishes
#$ -pe ompi 8
##$ -t 1
#$ -t 2001-2020                         #-- specify the number of tasks
#$ -N n82_c2
#########################################

: '#lyre usage : nohup ./job_test.sh 20000 output > job_test.log &
NSLOTS=1    ## Should be an "EVEN number" or 1
SGE_TASK_ID=1
'
# load MPI modules
#module load openmpi-1.6-nodlopen
#module load sali-libraries
#mpirun -V

export IMP=setup_environment.sh
MODELING_SCRIPT=modeling_allEM_except11_19.py
SAXS_FILE=SAXS.dat
XL_FILE=XL.csv
RMF_FILE=../data/modeling2_0.rmf3
#RMF_FILE=../data/rmfs_1st/dr30_0.rmf3
#RMF_FILE=../data/rmfs_1st/2_8583.rmf
RMF_FRAME=0
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
    #cp -pr template/representation_nup82.py $DIR
    #cp -pr template/crosslinking_nup82.py $DIR
fi
cd $DIR

PWD=$(pwd)
echo $PWD_PARENT : $PWD

if [ $PWD_PARENT != $PWD ]; then
    # run the job
    #mpirun -np $NSLOTS $IMP python ./$MODELING_SCRIPT -r $REPEAT -out $OUTPUT -em2d $EM2D_FILE -weight $EM2D_WEIGHT
    mpirun -np $NSLOTS $IMP python ./$MODELING_SCRIPT -r $REPEAT -out $OUTPUT -rmf $RMF_FILE -rmf_n $RMF_FRAME -em2d $EM2D_FILE -weight $EM2D_WEIGHT
    cd ..
fi

# done
hostname
date
#exit -1

################################################################################
#cd $OUTPUT

#process_output.py -f modeling1/output/stat.0.out -s ISDCrossLinkMS_Distance_interrb_13-State:0-127:Cul5_80:EloC-1-1-1.0_DSS ISDCrossLinkMS_Score_interrb_13-State:0-127:Cul5_80:EloC-1-1-1.0_DSS ISDCrossLinkMS_PriorSig_Score_DSS ISDCrossLinkMS_Sigma_1_DSS
#process_output.py -f modeling1/output/stat.1.out -s ISDCrossLinkMS_Distance_interrb_13-State:0-127:Cul5_80:EloC-1-1-1.0_DSS ISDCrossLinkMS_Score_interrb_13-State:0-127:Cul5_80:EloC-1-1-1.0_DSS ISDCrossLinkMS_PriorSig_Score_DSS ISDCrossLinkMS_Sigma_1_DSS

#process_output.py -f modeling1/output/stat.0.out -n 1
#process_output.py -f modeling1/output/stat.0.out -p

################################################################################
### searching for correlationo of the pdb file with rmf file / output log file
################################################################################
#process_output.py -f modeling1/output/stat.0.out -s SimplifiedModel_Total_Score_None ISDCrossLinkMS_Data_Score_scEDC ElectronMicroscopy2D_None Stopwatch_None_delta_seconds rmf_file rmf_frame_index
#process_output.py -f modeling1/output/stat.0.out -s SimplifiedModel_Total_Score_None ISDCrossLinkMS_Data_Score_scEDC ElectronMicroscopy2D_None Stopwatch_None_delta_seconds
#process_output.py -f modeling1/pre-EM2D_output/stat.0.out -s SimplifiedModel_Total_Score_None ISDCrossLinkMS_Data_Score_scEDC Stopwatch_None_delta_seconds
#process_output.py -f modeling1/output/stat.0.out -s SimplifiedModel_Total_Score_None rmf_file rmf_frame_index | grep 103.147


#rmf_slice 29.rmf3 n82_r29_f1445.rmf3 -f 1445 -s 1000000
#em2d_single_score model.0.pdb -r 20 -s 3.23 -n 400 -c *.pgm
