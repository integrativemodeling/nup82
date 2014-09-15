
mapname=$1
ngaussians=$2
voxelsize=3.23
~/imp-projects/imp-020914/imp-fast-mpich/setup_environment.sh python ~/imp-projects/imp-020914/imp/modules/isd_emxl/bin/create_gmm.py $mapname.mrc $ngaussians $mapname.mrc.gmm.$ngaussians.txt -s 0.00079 -m $mapname.mrc.gmm.$ngaussians.mrc  -a 3.23 -i 200
