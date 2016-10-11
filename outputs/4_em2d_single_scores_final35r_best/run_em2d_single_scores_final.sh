FILE=225-17r.pdb
#FILE=225-17.pdb

#class_start=0
class_start=23

class_end=0
#class_end=$class_start

#indi_runs_start=6
#indi_runs_end=$indi_runs_start
#indi_runs=3

#n_projection=100
n_projection=1000
n_resolution=35
n_pixel_size=3.23

for ((i=0; i<=10; i++)); do
    PGM+="../data/em2d/${i}.pgm "
done
for ((i=12; i<=18; i++)); do
    PGM+="../data/em2d/${i}.pgm "
done
for ((i=20; i<=22; i++)); do
    PGM+="../data/em2d/${i}.pgm "
done

##############################
## Repeat for class 0 to 22
##############################
#cd em2d_single_scores_final

for ((i=$class_start; i>=$class_end; i--)); do
    echo ""
    echo "######################"
    echo "## EM 2D class $i   ##"
    echo "######################"

    if [ "$i" != "$class_start" ]; then
        PGM=../data/em2d/${i}.pgm
    fi

    #DIR=../modeling1/output/pdbs
    #DIR=../analysis/kmeans_1000_1/cluster.0

    #echo "copying from $DIR/$FILE"
    #cp -pr $DIR/$FILE .
    
    python ./process_for_em.py -pdb $FILE -out temp_$FILE
    
    grep -vwE "(END|ENDMDL|MODEL      0)" temp_$FILE > tr_$FILE
    rm -rf chain*.pdb
    rm -rf temp_$FILE
    
    em2d_single_score tr_$FILE -r $n_resolution -s $n_pixel_size -n $n_projection -l 2 $PGM
    
    #mv best_projections.pgm ${i}_best_projections.pgm
    #mv images.pgm ${i}_images.pgm
    #mv $FILE ${i}_$FILE
    #mv tr_$FILE ${i}_tr_$FILE

    mv best_projections.pgm ${i}_best_projections.pgm
    mv images.pgm no_bg.${i}.pgm
    #mv $FILE ${i}_$FILE
    mv tr_$FILE ${i}_tr_$FILE
done

cd ..

