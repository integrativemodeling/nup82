NO_MODELS=500
#rm -rf B250
#mkdir B250

#for ((i=0; i<250; i++)); do
for ((i=0; i<10; i++)); do
    #echo $RANDOM % $NO_MODELS | bc
    r=$(( $RANDOM % $NO_MODELS )); echo $r

    mv A250/$r.rmf3 B250
done


