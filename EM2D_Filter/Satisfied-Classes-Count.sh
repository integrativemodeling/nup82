#!/bin/bash

if [ -f A-900.dat ]
then
    rm A-*.dat S-*.dat T-*.dat
fi

for j in 875 880 885 890 895 900;
do 
    awk '{print $1, $2, $3}' Statistics.dat > tmpS
    python Models-Selection-2Sigma.py tmpS  ${j} | sort -n -k1 -k2 -k3 | awk '{h[$1" "$2]=h[$1" "$2]" "$3; b[$1" "$2]++}END{for(i in h) print b[i], i,h[i]}' | sort -n -k1 -k2 -k3 > A-${j}.dat


    for i in `seq 4 26`; 
    do 
	awk -v col="${i}" '{if($col=="6") print $0}' A-${j}.dat >> tmp6
	awk -v col="${i}" '{if($col=="13") print $0}' A-${j}.dat >> tmp13
    done

    sort -n -k1 -k2 -k3 tmp6 > S-${j}.dat
    sort -n -k1 -k2 -k3 tmp13 > T-${j}.dat

    cp A-${j}.dat tmpN

    for i in `seq 4 26`; 
    do 
	awk -v col="${i}" '{if($col=="11" || $col=="19") print $0}' A-${j}.dat >> tmpN
    done

    cat tmpN | sort -n -k1 -k2 -k3 -k4 | uniq -c | awk '{if($1=="1") {for (i=2; i<NF; i++) printf $i " "; print $NF}}' > N-${j}.dat

    rm tmp6 tmp13 tmpN tmpS

done