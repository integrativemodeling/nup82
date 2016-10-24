#!/bin/bash

#This scripts requires Gnuplot. We load ImageMagick to convert to PDF but this is not necessary. Just comment out the convert command below.
module load ImageMagick

###Define the two statistics file: 1) histogram 2) statistics for next steps
histogram1D() {
    NPoints=$(wc -l $1 | awk '{print $1**(-1/3)}')
    ESigma=$(awk '{sum+=$4; sum2+=$4*$4; b++}END{print sqrt(sum2/b-(sum/b)^2)}' $1)
    Average=$(awk '{sum+=$4;; b++}END{print sum/b}' $1)
    MaxE=$(awk 'BEGIN{max=0}{if($4>max) max=$4}END{print max}' $1)
    MinE=$(awk 'BEGIN{min=100000000}{if(min > $4) min=$4}END{print min}' $1)
    BinSize=$(awk -v sigm="$ESigma" -v Np="$NPoints" 'BEGIN{print 3.49*sigm*Np}')
    NumberOfBin=$(awk -v max="$MaxE" -v min="$MinE" -v Nb="$BinSize" 'BEGIN{print (max-min)/Nb}')
    awk -v bine="$BinSize" -v bins="$6" '{print $4/(bine*bins)}' $1 | perl -nl -MPOSIX -a -n -e 'print floor($_);' | awk -v bine="$BinSize" -v bins="$6" '{print $1*(bine/bins)}' | sort -n -k1 | uniq -c | awk '{print $2, $1}' | sort -n -k1 > $2

    echo $3 $Average $ESigma $MinE $MaxE $BinSize $NPoints $5 >> $4
}

###Make sure to start from scratch.
if [ -f Scores.dat ]; 
then 
    rm Scores.dat
fi

if [ -f Statistics.dat ];
then
    rm Statistics.dat
fi


if [ -f Statistics-2.dat ];
then
    rm Statistics-2.dat
fi

##Create the file with all the scores
for j in `seq 0 22`;
do
    cd ./Res-Class-${j}
    for i in `seq 1 200`; 
    do 
	grep "^/scrapp/" Scoring.o*.${i} | awk -v RUN="${j}" 'BEGIN{FS=OFS="/"}{printf("%i %s\n", RUN, $5)}' | awk 'BEGIN{FS=OFS="."}{printf("%s %s.%s\n", $1, $2, $3)}' | awk '{print $1, $2, $4, $6}' | awk 'BEGIN{FS=OFS="_"}{printf("%s %s\n", $1, $2)}' | awk '{print $1, $3, $4, $5}' >> ../Scores.dat; 
    done
    cd ..
    echo "Done with Class $j"
done

###Create a binned Score file to plot an histogram.
histogram1D Scores.dat Scores-H.dat 26 Statistics.dat Total 1

MinX=$(cat Scores-H.dat | awk '{print $1}' | sort -r | tail -1 | awk '{printf("%.1f\n", int(10*$1)/10)}')
MaxY=$(cat Scores-H.dat | awk '{print $2}' | sort -n | tail -1 | awk '{printf("%i\n", int(10*(($1/10000)+1)*1000}')

gnuplot <<EOF
reset
set terminal postscript enhanced color "Helvetica, 24"
set xr[0.5:1]
set yr[0:$MaxY]
set xtics 0.5, 0.1, 1
set ytics 0, 1000, $MaxY
set xlabel "Cross Correlation Coefficient"
set ylabel "Number of Images"
set key top left
plot "Scores-H.dat" using 1:2 w p pt 7 ps 1.2 lc -1 title 'Total'
set output "Distribution-Class-Total.eps"
replot   
EOF

###Create a binned Score file to plot an histogram for each class averages scores.
for j in `seq 0 22`;
do
    awk -v Class="${j}" '{if($3 == Class) print $0}' Scores.dat > Scores-C-${j}.dat;
    histogram1D Scores-C-${j}.dat Scores-C-${j}-H.dat ${j} Statistics.dat ${j} 3
done

echo "Done with Statistics for all the classes"

MinTotalX=$(cat Scores-C-*-H.dat | awk '{print $1}' | sort -r | tail -1 | awk '{printf("%.1f\n", int(10*$1)/10)}')
MaxTotalY=$(cat Scores-C-*-H.dat | awk '{print $2}' | sort -n | tail -1)

for j in `seq 0 22`;
do
    gnuplot <<EOF
reset
set terminal postscript enhanced color "Helvetica, 24"
set xr[$MinTotalX:1]
set yr[0:$MaxTotalY]
set xtics $MinTotalX, 0.1, 1
set ytics 0, 250, $MaxTotalY
set xlabel "Cross Correlation Coefficient"
set ylabel "Number of Images"
set key top right
set key box linestyle -1 lw 2.5
plot "Scores-C-${j}-H.dat" using 1:2 w p pt 7 ps 1.2 lc -1 title "Class ${j}"
set output "Distribution-Class-${j}.eps"
replot
EOF
done;

sort -n Statistics.dat > SStatistics.dat

gnuplot <<EOF
set terminal postscript enhanced color "Arial Narrow, 24"
set xlabel "Class Average" offset graph 0,0.038 font ",24"
set ylabel "Average" offset graph 0.06,0 font ",24"
set boxwidth 0.8
set style fill solid 1.00 border lt -1
set xtics border in scale 0,0 mirror norotate offset graph 0, 0.03, 0 autojustify font ",18"
set ytics border in scale 0,0 mirror norotate offset graph 0.015, 0, 0 autojustify font ",24"
set xrange [ -1.000000 : 27.00000 ]
set yrange [ 0.00000 : 1.000 ]
unset key
plot 'SStatistics.dat' using 1:2:xtic(8) with boxes lc -1 notitle, \
'' using 1:2:4:5 with yerrorbars lc rgb 'red' pt 7 ps 0 notitle 

set output 'Average.eps'
replot
EOF
  

mkdir Plot-Distribution
mv Distribution-Class-*.eps Plot-Distribution
mv Average.eps Plot-Distribution

cd Plot-Distribution

for i in `seq 0 22`;
do
    convert -density 300 -compress none -rotate 90 Distribution-Class-${i}.eps Distribution-Class-${i}.pdf
done

convert -density 300 -compress none -rotate 90 Distribution-Class-Total.eps Distribution-Class-Total.pdf
convert -density 300 -compress none -rotate 90 Average.eps Average.pdf
cd ..
