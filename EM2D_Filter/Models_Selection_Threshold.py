#!/usr/bin/env python

from __future__ import print_function
import os
import sys
from math import sqrt

StatFile = sys.argv[1]
threshold_value = sys.argv[2]
stats = []

with open(StatFile, "r") as inn:
    for lines in inn:
        stats.append([float(n) for n in lines.strip('\n').split(' ')])


Threshold = []
for ClassNum in range (0,23):
    Threshold.append(stats[ClassNum][1]+2*stats[ClassNum][2])

for ClassNum in range (0,23):
    
    FileName = "Scores-C-%i.dat" % ClassNum
    classes = []
    
    with open(FileName, "r") as ins:
        print(ClassNum)
        for line in ins:
            classes.append([float(n) for n in line.strip('\n').split(' ')])
            
    for cls in classes:
        if(cls[3] > float(threshold_value)/1000):
            print(int(cls[0]), int(cls[1]), int(cls[2]))

