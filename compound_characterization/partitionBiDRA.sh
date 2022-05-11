#!/bin/bash

# Dataset to analyze
dt=$1

# Total number of Batch
B=20
# Nb of batch per run
b=$B
# Run we are currently at
R=1
# Starting and finishing batch numbers given run
FB=$b*$R
SB=$FB-$b


for ((i=$SB; i<$FB; i++)); 
do
    nohup julia bidra.jl $B $i $dt &
done

