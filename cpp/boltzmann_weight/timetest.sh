#! /bin/bash

for i in 10 100 1000 ; do
    (time ./Boltzmann_Weight -b -s -n $i -o $1.$i.out) &> $1.$i.time 
    done 
