#!/bin/bash

dt=0.0000001
declare -a arr=("Cython" "Interpreted" "Compiled")
for f in "${arr[@]}"
do

  printf "\n\nRunning $f Version: dt=$dt\n--------------------------------"
  time python Test_$f.py $dt

done
