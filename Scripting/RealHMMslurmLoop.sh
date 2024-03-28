#!/bin/bash -l

TimeArray=("25" "25" "30" "35" "40" "45" "50" "65" "75")
MemArray=("30" "30" "30" "40" "40" "50" "50" "60" "60")
for distribution in norm; do
    for size in 1; do
        for RE_num in 0; do
			sbatch --time="${TimeArray[$RE_num]}":00:00 --mem="${MemArray[$RE_num]}"gb RealHMMmulti.sh $RE_num $size  $distribution
        done
    done
done


