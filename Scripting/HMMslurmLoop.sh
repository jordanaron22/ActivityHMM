#!/bin/bash -l

TimeArray=("20" "40" "60" "80" "96" "96")
MemArray=("20" "20" "30" "30" "40" "40")
for distribution in norm; do
    for size in 1; do
        for RE_num in 2 3 4 5 6; do
			sbatch --time="${TimeArray[$RE_num - 1]}":00:00 --mem="${MemArray[$RE_num - 1]}"gb AHMMmulti.sh $RE_num $size  $distribution
        done
    done
done


