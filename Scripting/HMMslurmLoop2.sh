#!/bin/bash -l

TimeArrayNorm=("2" "5" "6" "10" "18" "22")
TimeArrayGamma=("3" "5" "9" "16" "24" "30")
TimeArrayStudent=("2" "5" "6" "10" "18" "22")
MemArray=("5" "10" "10" "20")
GDIST="gamma"
SDIST="student"
for distribution in norm gamma student; do
    for size in 2; do
        for RE_num in 1 2 3 4 5 6; do
			if [ "$distribution" = "$GDIST" ]
			then
				sbatch --time="${TimeArrayGamma[$RE_num -1]}":00:00 --mem="${MemArray[$size - 1]}"gb SimAHMMmulti.sh $RE_num $size  $distribution
			elif [ "$distribution" = "$SDIST" ]
			then
				sbatch --time="${TimeArrayStudent[$RE_num -1]}":00:00 --mem="${MemArray[$size - 1]}"gb SimAHMMmulti.sh $RE_num $size  $distribution
			else 
				sbatch --time="${TimeArrayNorm[$RE_num -1]}":00:00 --mem="${MemArray[$size - 1]}"gb SimAHMMmulti.sh $RE_num $size  $distribution
			fi
        done
    done
done

