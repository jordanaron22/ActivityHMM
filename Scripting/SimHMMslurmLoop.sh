#!/bin/bash -l

TimeArray1=("1" "1" "1" "2" "2" "2" "3" "3" "3")
TimeArray2=("3" "3" "3" "4" "4" "6" "6" "6" "6")
TimeArray3=("3" "3" "3" "4" "4" "6" "6" "6" "6")
TimeArray4=("8" "3" "8" "8" "12" "16" "20" "30" "30")
MemArray=("5" "10" "10" "20")
for distribution in stud3nt; do
    for size in 3 4; do
        for RE_num in 0 1 2 3 4 5 6 7 8; do
			if [[ "$size" == "1" ]];then
				sbatch --time="${TimeArray1[$RE_num]}":00:00 --mem="${MemArray[$size - 1]}"gb SimAHMMmulti.sh $RE_num $size  $distribution
			elif [[ "$size" == "2" ]];then
				sbatch --time="${TimeArray2[$RE_num]}":00:00 --mem="${MemArray[$size - 1]}"gb SimAHMMmulti.sh $RE_num $size  $distribution
			elif [[ "$size" == "3" ]];then
				sbatch --time="${TimeArray3[$RE_num]}":00:00 --mem="${MemArray[$size - 1]}"gb SimAHMMmulti.sh $RE_num $size  $distribution
			else 
				sbatch --time="${TimeArray4[$RE_num]}":00:00 --mem="${MemArray[$size - 1]}"gb SimAHMMmulti.sh $RE_num $size  $distribution
			fi
        done
    done
done

