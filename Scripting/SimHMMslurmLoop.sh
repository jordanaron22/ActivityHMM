#!/bin/bash -l

TimeArray1=("1" "1" "2" "2" "3" "3")
TimeArray2=("3" "3" "5" "5" "8" "8")
TimeArray3=("4" "4" "8" "10" "15" "15")
TimeArray4=("6" "6" "12" "15" "24" "24")
MemArray=("5" "10" "10" "20")
for distribution in student3; do
    for size in 3 4; do
        for RE_num in 1 2 3 4 5 6; do
			if [[ "$size" == "1" ]];then
				sbatch --time="${TimeArray1[$RE_num -1]}":00:00 --mem="${MemArray[$size - 1]}"gb SimAHMMmulti.sh $RE_num $size  $distribution
			elif [[ "$size" == "2" ]];then
				sbatch --time="${TimeArray2[$RE_num -1]}":00:00 --mem="${MemArray[$size - 1]}"gb SimAHMMmulti.sh $RE_num $size  $distribution
			elif [[ "$size" == "3" ]];then
				sbatch --time="${TimeArray3[$RE_num -1]}":00:00 --mem="${MemArray[$size - 1]}"gb SimAHMMmulti.sh $RE_num $size  $distribution
			else 
				sbatch --time="${TimeArray4[$RE_num -1]}":00:00 --mem="${MemArray[$size - 1]}"gb SimAHMMmulti.sh $RE_num $size  $distribution
			fi
        done
    done
done

