#!/bin/bash -l

TimeArray1=("1" "1" "1" "2" "2" "3" "3" "4" "4")
TimeArray2=("8" "3" "3" "5" "5" "8" "8" "10" "10")
#TimeArray3=("15" "4" "4" "8" "10" "15" "15" "20" "20")
TimeArray3=("3" "3" "3" "4" "4" "4" "5" "6" "6")
TimeArray4=("24" "6" "6" "12" "15" "24" "24" "30" "36")
MemArray=("5" "10" "10" "20")
for distribution in gamma norm; do
    for size in 1 3; do
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

