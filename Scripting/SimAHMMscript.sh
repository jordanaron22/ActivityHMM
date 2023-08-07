#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --array=1-400
#SBATCH --mem=40gb
#SBATCH -t 5:00:00
#SBATCH --mail-type=NONE
#SBATCH --mail-user=aron0064@umn.edu
#SBATCH -A mfiecas
#SBATCH -o /panfs/jay/groups/29/mfiecas/aron0064/ActHMM/LogFiles/%A_%a.out
#SBATCH -e /panfs/jay/groups/29/mfiecas/aron0064/ActHMM/LogFiles/%A_%a.err
date
path=/panfs/jay/groups/29/mfiecas/aron0064/ActHMM
cd $path/Routputs
module load R
Rscript $path/Rcode/AHMM.R