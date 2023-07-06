#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=50gb
#SBATCH -t 8:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aron0064@umn.edu
#SBATCH -A mfiecas
#SBATCH -o /panfs/jay/groups/29/mfiecas/aron0064/ActHMM/LogFiles/%A.out
#SBATCH -e /panfs/jay/groups/29/mfiecas/aron0064/ActHMM/LogFiles/%A.err
date
path=/panfs/jay/groups/29/mfiecas/aron0064/ActHMM
cd $path/Routputs
module load R
Rscript $path/Rcode/DataFormat.R