#!/bin/bash
#SBATCH --job-name=PPC_M21_2021_08_04_D9_TT4_01_3
#SBATCH --account=def-wilsyl
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32000#SBATCH --mail-user=ecarmichael@gmail.com
#SBATCH --mail-type=ALL
module load StdEnv/2020
module load matlab/2020a

matlab -nodisplay -batch "Beluga_PPC('/lustre04/scratch/ecar/M21/M21_2021-08-04_D9', 'TT4_01_3')
