#!/bin/bash
#SBATCH --job-name=PPC_matlab
#SBATCH --account=def-wilsyl
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16000
#SBATCH --mail-user=ecarmichael@gmail.com
#SBATCH --mail-type=ALL

module load StdEnv/2020
module load matlab/2020a

matlab -nodisplay -batch "Beluga_PPC_setup('M23_2021-08-05_D10')"
