#!/bin/sh 

#SBATCH --partition=mono-shared-EL7
#SBATCH --job-name=cfm-predict_npatlas
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8000
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --mail-user=user
#SBATCH --mail-type=FAIL
#SBATCH --clusters=baobab
#SBATCH --output=slurm-%J.out


module load foss/2018b GCC/7.3.0-2.30 OpenMPI/3.1.1 CFM-ID/r33


srun cfm-predict /home/allardp/data/npatlas/npatlas_split_0000002.txt 0.001 /home/allardp/CFM_models/metab_se_cfm/params_metab_se_cfm/param_output0.log /home/allardp/CFM_models/metab_se_cfm/param_config.txt 1 /home/allardp/CFM_results/npatlas/ 0
