#!/bin/sh
#SBATCH --partition=public-cpu
#SBATCH --time=4-00:00:00
#SBATCH --mail-user=adriano.rutz@unige.ch
#SBATCH --mail-type=ALL

ml GCC/9.3.0 Singularity/3.7.3-Go-1.14

#convert job array index to three digit padded with zeros
printf -v FILE_INDEX "%04d" ${SLURM_ARRAY_TASK_ID}

FILE=smiles/smiles-${FILE_INDEX}.txt

srun singularity run cfm-4/cfm.sif -c "cfm-predict $FILE 0.001 /trained_models_cfmid4.0/[M+H]+/param_output.log /trained_models_cfmid4.0/[M+H]+/param_config.txt 1 posout 0 0"