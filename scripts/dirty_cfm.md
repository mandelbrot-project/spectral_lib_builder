# cfm_id

```
mkdir cfm-4
```
```
module load GCC/6.3.0-2.27 Singularity/2.4.2
```
```
singularity build cfm-4/cfm.sif docker://wishartlab/cfmid
```
```
scp Downloads/tmp/lotus/smiles4cfm.txt rutza@login2.baobab.hpc.unige.ch:smiles.txt
```
```
head smiles.txt -n 10 > test.txt
```
```
mkdir test
```
```
split --lines=1 --numeric-suffixes=1 --suffix-length=4 --additional-suffix=.txt test.txt test/test-
```

```
mkdir testout
```

```{run_cfm_test.sh}
#!/bin/sh
#SBATCH --partition=public-cpu
#SBATCH --time=00:30
#SBATCH --mail-user=adriano.rutz@unige.ch
#SBATCH --mail-type=ALL

ml GCC/9.3.0 Singularity/3.7.3-Go-1.14

#convert job array index to three digit padded with zeros
printf -v FILE_INDEX "%04d" ${SLURM_ARRAY_TASK_ID}

FILE=test/test-${FILE_INDEX}.txt

srun singularity run cfm-4/cfm.sif -c "cfm-predict $FILE 0.001 /trained_models_cfmid4.0/[M+H]+/param_output.log /trained_models_cfmid4.0/[M+H]+/param_config.txt 1 testout" 
```
```
sbatch --array=1-10 run_cfm_test.sh
```
```
mkdir smiles
```
```
split --lines=100 --numeric-suffixes=1 --suffix-length=4 --additional-suffix=.txt smiles.txt smiles/smiles-
```

```
mkdir posout
mkdir negout
```

```{run_cfm.sh}
#!/bin/sh
#SBATCH --partition=public-cpu
#SBATCH --time=4-00:00:00
#SBATCH --mail-user=adriano.rutz@unige.ch
#SBATCH --mail-type=ALL

ml GCC/9.3.0 Singularity/3.7.3-Go-1.14

#convert job array index to three digit padded with zeros
printf -v FILE_INDEX "%04d" ${SLURM_ARRAY_TASK_ID}

FILE=smiles/smiles-${FILE_INDEX}.txt

srun singularity run cfm-4/cfm.sif -c "cfm-predict $FILE 0.001 /trained_models_cfmid4.0/[M+H]+/param_output.log /trained_models_cfmid4.0/[M+H]+/param_config.txt 1 posout"
```

```
sbatch --array=1-1461 run_cfm.sh
```

```
sbatch --array=1-1461 run_cfm_neg.sh
```
