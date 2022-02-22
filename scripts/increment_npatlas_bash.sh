for ((i=0; i<=25; i++)); do
   printf -v fn '%07d' $i
   sed "s/npatlas_split_0000000/npatlas_split_${fn}/g" ../npatlas_data/Bash_npatlas_template.sh > "../npatlas_data/bash_files/Bash_npatlas_${fn}.sh" 
done
