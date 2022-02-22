for ((i=0; i<=412; i++)); do
   printf -v fn '%04d' $i
   test
   sed "s/coconut_0000/coconut_${fn}/g" ../coconut_data/Bash_Coconut_template.sh > "../coconut_data/bash_files/Bash_Coconut_${fn}.sh" 
done