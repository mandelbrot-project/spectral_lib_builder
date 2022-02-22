for ((i=0; i<=429; i++)); do
   printf -v fn '%05d' $i
   test
   sed "s/lotus_00000/lotus_${fn}/g" ../lotus_data/bash_lotus_template.sh > "../lotus_data/bash_files/bash_lotus_${fn}.sh" 
done