for ((i=0; i<=245; i++)); do
   printf -v fn '%05d' $i
   sed "s/opennpdb_00000/opennpdb_${fn}/g" ../open_np_db_data/bash_opennpdb_template.sh > "../open_np_db_data/bash_files/bash_opennpdb_${fn}.sh"
done
