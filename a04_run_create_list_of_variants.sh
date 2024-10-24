#!/bin/bash
#$ -cwd
#$ -N create_var_list 
#$ -l h_vmem=5G
#$ -l h_rt=00:10:00

. /etc/profile.d/modules.sh
module load igmm/apps/R/4.0.2



run_folder=$1

cd ${run_folder}


R --no-save -f scripts/a05_create_list_of_variants.r --args ${run_folder}

mv scripts/create_var_list.* logs/
