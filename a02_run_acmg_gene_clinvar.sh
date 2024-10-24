#!/bin/bash
#$ -cwd
#$ -N extract_clinical
#$ -l h_vmem=10G
#$ -l h_rt=00:15:00

. /etc/profile.d/modules.sh
module load igmm/apps/R/4.0.2


line_num=${SGE_TASK_ID}

run_folder=$1

cd ${run_folder}

# here we are opening the gene_list file and taking gene based on the number of row in the file (and the number of row comes from line_num
gene_list=data/acmg_v3.1_gene_list.txt
gene=$(head -n ${line_num} ${gene_list} | tail -n 1)

# now when we have the name of the gene, we run the R script that does everything that's defined in the acmg_gene_clinvar_freqs.r script,
#  but instead of KCNH2 or BRCA1 using the current gene as the gene of interest
R --no-save -f scripts/a03_acmg_gene_clinvar.r --args ${gene} ${run_folder}

mv scripts/extract_clinical.* logs/
