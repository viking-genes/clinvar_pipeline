#!/bin/bash
#$ -cwd
#$ -N clinvar_extract
#$ -l h_vmem=5G
#$ -l h_rt=00:10:00

. /etc/profile.d/modules.sh

line_num=${SGE_TASK_ID}

run_folder=$1

cd ${run_folder}


gene_list=data/gene_list.txt

gene=$(head -n ${line_num} ${gene_list} | tail -n 1)

echo "extracting "${gene}" from data/clinvar/variants_summary.txt.gz file using zgrep  \n"

zgrep ${gene} data/clinvar/variant_summary.txt.gz > data/clinvar/per_gene_clinvar/${gene}_clinvar.txt


header=$(zcat data/clinvar/variant_summary.txt.gz | head -n1)

echo "adding header to the file "
echo ${header}
sed -i "1i${header}" data/clinvar/per_gene_clinvar/${gene}_clinvar.txt

mv scripts/clinvar_extract.* logs/


