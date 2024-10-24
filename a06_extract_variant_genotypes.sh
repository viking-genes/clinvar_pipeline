#!/bin/bash
#$ -N regen_extract_var
#$ -l h_vmem=2G
#$ -l h_rt=00:05:00

. /etc/profile.d/modules.sh
module load igmm/apps/R/4.0.2
module load igmm/apps/bcftools/1.9

run_folder=$1

cd ${run_folder}


# extracting which variant to look for
line_num=${SGE_TASK_ID}
variant_list=scripts/variants_to_extract.txt
variant_input=$(head -n ${line_num} ${variant_list} | tail -n 1)


chr=$(echo ${variant_input} | awk '{print $1}')
variant_pos=$(echo ${variant_input} | awk '{print $2}')
gene=$(echo ${variant_input} | awk '{print $3}')
variant="${chr}:${variant_pos}"
file_name="variant_${gene}_${variant}"

folder="${gene}_${chr}_${variant_pos}"

# make the variant-specific-folder
mkdir results/genotypes/chr_pos_allele_merge/${folder}

echo "extracting data for "${variant_input}
echo "saving results in the "${folder}" folder"
# extract genotypes from the exome file and save in the variant-specific-folder
exome_path=/exports/igmm/eddie/qtl-viking-genetic-data/orcades/whole_exome/regeneron/v2/raw_data/data/pVCF/GL_by_chrom/
bcftools view -r ${variant} ${exome_path}EDINBURGH_Freeze_Two_remake.${chr}.GL.vcf.gz > results/genotypes/chr_pos_allele_merge/${folder}/${file_name}.vcf

echo "command run: bcftools view -r "${variant}" /exports/igmm/eddie/qtl-viking-genetic-data/orcades/whole_exome/regeneron/v2/raw_data/data/pVCF/GL_by_chrom/EDINBURGH_Freeze_Two_remake."${chr}".GL.vcf.gz > results/genotypes/chr_pos_allele_merge/"${folder}"/"${file_name}".vcf"

# reformat vcf to .txt
grep -v "##" results/genotypes/chr_pos_allele_merge/${folder}/${file_name}.vcf > results/genotypes/chr_pos_allele_merge/${folder}/${file_name}.txt

echo "reformatting vcf to txt"
echo "grep -v "##" results/genotypes/chr_pos_allele_merge/"${folder}"/"${file_name}".vcf > results/genotypes/chr_pos_allele_merge/"${folder}/${file_name}".txt"

# running R script to extract hetero and homozygotes
echo "extracting hetero and homozygotes"
echo "R --no-save -f scripts/a07_extract_carriers.r --args "${folder}" "${file_name}

R --no-save -f scripts/a07_extract_carriers.r --args ${folder} ${file_name} ${run_folder}

mv ~/regen_extract_var.* logs/
