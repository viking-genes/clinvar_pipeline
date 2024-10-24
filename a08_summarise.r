#### summarising the number of carriers for each variant and checking if all variants were extracted
#### script written lklaric
#### script started 17032023
rm(list = ls(all = TRUE))
library(data.table)
library(plyr)
library(dplyr)
library(tidyr)

setwd(run_path)

# merging all carriers into one table
gt_files_hets = list.files("results/genotypes/chr_pos_allele_merge/", pattern = "_hets")
gt_files_hom = list.files("results/genotypes/chr_pos_allele_merge/", pattern = "_hom")

### heterozygotes
hets = lapply(gt_files_hets, function(f) {
	variant = gsub("variant_.*_(\\d+:\\d+.*)_hets.txt", "\\1", f)
	gene = gsub("variant_(.*)_\\d+:\\d+.*_hets.txt", "\\1", f)
	d = fread(sprintf("results/genotypes/chr_pos_allele_merge/%s", f))
#	d = d %>% select(id, gt)
	d = d %>% mutate(variant = variant, gene = gene)
})
hets = do.call(rbind, hets)

### homozygotes

homs = lapply(gt_files_hom, function(f) {
	variant = gsub("variant_.*_(\\d+:\\d+.*)_homo.txt", "\\1", f)
	gene = gsub("variant_(.*)_\\d+:\\d+.*_homo.txt", "\\1", f)
	d = fread(sprintf("results/genotypes/chr_pos_allele_merge/%s", f))
	#d = d %>% select(id, gt)
	d = d %>% mutate(variant = variant, gene = gene)
})
homs = do.call(rbind, homs)

### merging into one table
carriers = rbind(hets, homs)

### add stars to variants
## need to extract chromosome, position, ref and alt for variants
carriers = carriers %>% mutate(variant = gsub("_", ":", variant))
carriers = carriers %>% 
	mutate(Chromosome = sapply(strsplit(variant, split = ":"), function(d) d[1])) %>%
	mutate(Position = sapply(strsplit(variant, split = ":"), function(d) d[2])) 

# renaming variant to match ref and alt observed in the data
carriers = carriers %>% mutate(variant_observed = paste(Chromosome, Position, REF, ALT, sep = ":"))

# variant_exome - the name of the variant in the exome gt files
carriers = carriers %>% rename(variant_exome = variant)

## merge with pathogenicity, stars etc
all_variants_files = list.files("results/variants/chr_pos_allele_merge/", pattern = "clin_signif.txt")
all_variants = lapply(all_variants_files, function(f) {
	d = fread(paste0("results/variants/chr_pos_allele_merge/", f))
})
all_variants = do.call(rbind, all_variants)

names(all_variants)[grep("AlleleID", names(all_variants))] = "AlleleID"

all_variants = all_variants %>% mutate(Position = as.character(Position))
carriers = carriers %>% mutate(Position = as.character(Position))

carriers = left_join(carriers, 
				select(all_variants, Chromosome, Position, Name, REF, ALT, #AlleleID, ClinSigSimple, 
						rsid, Type, funct_conseq, isLof, hgvs.c, hgvs.p, 
						ClinicalSignificance, ClinSigSimple, clinical_significance_categories, ReviewStatus, NumberSubmitters, star_status), 
				by = c("Chromosome", "Position", "REF", "ALT"))

# keeping only those that have non-missing star status: since my extraction of genotypes depended soley on chromosome and position some extracted
# variants would not have necessarily been in the clinvar - this happens for multiallelic variants. these can be removed

carriers = carriers %>% filter(!is.na(star_status))

carriers_variants = carriers %>% select(gene, variant_exome, variant_observed, Chromosome, Position, Name, REF, ALT, 
						rsid, Type, funct_conseq, isLof, hgvs.c, hgvs.p, #AlleleID, 
						ClinicalSignificance, ClinSigSimple, clinical_significance_categories, ReviewStatus, NumberSubmitters, star_status) %>%
					distinct()

n_carriers = carriers %>% group_by(variant_observed, Chromosome, Position, REF, ALT, gt) %>% summarise(n = n()) %>% data.frame()
n_carriers = n_carriers %>% mutate(type = ifelse(gt == "0/1", "n_hets", "n_homs"))
n_carriers = n_carriers %>%  select(variant_observed, n, type) %>% spread(type, n) 

n_carriers_orca = carriers %>% filter(grepl(orca_id, id)) %>% group_by(variant_observed, Chromosome, Position, REF, ALT, gt) %>% summarise(n = n()) %>% data.frame()
n_carriers_orca = n_carriers_orca %>% mutate(type = ifelse(gt == "0/1", "n_hets_orca", "n_homs_orca"))
n_carriers_orca = n_carriers_orca %>%  select(variant_observed, n, type) %>% spread(type, n) 

n_carriers_viki = carriers %>% filter(grepl(viki_id, id)) %>% group_by(variant_observed, Chromosome, Position, REF, ALT, gt) %>% summarise(n = n()) %>% data.frame()
n_carriers_viki = n_carriers_viki %>% mutate(type = ifelse(gt == "0/1", "n_hets_viki", "n_homs_viki"))
n_carriers_viki = n_carriers_viki %>%  select(variant_observed, n, type) %>% spread(type, n) 


carriers_variants = left_join(carriers_variants, n_carriers, by = "variant_observed")
carriers_variants = left_join(carriers_variants, n_carriers_orca, by = "variant_observed")
carriers_variants = left_join(carriers_variants, n_carriers_viki, by = "variant_observed")

### adding nfe frequencies
# this file was created here: /gpfs/igmmfs01/eddie/QTL-Regeneron/users/lklaric/rare_variants/scripts/orca_viki_variants_frequencies.r
freqs = fread("/exports/igmm/eddie/QTL-Regeneron/users/lklaric/rare_variants/results/orcades_viking1_exomes_frequencies_annotations.txt")
freqs = freqs %>% rename(variant_observed = SNP)

carriers_variants = left_join(carriers_variants, select(freqs, variant_observed, maf.orca:shet_ukb), by = "variant_observed")
carriers_variants = carriers_variants %>% distinct()

### adding hfe manually
hfe_all_var = fread(paste0(run_path, "results/variants/chr_pos_allele_merge/HFE_ni_exomes_clinvar_full.txt"))
hfe_var = hfe_all_var %>% filter(rsid == "rs1800562")
hfe_var = hfe_var %>% rename(variant_observed = SNP)

hfe_freq = freqs %>% filter(rsid == "rs1800562")
hfe_var = left_join(hfe_var, select(hfe_freq, variant_observed, maf.orca:shet_ukb), by = "variant_observed")
hfe_var = hfe_var %>% 
	mutate(n_hets = 214+254 , n_homs = 12) %>% 
	mutate(n_hets_orca = 214, n_homs_orca = 6, n_hets_viki = 254, n_homs_viki = 6) %>% mutate(gene = "HFE")
hfe_var = hfe_var %>% mutate(variant_exome = variant_observed)
hfe_var = hfe_var %>% select("gene", "variant_exome", "variant_observed", "Chromosome", "Position", "Name","REF", "ALT", # "AlleleID", 
							"rsid" , "Type", "funct_conseq", "isLof", "hgvs.c" , "hgvs.p", 
							"ClinicalSignificance", "ClinSigSimple", "clinical_significance_categories", "ReviewStatus", "NumberSubmitters", "star_status", 
							"n_hets", "n_homs", "n_hets_orca" , "n_hets_viki", "n_homs_viki", "maf.orca", "mac.orca", "maf.shet", "mac.shet", 
							"maf.ukb", "mac.ukb", "maf.nfe", "orca_nfe", "orca_ukb", "shet_nfe", "shet_ukb")

carriers_variants = rbind(carriers_variants, hfe_var)

carriers_variants = carriers_variants %>% arrange(gene, desc(star_status))

write.table(carriers_variants, 
	"results/acmg_v3.1_ni_variants_23052023.txt", 
	sep = "\t", dec = ".", quote = FALSE, row.names = FALSE)


carriers_out = carriers %>% select(gene, variant_observed, gt, rsid, Name, star_status, ClinicalSignificance, funct_conseq, isLof, id) %>% distinct()

### adding hfe manually; we're only interested in homozygotes
hfe_carriers = fread("data/variant_HFE_6:26092913_G_A_homo.txt")
hfe_carriers = hfe_carriers %>% mutate(variant_observed = unique(hfe_var$variant_observed))

hfe_carriers = left_join(select(hfe_carriers, variant_observed, gt, id), hfe_var, by = "variant_observed", multiple = "all")
hfe_carriers = hfe_carriers %>% select(gene, variant_observed, gt, rsid, Name, star_status, ClinicalSignificance, funct_conseq, isLof, id) %>% distinct()

carriers_out = rbind(carriers_out, hfe_carriers)
carriers_out = carriers_out %>% arrange(gene, desc(star_status))

write.table(carriers_out, 
	"results/acmg_v3.1_ni_carriers.txt", 
	sep = "\t", dec = ".", quote = FALSE, row.names = FALSE)

