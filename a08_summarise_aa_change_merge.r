#### summarising the number of carriers for each variant and checking if all variants were extracted
#### script written lklaric
#### script started 17032023 (based on summarise_carriers.r)
#### script adjusted to work with carriers of the variants that were merged based on amino acid change
rm(list = ls(all = TRUE))
library(data.table)
library(plyr)
library(dplyr)
library(tidyr)

setwd(run_path)

# merging all carriers into one table
gt_files_hets = list.files("results/genotypes/aa_change_merge/", pattern = "_hets")
gt_files_hom = list.files("results/genotypes/aa_change_merge/", pattern = "_hom")

### heterozygotes
hets = lapply(gt_files_hets, function(f) {
	variant = gsub("variant_.*_(\\d+:\\d+.*)_hets.txt", "\\1", f)
	gene = gsub("variant_(.*)_\\d+:\\d+.*_hets.txt", "\\1", f)
	d = fread(sprintf("results/genotypes/aa_change_merge/%s", f))
#	d = d %>% select(id, gt)
	d = d %>% mutate(variant = variant, gene = gene)
})
hets = do.call(rbind, hets)


### homozygotes

homs = lapply(gt_files_hom, function(f) {
	variant = gsub("variant_.*_(\\d+:\\d+.*)_homo.txt", "\\1", f)
	gene = gsub("variant_(.*)_\\d+:\\d+.*_homo.txt", "\\1", f)
	d = fread(sprintf("results/genotypes/aa_change_merge/%s", f))
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

## merge with pathogenicity, stars etc
all_variants_files = list.files("results/variants/aa_change_merge", pattern = "clinvar_full_aa_change.txt")
all_variants = lapply(all_variants_files, function(f) {
	d = fread(paste0("results/variants/aa_change_merge/", f))
})
all_variants = do.call(rbind, all_variants)

names(all_variants)[grep("AlleleID", names(all_variants))] = "AlleleID"


# merging by exome position because that's in our data; clinvar position/alleles expectedly won't match 
# because that's the whole point of mergeing by aa change
carriers = carriers %>% rename(Position.exome = Position)
carriers = carriers %>% rename(REF.exome = REF)
carriers = carriers %>% rename(ALT.exome = ALT)
all_variants = all_variants %>% mutate(Position.exome = as.character(Position.exome))
all_variants = all_variants %>% mutate(Chromosome = as.character(Chromosome))
carriers = carriers %>% mutate(Position.exome = as.character(Position.exome))
carriers = carriers %>% mutate(Chromosome = as.character(Chromosome))

carriers = left_join(carriers, 
				select(all_variants, Chromosome, Position.exome, Position.clinvar, Name, REF.exome, REF.clinvar, ALT.exome, ALT.clinvar, 
						aa_change, rsid, Type, funct_conseq, isLof,  
						ClinicalSignificance, clinical_significance_categories, ReviewStatus, NumberSubmitters, star_status), 
				by = c("Chromosome", "Position.exome", "REF.exome", "ALT.exome"))

carriers = carriers %>% distinct()


# keeping only those that have non-missing star status: since my extraction of genotypes depended soley on chromosome and position some extracted
# variants would not have necessarily been in the clinvar - this happens for multiallelic variants. these can be removed

carriers = carriers %>% filter(!is.na(star_status))

carriers_variants = carriers %>% select(gene, variant, aa_change, rsid, Chromosome, Position.exome, Position.clinvar, Name, 
										REF.exome, ALT.exome, REF.clinvar, ALT.clinvar, 
										Type, funct_conseq, isLof, 
										ClinicalSignificance, clinical_significance_categories, ReviewStatus, NumberSubmitters, star_status) %>%
					distinct()
carriers_variants = carriers_variants %>% mutate(variant = gsub("\\:aa\\:change\\:merge", "", variant))

n_carriers = carriers %>% 
	group_by(variant, Chromosome, Position.exome, REF.exome, ALT.exome, gt) %>% 
	select(id, variant, Chromosome, Position.exome, REF.exome, ALT.exome, gt) %>%  
	distinct() %>%
	summarise(n = n()) %>% data.frame()
n_carriers = n_carriers %>% mutate(type = ifelse(gt == "0/1", "n_hets", "n_homs"))
n_carriers = n_carriers %>%  select(variant, n, type) %>% spread(type, n) 
n_carriers = n_carriers %>% mutate(variant = gsub("\\:aa\\:change\\:merge", "", variant))

n_carriers_orca = carriers %>% filter(grepl(orca_id, id)) %>% 
	group_by(variant, Chromosome, Position.exome, REF.exome, ALT.exome, gt) %>% 
	select(id, variant, Chromosome, Position.exome, REF.exome, ALT.exome, gt) %>%  
	distinct() %>% 
	summarise(n = n()) %>% data.frame()
n_carriers_orca = n_carriers_orca %>% mutate(type = ifelse(gt == "0/1", "n_hets_orca", "n_homs_orca"))
n_carriers_orca = n_carriers_orca %>%  select(variant, n, type) %>% spread(type, n) 
n_carriers_orca = n_carriers_orca %>% mutate(variant = gsub("\\:aa\\:change\\:merge", "", variant))

n_carriers_viki = carriers %>% filter(grepl(viki_id, id)) %>% 
	group_by(variant, Chromosome, Position.exome, REF.exome, ALT.exome, gt) %>% 
	select(id, variant, Chromosome, Position.exome, REF.exome, ALT.exome, gt) %>% 
	distinct() %>% 
	summarise(n = n()) %>% data.frame()
n_carriers_viki = n_carriers_viki %>% mutate(type = ifelse(gt == "0/1", "n_hets_viki", "n_homs_viki"))
n_carriers_viki = n_carriers_viki %>%  select(variant, n, type) %>% spread(type, n) 
n_carriers_viki = n_carriers_viki %>% mutate(variant = gsub("\\:aa\\:change\\:merge", "", variant))


carriers_variants = left_join(carriers_variants, n_carriers, by = "variant")
carriers_variants = left_join(carriers_variants, n_carriers_orca, by = "variant")
carriers_variants = left_join(carriers_variants, n_carriers_viki, by = "variant")

### adding nfe frequencies
# this file was created here: /gpfs/igmmfs01/eddie/QTL-Regeneron/users/lklaric/rare_variants/scripts/orca_viki_variants_frequencies.r
freqs = fread("data/orcades_viking1_exomes_frequencies_annotations.txt")
freqs = freqs %>% rename(variant = SNP)

carriers_variants = left_join(carriers_variants, select(freqs, variant, maf.orca:shet_ukb), by = "variant")
carriers_variants = carriers_variants %>% distinct()


carriers_variants = carriers_variants %>% arrange(gene, desc(star_status))
# the rsids here are likely not to be trusted because they're coming from clinvar entry, meaning they're not exactly the same as our variants
carriers_variants = carriers_variants %>% select(-rsid)
carriers_variants = carriers_variants %>% distinct()

write.table(carriers_variants, 
	"results/acmg_v3.1_ni_variants_aa_change_merge.txt", 
	sep = "\t", dec = ".", quote = FALSE, row.names = FALSE)


carriers_out = carriers %>% select(gene, variant, gt, aa_change, Name, star_status, ClinicalSignificance, funct_conseq, isLof, id) %>% distinct()

carriers_out = carriers_out %>% arrange(gene, desc(star_status))

write.table(carriers_out, 
	"results/acmg_v3.1_ni_carriers_aa_change_merge.txt", 
	sep = "\t", dec = ".", quote = FALSE, row.names = FALSE)
