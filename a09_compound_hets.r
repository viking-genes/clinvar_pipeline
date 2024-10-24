rm(list = ls(all = TRUE))
library(data.table)
library(plyr)
library(dplyr)
library(tidyr)

chr_carr = fread("results/acmg_v3.1_clinvar_ni_carriers.txt")
aa_carr = fread("results/acmg_v3.1_clinvar_ni_carriers_aa_change_merge.txt")


chr_carr = chr_carr %>% mutate(type = "chr_pos")
aa_carr = aa_carr %>% mutate(type = "aa")
aa_carr = aa_carr %>% mutate(variant_observed = gsub("\\:aa\\:change\\:merge", "", variant))

### compound hets - within the gene people with multiple variants
gene_n_var = chr_carr %>% group_by(id, gene) %>% select(variant_observed) %>% distinct() %>% summarise(n_var = n())
# one compound het, for MUTYH gene

chr_comp = gene_n_var %>% filter(n_var >= 2)
chr_comp = left_join(chr_comp, chr_carr, by = c("id", "gene"))

write.table(chr_comp, 
	"results/acmg_v3.1_clinvar_ni_compound_hets.txt", 
	sep = "\t", dec = ".", quote = FALSE, row.names = FALSE)

### adding in the amino-acid change variants
chr_carr = chr_carr %>% mutate(aa_change = NA)
aa_carr = aa_carr %>% mutate(rsid = NA)

chr_carr = chr_carr %>% select(id, gene, variant_observed, gt, aa_change, rsid, Name, star_status, ClinicalSignificance, funct_conseq, isLof, type)
aa_carr = aa_carr %>% select(id, gene, variant_observed, gt, aa_change, rsid, Name, star_status, ClinicalSignificance, funct_conseq, isLof, type)

all_carr = rbind(chr_carr, aa_carr)
all_carr = all_carr %>% distinct()

all_gene_n_var = all_carr %>% group_by(id, gene) %>% select(variant_observed) %>% distinct() %>% summarise(n_var = n())


all_comp = all_gene_n_var %>% filter(n_var >= 2)
all_comp = left_join(all_comp, all_carr, by = c("id", "gene"))


write.table(all_comp, 
	"results/acmg_v3.1_clinvar_ni_compound_hets_aa_change_included.txt", 
	sep = "\t", dec = ".", quote = FALSE, row.names = FALSE)


### looking at individuals that have mutliple mutations

multi_n = all_carr %>% group_by(id) %>% select(gene) %>% distinct() %>% summarise(n_gene = n())
multi_n = multi_n %>% filter(n_gene >= 2) %>% arrange(n_gene)
multi = left_join(multi_n, all_carr, by = "id")


write.table(multi, 
	"results/acmg_v3.1_clinvar_ni_carriers_multi_genes.txt", 
	sep = "\t", dec = ".", quote = FALSE, row.names = FALSE)


#### there are two BRCA2 variants (13:32340629:CTT:C) that are the same, but one was merged by chr_pos_allele and the other by aa_change. 
## the reason this happens is that there are actually 2 variants in the clinvar - one matches our one by chr,pos,allele and the other one where alleles don't match
## but it's the same amino acid change. therefore i would keep it as aa_change variant, while in fact it is the same as we already have. the pipeline should therefore 
## be improved by double-checking whether aa change variant already exists in chr_pos_allele list and if yes, it should be excluded from the aa_change list.

