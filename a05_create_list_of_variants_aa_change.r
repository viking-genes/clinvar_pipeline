### extracting the list of 2* pathogenic variants to extract
### script started lklaric
### script started 16032023
rm(list = ls(all = TRUE))
library(data.table)
library(plyr)
library(dplyr)
library(tidyr)

args = commandArgs(TRUE)
run_path = args[1]


setwd(run_path)

gene_files = list.files("results/variants/aa_change_merge/", pattern = "_clin_signif_aa_change.txt")

# the file needs to be in the format: CHR POS GENE; no header, tab seperated
vars = lapply(gene_files, function(f) {
	gene = gsub("(.*)_ni_exomes.*", "\\1", f)
	d = fread(sprintf("results/variants/aa_change_merge/%s", f))
})
vars = do.call(rbind, vars)

# the file needs to be in the format: CHR POS GENE; no header, tab seperated
vars_list = vars %>% select(Chromosome, Position.exome, hgnc_name) %>% rename(gene = hgnc_name) %>% distinct()

write.table(vars_list, 
	"scripts/variants_to_extract_aa_change_merge.txt", 
	sep = "\t", dec = ".", quote = FALSE, row.names = FALSE, col.names = FALSE)
