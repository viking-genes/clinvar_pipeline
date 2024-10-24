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

gene_files = list.files("results/variants/chr_pos_allele_merge/", pattern = "_clin_signif.txt")

# the file needs to be in the format: CHR POS GENE; no header, tab seperated
vars = lapply(gene_files, function(f) {
	gene = gsub("(.*)_ni_exomes.*", "\\1", f)
	d = fread(sprintf("results/variants/chr_pos_allele_merge/%s", f))
	d = d %>% select(Chromosome, Position) %>% mutate(gene = gene)
})
vars = do.call(rbind, vars)

write.table(vars, 
	"scripts/variants_to_extract.txt", 
	sep = "\t", dec = ".", quote = FALSE, row.names = FALSE, col.names = FALSE)
