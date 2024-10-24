### counting the number of variant carriers
### script started 18/11/2019
### script edited 27032023 to deal with multiallelic variants
### script written lklaric
library(data.table)
library(plyr)
library(dplyr)
library(tidyr)

				
args = commandArgs(TRUE)
folder = args[1]
file_name = args[2]				
run_path = args[3]

				
d = fread(paste0(run_path, "/results/genotypes/chr_pos_allele_merge/", folder, "/", file_name, ".txt"))
d_t = data.frame(t(d))
ref_a = as.character(d_t["REF",])
alt_a = as.character(d_t["ALT",])

d_t$id = toupper(row.names(d_t))

d_t = d_t %>% mutate(id = gsub(".*(ORCA\\d+)", "\\1", id))
d_t = d_t %>% mutate(id = gsub(".*(VIKI\\d+)", "\\1", id))


######## new code, designed to deal with multiallelic variants: 
### from vcf descriptors
# AF="Allele Frequency estimate for each alternate allele"
# AQ="Allele Quality score reflecting evidence for each alternate allele (Phred scale)"
# GT="Genotype"
# RNC="Reason for No Call in GT: . = n/a, M = Missing data, P = Partial data, D = insufficient Depth of coverage, - = unrepresentable overlapping deletion, L = Lost/unrepresentable allele (other than deletion), U = multiple Unphased variants present, O = multiple Overlapping variants present, 1 = site is Monoallelic, no assertion about presence of REF or ALT allele">
# DP="Read Depth"
# AD="Allelic depths for the ref and alt alleles in the order listed"
# SBPV="Strand bias P-value; probability that the fraction of forward reads (VCF) amongst reads supporting alt allele (VC) is more extreme than expected assuming a beta-binomial distribution.">
# GQ="Genotype Quality">
# PL="Phred-scaled genotype Likelihoods"
# FT="FILTER field from sample gVCF">
# FORMAT=GT:DP:AD:SBPV:GQ:PL:FT:RNC

d_t$gt = sapply(strsplit(d_t$t.d., split = ":"), function(k) k[1])
d_t$dp = sapply(strsplit(d_t$t.d., split = ":"), function(k) k[2])
d_t$ad = sapply(strsplit(d_t$t.d., split = ":"), function(k) k[3])
d_t$sbpv = sapply(strsplit(d_t$t.d., split = ":"), function(k) k[4])
d_t$gq = sapply(strsplit(d_t$t.d., split = ":"), function(k) k[5])
d_t$pl = sapply(strsplit(d_t$t.d., split = ":"), function(k) k[6])
d_t$ft = sapply(strsplit(d_t$t.d., split = ":"), function(k) k[7])
d_t$rnc = sapply(strsplit(d_t$t.d., split = ":"), function(k) k[8])

d_out = ddply(d_t, .(id), function(k) {
	data.frame(alleles = c(ref_a, unlist(strsplit(alt_a, split = ","))), 
		allele_type = c("REF", rep("ALT", length(unlist(strsplit(alt_a, split = ","))))),
		gt = k$gt, 
		dp = k$dp, 
		ad = unlist(strsplit(k$ad, split = ",")),
		sbpv = k$sbpv,
		gq = k$gq, 
		pl = k$pl, 
		ft = k$ft, 
		rnc = k$rnc, 
		full_info = k$t.d.)
})

# figuring out which alternate allele is present in the case of multiallelic variants - its allelic depth should be different than 0
hets = d_out %>% filter(gt == "0/1")
hets = hets %>% filter(ad != 0)
hets = hets %>% spread(allele_type, alleles)
# making sure I have unified columns with REF and ALT
hets = ddply(hets, .(id), function(k) {
	k$new_alt = unique(k$ALT[!is.na(k$ALT)])
	k$new_ref = unique(k$REF[!is.na(k$REF)])
	k
})
# clearing up old columns; i can safely remove ad because I still keep this info in the full_info column
hets = hets %>% select(-ALT, -REF, -ad) %>% distinct() 
hets = hets %>% rename(REF = new_ref, ALT = new_alt)
# keeping info on all alt alleles
hets = hets %>% mutate(all_ALT = alt_a)

write.table(hets, 
	paste0(run_path, "/results/genotypes/chr_pos_allele_merge/", file_name, sprintf("_%s_%s_hets.txt", ref_a, alt_a)), 
	sep = "\t", dec = ".", quote = FALSE, row.names = FALSE)


### homozygotes. in theory I should just do this thing on the d_out and then split by genotype but that requires more checking
### and i don't currently have time for that
homo = d_out %>% filter(gt == "1/1")
if (nrow(homo) >=1) {
	homo = homo %>% filter(ad != 0)
	homo = homo %>% spread(allele_type, alleles)
	# making sure I have unified columns with REF and ALT
	homo = ddply(homo, .(id), function(k) {
		if (length(unique(k$ALT[!is.na(k$ALT)]) >= 1)) {
			k$new_alt = unique(k$ALT[!is.na(k$ALT)])	
		} else {
			k$new_alt = NA
		}

		if (length(unique(k$REF[!is.na(k$REF)]) >= 1)) {
			k$new_ref = unique(k$REF[!is.na(k$REF)])	
		} else {
			k$new_ref = NA
		}
		k
	})
	# clearing up old columns; i can safely remove ad because I still keep this info in the full_info column
	homo = homo %>% select(-ALT, -REF, -ad) %>% distinct() 
	homo = homo %>% rename(REF = new_ref, ALT = new_alt)
	# keeping info on all alt alleles
	homo = homo %>% mutate(all_ALT = alt_a)

	write.table(homo, 
		paste0(run_path, "/results/genotypes/chr_pos_allele_merge/", file_name, sprintf("_%s_%s_homo.txt", ref_a, alt_a)), 
		sep = "\t", dec = ".", quote = FALSE, row.names = FALSE)
}

######### \new code



