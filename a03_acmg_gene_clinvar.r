#### merging exome variants with clinvar data and NI frequencies
rm(list = ls(all = TRUE))
library(data.table)
library(plyr)
library(dplyr)
library(tidyr)

args = commandArgs(TRUE)
gene = args[1]
run_path = args[2]

### getting new exome annotations (updated to dbsnp151)
exome_dbsnp151 = fread("data/exomes/EDINBURGH_Freeze_Two_remake.rep.dbsnp151.vcf.gz")
names(exome_dbsnp151)[1] = "CHROM"
# to save time extracting only the gene of interest
exome_dbsnp151 = exome_dbsnp151 %>% filter(grepl(gene, INFO))
# i just want to remap my variants to newer rsids, so keeping only necessary columns
exome_dbsnp151 = exome_dbsnp151 %>% select(CHROM, POS, ID, REF, ALT)

# multialleles are all in one row; need to split them into separate rows
exome_dbsnp151 = ddply(exome_dbsnp151, .(CHROM, POS, ID, REF), function(d) {
	data.frame(ALT = unlist(strsplit(d$ALT, split = ",")))
})
# also there are multiple rsids per row sometimes, splitting them up too
exome_dbsnp151 = ddply(exome_dbsnp151, .(CHROM, POS, ID, REF, ALT), function(d) {
	data.frame(rsid = unlist(strsplit(d$ID, split = ";")))
})
exome_dbsnp151 = exome_dbsnp151 %>% select(-ID) %>% rename(ID = rsid)

# set working directory 
setwd(run_path)

# load unfiltered data
d_exome = fread("data/regeneron_variants_gene_annotations_functional_consequence_new_parsing.txt")
d_exome = d_exome %>% filter(hgnc_name == gene)
# there can be multiple rsids per row here; 
# splitting them into individual rows because I don't necessarily know which of rsids will match with clinvar, so keeping them all
d_exome = ddply(d_exome, .(CHROM, POS, REF, ALT), function(d) {
	#writeLines(sprintf("%s", d$SNP))
	data.frame(
		SNP = d$SNP,
		hgnc_name = d$hgnc_name,
		ID = unlist(strsplit(d$ID, split = ";")),
		funct_conseq = d$funct_conseq, 
		isLof = d$isLof, 
		isAncestralAllele = d$isAncestralAllele, 
		deleteriousMissenseCount = d$deleteriousMissenseCount, 
		geneID = d$geneID, 
		fractionOfTranscriptsAffected = d$fractionOfTranscriptsAffected, 
		hgvs.c = d$hgvs.c, 
		hgvs.p = d$hgvs.p,
		INFO = d$INFO)
})

# merge with new dbsnp annotation
d_exome = left_join(d_exome, exome_dbsnp151, by = c("CHROM", "POS", "REF", "ALT"))

# again keep all possible rsids because I don't know which will fit with the clinvar
# first replace all "." with rsid if it exists in dbsnp151
d_exome$ID.x[which(d_exome$ID.x == ".")] = NA
d_exome$ID.y[which(d_exome$ID.y == ".")] = NA

d_exome = d_exome %>% unite("ID", ID.x,ID.y, na.rm = TRUE, remove = FALSE, sep = ";")

# replacing missing values back with .
d_exome$ID[which(d_exome$ID == "")] = "."


d_exome = ddply(d_exome, .(CHROM, POS, REF, ALT), function(d) {
	#writeLines(sprintf("%s", d$SNP))
	data.frame(
		SNP = d$SNP,
		hgnc_name = d$hgnc_name,
		ID = paste(unique(unlist(strsplit(d$ID, split = ";"))), collapse = ";"),
		funct_conseq = d$funct_conseq, 
		isLof = d$isLof, 
		isAncestralAllele = d$isAncestralAllele, 
		deleteriousMissenseCount = d$deleteriousMissenseCount, 
		geneID = d$geneID, 
		fractionOfTranscriptsAffected = d$fractionOfTranscriptsAffected, 
		hgvs.c = d$hgvs.c, 
		hgvs.p = d$hgvs.p,
		INFO = d$INFO)
})
d_exome = d_exome %>% distinct()

# keeping only coding variants
d_exome = d_exome %>% 
	filter(funct_conseq != "intronic" & funct_conseq != "downstream" & funct_conseq != "upstream" & funct_conseq != "3_prime_UTR" & funct_conseq != "5_prime_UTR")

# extract gene from clinvar table
d_clinvar = fread(sprintf("data/clinvar/per_gene_clinvar/%s_clinvar.txt", gene))

# turns out that for indels PositionVCF corresponds to the position in exome data. for a large majority of positions in the genome
# PositionVCF == Start. So I want to keep PositionVCF unless it's -1: then I want to keep start
d_clinvar = d_clinvar %>% mutate(newPos = PositionVCF)
d_clinvar$newPos[which(d_clinvar$PositionVCF == -1)] = d_clinvar$Start[which(d_clinvar$PositionVCF == -1)]
d_clinvar = d_clinvar %>% rename(Position = newPos)

# rename columns
d_exome <- d_exome %>%
	rename(Chromosome = CHROM) %>%
	rename(Position = POS) %>%
	rename(rsid = ID)

d_exome = d_exome %>% mutate(Chromosome = as.character(Chromosome))

d_clinvar <- d_clinvar %>%
	rename(rsid = `RS# (dbSNP)`) %>%
	rename(REF = ReferenceAlleleVCF) %>% 
	rename(ALT = AlternateAlleleVCF)
# filtering out build 38
d_clinvar = d_clinvar %>% filter(Assembly == "GRCh38")
# adding rs to rsid
d_clinvar = d_clinvar %>% mutate(rsid = paste0("rs", rsid))
d_clinvar = d_clinvar %>% mutate(Chromosome = as.character(Chromosome))

# merging by chromosome, position and alleles 
d_merge = left_join(d_exome, d_clinvar, by = c("Chromosome", "Position", "REF", "ALT"))

### clinical significance according to clinvar: https://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig/#:~:text=The%20primary%20goal%20of%20ClinVar,without%20interpreting%20the%20clinical%20significance
#Benign
#Likely benign
#Uncertain significance
#Likely pathogenic
#Pathogenic
#Likely pathogenic, low penetrance - As recommended by ClinGen for variants with decreased penetrance for Mendelian conditions.
#Pathogenic, low penetrance - As recommended by ClinGen for variants with decreased penetrance for Mendelian conditions.
#Uncertain risk allele - As recommended by ClinGen for variants with decreased penetrance for Mendelian conditions.
#Likely risk allele - As recommended by ClinGen for variants with decreased penetrance for Mendelian conditions.
#Established risk allele - As recommended by ClinGen for variants with decreased penetrance for Mendelian conditions.
#drug response - A general term for a variant that affects a drug response, not a disease.
#association - For variants identified in a GWAS study and further interpreted for their clinical significance.
#protective - For variants that decrease the risk of a disorder, including infections.
#Affects - For variants that cause a non-disease phenotype, such as lactose intolerance.
#conflicting data from submitters - Only for submissions from a consortium, where groups within the consortium have conflicting intepretations of a variant but provide a single submission to ClinVar.
#other - 	If ClinVar does not have the appropriate term for your submission, we ask that you submit "other" as clinical significance and contact us to discuss if there are other terms we should add.
#not provided


# filter merged data by creating new categories of clinical significance
d_merge = d_merge %>% mutate(clinical_significance_categories = "")

d_merge = d_merge %>% 
	mutate(clinical_significance_categories = ifelse(grepl("likely pathogenic", ClinicalSignificance, ignore.case = T), "likely pathogenic", clinical_significance_categories)) %>%
	mutate(clinical_significance_categories = ifelse(grepl("likely benign", ClinicalSignificance, ignore.case = T), "likely benign", clinical_significance_categories)) %>%
	mutate(clinical_significance_categories = ifelse(grepl("conflicting", ClinicalSignificance, ignore.case = T), "conflicting", clinical_significance_categories)) %>%
	mutate(clinical_significance_categories = ifelse(grepl("uncertain", ClinicalSignificance, ignore.case = T), "uncertain", clinical_significance_categories)) %>%
	mutate(clinical_significance_categories = ifelse(ClinicalSignificance == "Pathogenic", "pathogenic", clinical_significance_categories)) %>% 
	mutate(clinical_significance_categories = ifelse(ClinicalSignificance == "Benign", "benign", clinical_significance_categories)) %>% 
	mutate(clinical_significance_categories = ifelse(grepl("Benign,", ClinicalSignificance, ignore.case=T), "benign", clinical_significance_categories)) %>%
	mutate(clinical_significance_categories = ifelse(grepl("Pathogenic,", ClinicalSignificance, ignore.case=T), "pathogenic", clinical_significance_categories))

d_merge = d_merge %>%
	mutate(clinical_significance_categories = ifelse(clinical_significance_categories == "", "other", clinical_significance_categories))
	
#### create new categories of review status
### from: https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/
#four	practice guideline	practice guideline
#three	reviewed by expert panel	reviewed by expert panel
#two	criteria provided, multiple submitters, no conflicts	Two or more submitters with assertion criteria and evidence (or a public contact) provided the same interpretation.
#one	criteria provided, conflicting interpretations	Multiple submitters provided assertion criteria and evidence (or a public contact) but there are conflicting interpretations. The independent values are enumerated for clinical significance.
#one	criteria provided, single submitter	One submitter provided an interpretation with assertion criteria and evidence (or a public contact).
#none	no assertion for the individual variant	The allele was not interpreted directly in any submission; it was submitted to ClinVar only as a component of a haplotype or a genotype.
#none	no assertion criteria provided	The allele was included in a submission with an interpretation but without assertion criteria and evidence (or a public contact).
#none	no assertion provided	The allele was included in a submission that did not provide an interpretation.

d_merge = d_merge %>% mutate(star_status = "")

d_merge = d_merge %>% 
	mutate(star_status = ifelse(grepl("no assertion criteria provided", ReviewStatus, ignore.case = T), "0", star_status)) %>%
	mutate(star_status = ifelse(grepl("no assertion provided", ReviewStatus, ignore.case = T), "0", star_status)) %>%
	mutate(star_status = ifelse(grepl("no assertion for the individual variant", ReviewStatus, ignore.case = T), "0", star_status)) %>%
	mutate(star_status = ifelse(grepl("criteria provided, conflicting interpretations", ReviewStatus, ignore.case = T), "1", star_status)) %>%
	mutate(star_status = ifelse(grepl("criteria provided, single submitter", ReviewStatus, ignore.case = T), "1", star_status)) %>%
	mutate(star_status = ifelse(grepl("criteria provided, multiple submitters, no conflicts", ReviewStatus, ignore.case = T), "2", star_status)) %>%
	mutate(star_status = ifelse(grepl("reviewed by expert panel", ReviewStatus, ignore.case = T), "3", star_status)) %>%
	mutate(star_status = ifelse(grepl("practice guideline", ReviewStatus, ignore.case = T), "4", star_status)) 

# for those that have star_status == "" it's likely that they're not in Clinvar (not the exact chr:pos:ref:alt) - I checked manually for 
# ~10 missense variants

# fixing up rsids
d_merge = d_merge %>% mutate(rsid = gsub("\\.;", "", rsid.x)) 
d_merge$rsid[which(d_merge$rsid.x == "." & d_merge$rsid.y != "." & d_merge$rsid.y != "rs-1")] = d_merge$rsid.y[which(d_merge$rsid.x == "." & d_merge$rsid.y != "." & d_merge$rsid.y != "rs-1")]
d_merge = d_merge %>% select(-rsid.x, -rsid.y)

### saving the data: full version and only pathogenic/likely pathogenic with star >=1
d_merge_filter = d_merge %>% 
	filter(grepl("pathogenic", clinical_significance_categories)) %>%
	filter(star_status >= 1) %>%
	select(Chromosome, Position, REF, ALT, SNP, hgnc_name, rsid, funct_conseq, isLof, geneID, fractionOfTranscriptsAffected, hgvs.c, hgvs.p, 
		Type, Name, ClinicalSignificance, ClinSigSimple, PhenotypeList, ReviewStatus, NumberSubmitters, clinical_significance_categories, star_status)

write.table(d_merge_filter, sprintf("results/variants/chr_pos_allele_merge/%s_exomes_clin_signif.txt", gene), sep = "\t", dec = ".", quote = FALSE, row.names = FALSE)

# saving full data as a txt
d_merge = d_merge %>% select(-INFO)
write.table(d_merge, sprintf("results/variants/chr_pos_allele_merge/%s_exomes_clinvar_full.txt", gene), sep = "\t", dec = ".", quote = FALSE, row.names = FALSE) 

