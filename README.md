Pipeline for extracting carriers from the sequencing data for variants of interest.

If using these scripts or parts of these scripts please cite XXX

Many thanks to Ashwini Shanmugam for helping to put this documentation in place.

# Pipeline description

This pipeline was originally built to extract carriers of ACMG actionable variant carriers from the exome sequencing data, but can be used to extract carriers of any genes. 

The pipeline primarily uses bash scripting to extract gene-based information from large WES/WGS and ClinVar files and to run various R scripts. R scripts are used to merge the files together and ensure that only variants of interest are extracted. Briefly, the pipeline takes as an input a list of genes, one at a time, extracts their data from the ClinVar and whole exome/genome sequence annotation files. The two are then merged and filtered to find which pathogenic/likely pathogenic variants are present in the WES/WGS data. The IDs of carriers of these variants are then extracted from the WES/WGS data. As the last step of the pipeline, carriers of these variants are then extracted from the whole exome/genome sequencing files. 
The pipeline is designed to take one gene from the gene_list.txt file at a time, extract all the info on that gene from the ClinVar and whole exome/genome annotation file, followed by mergeing of the gene-based files and extracting overlapping variants that satisfy certain criteria. As the last step of the pipeline, carriers of these variants are then extracted from the whole exome/genome sequencing files.

The pipeline was built to run one gene at a time, enabling parallelisation of the process. On our high-performance computing clusters (HPC) this is achieved by running array of jobs at once, where each job index (numbers from 1-number of genes) is used to extract corresponding line from a file containing a list of genes of interest. This is reflected in all bash (.sh) scripts in this pipeline, where SGE_TASK_ID is used to define the number of the line that needs to be extracted from the file listing all the genes. Different HPC set-ups might require different handling of the parallelisation procedure.

Scripts a01-a07 are needed for the purpose of running the pipeline. Scripts a08 and a09 were used to summarise and further annotate the data for the paper.


# Brief description of scripts: 

### a01_extract_gene_variants_clinvar.sh
  > Splits the ClinVar file into per-gene files for easier handling.
### a02_run_acmg_gene_clinvar.sh
  > Runs a03_acmg_gene_clinvar.r script for each gene of interest.
### a02_run_acmg_gene_clinvar_merge_by_aa_change.sh
  > Runs a03_acmg_gene_clinvar_merge_by_aa_change.r script for each gene of interest.
### a03_acmg_gene_clinvar.r
  > Creates a list of variants from the CLINVAR data which are in genes from a given gene list (eg. ACMG).
  > 
  > Extracts all variants in a given gene from the exome/genome variant annotation file, merges them using chromosome, position and alleles with the ClinVar annotations and filters based on pathogenicity/clinical significance.
  >
  > Additionally, filters for variants that are clinically relevant:
  > - Annotated as “pathogenic” or “likely pathogenic”
  > - Variant submitted by multiple submitters or reviewed by an expert panel or has practical guidelines according to ClinVar
  > - Assigns stars based on review status in ClinVar
  >   - four stars: practice guideline practice guideline
  >   - three stars: reviewed by expert panel reviewed by expert panel
  >   - two stars: criteria provided, multiple submitters, no conflicts Two or more submitters with assertion criteria and evidence (or a public contact) provided the same interpretation
  >   - one star: criteria provided, conflicting interpretations Multiple submitters provided assertion criteria and evidence (or a public contact) but there are conflicting interpretations. The independent values are enumerated for clinical significance. Alternatively, criteria provided, single submitter One submitter provided an interpretation with assertion criteria and evidence (or a public contact)
  > - Creates a complete list of variants present in Clinvar and the dataset. Creates a subset of variants within the gene which are clinically significant; i.e. 2 star review status and above, and are pathogenic or likely pathogenic
     
### a03_acmg_gene_clinvar_merge_by_aa_change.r
  > The same as above a03 script, but instead of mergeing variants based on chromosome, position and alleles, it merges variants based on the amino acid changed caused by the variant

### a04_run_create_list_of_variants.sh
  > Runs a05_create_list_of_variants.r

### a04_run_create_list_of_variants_aa_change.sh
  > Runs a05_create_list_of_variants_aa_change.r

### a05_create_list_of_variants.r
  > Extracts the list of all pathogenic and clinically significant variants in your data and prepares a list in the format chromosome; position; HGNC gene name which is needed as input for extracting carriers in the data.
  >
  > Basically, it concatenates a list of variants that are in the clinically significant variant files.
  >
  > These variants will be extracted from the cohort VCF files based on chromosome coordinates.

### a05_create_list_of_variants_aa_change.r
  > The same as the above a05 script, but with variants merged based on the amino acid change.

### a06_extract_variant_genotypes.sh
  > Creates directories in results directory for each variant in the shortlisted variant list from the previous step.
  >
  > Extracts variant from cohort vcf files.
  >
  > Identifies carriers of variants in cohort in homozygous and heterozygous states (output as two separate files).

### a06_extract_variant_genotypes_aa_change_merge.sh
  > The same as the above a06 script, but with variants merged based on the amino acid change.

### a07_extract_carriers.r
  > Reads the files that contain list of all carriers for a given variant (output of a06 scripts) and reformats it from wide into a long format, with each row being one ID.
  
### a07_extract_carriers_aa_change_merge.r
  > The same as the above a07 script, but with variants merged based on the amino acid change.

### a08_summarise.r
  > Summarises the number of carriers for each variant and checks if all variants were extracted; it also annotates all the variants with frequencies.

### a08_summarise_aa_change_merge.r 
  > The same as the above a08 script, but with variants merged based on the amino acid change.
  >
  > Compiles data from the heterozygous and homozygous carriers lists output.
  >
  > Creates a summary of the number of het and hom carriers overall and within each dataset.
  >
  > Annotated with gnomAD frequencies.

### a09_compound_hets.r
  > Loads the carrier data saved in a08 and checks if any of the carriers carry two different variants for the same gene (compound hets).

# Requirements 

## Data

> The following data are needed to run this pipeline:
> 1. genetic data - whole exome or genome sequencing data.
> 2. ClinVar data - a file containing annotations of pathogenicity of variants. ClinVar database is regularly updated, so it’s important to keep a record of the version of the data and the date of accession, since variant’s pathogenicity and star status can be significantly changed between two data releases. Location of the file: data/clinvar/variants_summary.txt.gz The file can be downloaded from: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/
> 3. list of all variants in the genetic data annotated with the following info: chromosome, position, reference allele, alternate allele, functional consequence
> 4. List of genes of interest - This list needs to contain the HUGO name of the gene, each gene in one row and should not contain the header. 







