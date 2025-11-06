# CHASE (Compound Heterozygosity in Autism Spectrum Etiology) 

This study has proposed three analytical procedures to investigate the role of deletion-SNV compound heterozygosity (DelCH) in the ASD etiology as following:

<img width="1479" height="1446" alt="CHASEmainschemeLastVersion pptxcut" src="https://github.com/user-attachments/assets/6fa8ed63-486d-42fc-b16d-740ec667cb8a" />

## Strategy 1: Deletion matched
Burden analysis of SNVs in probands and deletion-transmitting parents considering only inherited deletions (via conditional logistic regression stratified by deletion) - /1_inherited_del/
## Strategy 2 All deletions and all SNVs
Burden analysis of SNVs in probands and both parents considering all deletions (via conditional logistic regression stratified by family) - /2_sample_based/
## Strategy 3: TDT
Transmission Disequilibrium Test (TDT) of SNVs in deletion-non-transmitting parents where SNVs would make up DelCH events (via Fisher's exact test) - /3_TDT_analysis/

**Please see README file in each subdirectory for more information about each analysis.**

## Control Quality Steps
Before running any of the analysis, preprocessing steps for sample QCs and variants QCs are required, which can be done using scripts located in /prerun_family_QC/ folder.

## Miscellaneous files
/recessive_genes contains list of recessive genes obtained from GenomicsEngland gene panel for NDD genes, and our in-house manual curation as described in the paper.

/gene_info contains gnomAD gene contraint information and gene definition from RefSeq release 200 (GRCh38). 
