# CHASE (Compound Heterozygosity in Autism Spectrum Etiology) 

This study has proposed three analytical procedures to investigate the role of deletion-SNV compound heterozygosity (DelCH) in the ASD etiology as following:
1) Burden analysis of SNVs in probands and deletion-transmitting parents considering only inherited deletions (via conditional logistic regression stratified by deletion) - /1_inherited_del/
2) Burden analysis of SNVs in probands and both parents considering all deletions (via conditional logistic regression stratified by family) - /2_sample_based/
3) Transmission Disequilibrium Test (TDT) of SNVs in deletion-non-transmitting parents where SNVs would make up DelCH events (via Fisher's exact test) - /3_TDT_analysis/

Please see README file in each subdirectory for more information about each analysis.

Before running any of the analysis, preprocessing steps for sample QCs and variants QCs are required, which can be done using scripts located in /prerun_family_QC/ folder.
