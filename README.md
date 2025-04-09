# CHASE (Compound Heterozygosity in Autism Spectrum Etiology) 

This study has proposed three analytical procedures to investigate the role of deletion-SNV compound heterozygosity (DelCH) in the ASD etiology as following:

![image](https://github.com/user-attachments/assets/7c8f3cc1-9fa7-4188-a618-37f1c7a6ddb9)

## Strategy 1: Deletion matched
1) Burden analysis of SNVs in probands and deletion-transmitting parents considering only inherited deletions (via conditional logistic regression stratified by deletion) - /1_inherited_del/
## Strategy 2 All deletions and all SNVs
3) Burden analysis of SNVs in probands and both parents considering all deletions (via conditional logistic regression stratified by family) - /2_sample_based/
## Strategy 3: TDT
5) Transmission Disequilibrium Test (TDT) of SNVs in deletion-non-transmitting parents where SNVs would make up DelCH events (via Fisher's exact test) - /3_TDT_analysis/

Please see README file in each subdirectory for more information about each analysis.

## Control Quality Steps
Before running any of the analysis, preprocessing steps for sample QCs and variants QCs are required, which can be done using scripts located in /prerun_family_QC/ folder.
