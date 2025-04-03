This is the data preprocessing pipeline

#### perform sample QCs based on deletions detected by CNV pipeline (3IQR+Q3 | 3SD+mean), we didn't use deletions from SVs, as those are found to be mosaic and will be removed later
#### identify SVs/CNVs that are genomic disorders, impact ASD candidate genes, large >3MB and exclude those samples from the analysis
#### extract one proband and one unaffected sib per family
#### filter CNVs to CDS
#### filter SVs to those called by Manta and DELLY

To run this > submit 1_run_qc.sh

#### perform mosaic CNVs/SVs tagging
#### for each sample, find deletions and their SNVs within the deletion boundaries. Tag deletion with heterozygous SNVs as potential mosaic deletion and exclude them from the subsequent analysis

To run this > submit submit_jobs.sh

