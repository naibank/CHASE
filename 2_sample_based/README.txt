This is the pipeline for CHASE part 2

To run this, you need input from ../prerun_family_QCs/

1. run 1_run.R to generate job submission files
2. submit jobs using sh 1_submit_jobs.sh to identify all potential DelCH events
3. run 2_getSNVinCHEvent.R to identify SNVs make up DelCH events
4. run 3_analyzeCHEvent.R to perform statistical analysis

Description for additional files
1. 1_getCHEvent.R - an R script used by 1_run.R and 1_submit_jobs.sh to obtain all posible CH events
2. 2_submit_job.sh - a shell script to run 2_getSNVinCHEvent.R 
3. functions.R - R script contains required function for the analysis
4. run_template.sh - a shell script used as a template to generate multiple shell script files for batch job submission
