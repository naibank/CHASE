# This is the pipeline for CHASE strategy 1- Deletion matched:
Burden Analysis of DelCH events in inherited deletions: A group-level comparison between probands and their deletion-transmitting parents. 

To run this, you need input from ../prerun_family_QCs/

1. run 1_run.R to generate job submission files
2. submit jobs using sh 1_submit_jobs.sh
3. run 2_getSNVinCHEvent.R to identify SNVs make up DelCH events
4. run 3_analyzeCHEvent.R to perform statistical analysis
