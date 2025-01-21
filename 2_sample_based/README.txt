This is the pipeline for CHASE part 2

To run this, you need input from ../prerun_family_QCs/

1. run 1_run.R to generate job submission files
2. submit jobs using sh 1_submit_jobs.sh to identify all potential DelCH events
3. run 2_getSNVinCHEvent.R to identify SNVs make up DelCH events
4. run 3_analyzeCHEvent.R to perform statistical analysis
