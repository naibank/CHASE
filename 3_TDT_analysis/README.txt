This is the pipeline for CHASE part 3

To run this, you need input from ../prerun_family_QCs/

1. submit jobs using 1_get_all_SNVs.sh to obtain SNVs from your input dataset(s) using 1_get_SNVs.R
2. submit jobs using 2_run_get_CNV_SNV_table.sh to identify all DelCH events using 2_get_CNV_SNV_table.R
3. run 3_process_CNV_SNV_table.R to add event frequency to the DelCH events
4. run 4_get_FisherTest_eventfreq_combine_changeFilter.R to perform statistical analysis
