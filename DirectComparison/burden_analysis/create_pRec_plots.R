library(yaml)
library(survival)
library(data.table)
library(dplyr)
library(ggplot2)

setwd("/Users/faraz/Documents/Work/Workterm 5 - SickKids/CHASE Project/Faraz/Data")

#Load("../Unbiased/unbiased_synchevent_covariate.RData")

Generate_pRec_Distribution_Plot <- function(meta,
                                            SNV_data,
                                            name,
                                            title,
                                            proband_IDs,
                                            father_IDs,
                                            mother_IDs,
                                            save_path='../DT/pRec Plots/',
                                            plot_type='Box') {
  sj_df <- SNV_data[SNV_data$X.Sample %in% meta$Sample.ID & !is.na(SNV_data$gnomAD_pRec)
                    , c('X.Sample', 'gnomAD_pRec')]
  sj_df <- sj_df %>% group_by(X.Sample) %>% summarize(X.Sample=X.Sample, gnomAD_pRec=max(gnomAD_pRec)) %>% distinct(X.Sample, .keep_all=T)
  sj_df <- data.frame(sj_df)
  sj_df$Group <- 'NA'
  sj_df$Group[sj_df$X.Sample %in% proband_IDs] <- 'Case'
  sj_df$Group[sj_df$X.Sample %in% c(father_IDs, mother_IDs)] <- 'Control'
  sj_df <- sj_df[sj_df$Group != 'NA',]
  
  sj_df$Relation <- 'NA'
  sj_df$Relation[sj_df$X.Sample %in% proband_IDs] <- 'Proband'
  sj_df$Relation[sj_df$X.Sample %in% father_IDs] <- 'Father'
  sj_df$Relation[sj_df$X.Sample %in% mother_IDs] <- 'Mother'
  sj_df <- sj_df[sj_df$Relation != 'NA',]
  
  summary_df <- sj_df %>% group_by(Group) %>% summarize(m=mean(gnomAD_pRec))
  
  df = count(sj_df, Group)
  
  sj_df <- sj_df[sj_df$gnomAD_pRec >= 0.626,] ## Cut-off
  
  if (plot_type == 'Box') {
    ggplot(sj_df, aes(x=Group, y=gnomAD_pRec, color=Group)) +
      geom_boxplot(notch = T, outlier.shape = NA) +
      stat_summary(fun=mean, color='black') +
      scale_x_discrete(limits=c('Case','Control')) +
      coord_cartesian(ylim=quantile(sj_df$gnomAD_pRec, c(0.1, 0.9))) +
      coord_cartesian(ylim=c(0.6,1.05)) +
      geom_text(data=df,aes(y=0.6, label=paste('n =',n)),color='black') +
      geom_text(data=summary_df,aes(y=1.05, label=paste('Mean =',round(m,2))), color='black') +
      geom_jitter() +
      labs(y='pRec', x='Group', title=title) +
      theme(plot.title=element_text(hjust=0.5))
    
    ggsave(paste(save_path, paste(name,'.png')))
  }
  else if (plot_type == 'Scatter') {
    sj_df_scatter <- data.frame()
    sj_df_scatter$proband_pRec <- sj_df$gnomAD_pRec[sj_df$Relation=='Proband']
    sj_df_scatter$father_pRec <- sj_df$gnomAD_pRec[sj_df$Relation=='Father']
    sj_df_scatter$mother_pRec <- sj_df$gnomAD_pRec[sj_df$Relation=='Mother']
    
    ggplot(sj_df_scatter, aes(x=father_pRec, y=proband_pRec)) +
      geom_point() +
      labs(y='Proband pRec', x='Father pRec', title=title) +
      theme(plot.title=element_text(hjust=0.5))
    
    ggsave(paste(save_path, paste(name,'_father.png')))
    
    ggplot(sj_df, aes(x=mother_pRec, y=proband_pRec)) +
      geom_point() +
      labs(y='Proband pRec', x='Mother pRec', title=title) +
      theme(plot.title=element_text(hjust=0.5))
    
    ggsave(paste(save_path, paste(name,'_mother.png')))
  }
  else if (plot_type == 'Density') {
    ggplot(sj_df, aes(x=gnomAD_pRec, color=Group)) +
      geom_density() +
      geom_vline(data=summary_df, aes(xintercept=m, color=Group), linetype='dashed') +
      labs(y='Density', x='pRec', title=title) +
      theme(plot.title=element_text(hjust=0.5))
    
    ggsave(paste(save_path, paste(name,'_density.png')))
  }
}

Generate_pRec_Distribution_Plot(MSSNG_meta,
                                MSSNG_parent_proband_proc_SNVs,
                                'MSSNG_parent_proband_proc_SNVs',
                                "MSSNG parent-proband",
                                MSSNG_proband_IDs,
                                MSSNG_father_IDs,
                                MSSNG_mother_IDs)
Generate_pRec_Distribution_Plot(SSC_meta,
                                SSC_parent_proband_SNVs,
                                'SSC_parent_proband_proc_SNVs',
                                "SSC parent-proband",
                                SSC_proband_IDs,
                                SSC_father_IDs,
                                SSC_mother_IDs)

Generate_pRec_Distribution_Plot(MSSNG_meta,
                                MSSNG_parent_proband_proc_SNVs,
                                'MSSNG_parent_proband_proc_SNVs_0.9',
                                "MSSNG parent-proband",
                                MSSNG_proband_IDs,
                                MSSNG_father_IDs,
                                MSSNG_mother_IDs)
Generate_pRec_Distribution_Plot(SSC_meta,
                                SSC_parent_proband_SNVs,
                                'SSC_parent_proband_proc_SNVs_0.9',
                                "SSC parent-proband",
                                SSC_proband_IDs,
                                SSC_father_IDs,
                                SSC_mother_IDs)

Generate_pRec_Distribution_Plot(MSSNG_meta,
                                MSSNG_parent_proband_proc_SNVs,
                                'MSSNG_parent_proband_proc_SNVs_0.626',
                                "MSSNG parent-proband",
                                MSSNG_proband_IDs,
                                MSSNG_father_IDs,
                                MSSNG_mother_IDs)
Generate_pRec_Distribution_Plot(SSC_meta,
                                SSC_parent_proband_SNVs,
                                'SSC_parent_proband_proc_SNVs_0.626',
                                "SSC parent-proband",
                                SSC_proband_IDs,
                                SSC_father_IDs,
                                SSC_mother_IDs)
#### Scatter plots by parent ####
Generate_pRec_Distribution_Plot(MSSNG_meta,
                                MSSNG_parent_proband_proc_SNVs,
                                'MSSNG_parent_proband_proc_SNVs_0.626',
                                "MSSNG parent-proband",
                                MSSNG_proband_IDs,
                                MSSNG_father_IDs,
                                MSSNG_mother_IDs,
                                plot_type='Scatter')
Generate_pRec_Distribution_Plot(SSC_meta,
                                SSC_parent_proband_SNVs,
                                'SSC_parent_proband_proc_SNVs_0.626',
                                "SSC parent-proband",
                                SSC_proband_IDs,
                                SSC_father_IDs,
                                SSC_mother_IDs,
                                plot_type='Scatter')

#### Density plots ####
Generate_pRec_Distribution_Plot(MSSNG_meta,
                                MSSNG_parent_proband_proc_SNVs,
                                'MSSNG_parent_proband_proc_SNVs_0.626',
                                "MSSNG parent-proband",
                                MSSNG_proband_IDs,
                                MSSNG_father_IDs,
                                MSSNG_mother_IDs,
                                plot_type='Density')
Generate_pRec_Distribution_Plot(SSC_meta,
                                SSC_parent_proband_SNVs,
                                'SSC_parent_proband_proc_SNVs_0.626',
                                "SSC parent-proband",
                                SSC_proband_IDs,
                                SSC_father_IDs,
                                SSC_mother_IDs,
                                plot_type='Density')
