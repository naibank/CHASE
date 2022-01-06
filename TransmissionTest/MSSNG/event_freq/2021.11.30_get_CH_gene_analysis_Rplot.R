library(data.table)

for (group in c("MSSNG_CH_all", "MSSNG_CH_MPX", "MSSNG_CH_SPX", "SSC_CH_proband", "SSC_CH_unaffected_sibling")){
  ## import gene-set enrichment analysis results
  res <- fread(sprintf("/Users/shaniawu/SickKids CHASE/CHASE Meeting/2021.12.01 data/%s_gene_enrichment_analysis_gsMain.tsv", group), 
                   data.table = F)
  
  ## remove insignificant results (p >= 0.05)
  sig_res <- res[which(res$P < 0.05),]
  
  ## filter by pRec and event_freq
  for(pRecx in c(0, 0.9)){ # pRec cut-offs
    for(freq in c(0.01, 0.005, 0.0001, 0.00001)){ # event freq cut-offs; no cut-off = 0.005 (MSSNG), 0.01 (SSC)
      if ((grepl("MSSNG", group) & freq == 0.01)|(grepl("SSC", group) & freq == 0.0005)){
        next
      }
      
      sig_res_tmp <- sig_res[which(sig_res$event_freq == freq &
                                 sig_res$pRec == pRecx),]
      
      if (nrow(sig_res_tmp) == 0){ # skip empty results
        next
      }
      
      ## plot OR vs. gene-set for pRec and event_freq
      plot_name =  sprintf("%s_eventfreq%s_pRec%s", group, signif(freq, digits = 3), signif(pRecx, digits = 3))
      
      plot <- ggplot(data = sig_res_tmp, aes(x = OR, y = reorder(gene_set, -OR), fill = variant_type))
      plot + geom_bar(stat = "identity", position = position_dodge(width = 1)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x = "odds ratio (OR)", y = "gene_set") +
        ggtitle(plot_name)
      
      ggsave(sprintf("%s.png", plot_name), width = 10)
    }
  }
}