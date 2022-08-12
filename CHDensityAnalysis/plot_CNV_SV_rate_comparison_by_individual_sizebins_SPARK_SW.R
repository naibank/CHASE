################################################################################################################################################################################
# 
# plot_CNV_SV_rate_comparison_by_individual_sizebins_SPARK_SW.R
# purpose: - outputs SPARK sliding window cut-offs (median, 80, 90 perc.) spline tables
#          - plots SPARK sliding window cut-offs (median, 80, 90 perc.) spline
# input: SPARK CH_hits tables
#        /hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/CH_count_del_size/data/CH_hits/
# output: SPARK.2.5kb.n7.pRec.window.cutoffs.df, SPARK.all.cutoff_2.5kb window.n=7.pdf
#         /hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/CH_count_del_size/figures/ 
# 
# notes: incomplete SPARK deletion data (only CNVs, no SVs)
#
##############################################################################################################################################################################


setwd("/Users/shaniawu/hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/CH_count_del_size")

default <- c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF", "#0072B2", "#009E73", "#E69F00", "#CC79A7",  "#999999")

plot(NULL, xlim=c(0,length(default)), ylim=c(0,1), 
     xlab="", ylab="", xaxt="n", yaxt="n")
rect(0:(length(default)-1), 0, 1:length(default), 1, col=default)


## Read CH hits 
SPARK.CH.hits <- fread("./data/CH_hits/SPARK_CH_hits_nofilt_sex.tsv", data.table=F)
SPARK.CH.hits.pRec <- fread("./data/CH_hits/SPARK_CH_hits_pRec_sex.tsv", data.table=F)


########################################################################################################################################################################################

#### Get CUT-OFFs for each size bin function ####

Get_Sliding_CutOffs_Table <- function(CH_hits, bin.end=F, bin.size, window.bin.n){
  ## remove (0,0) individuals from CH_hits (outliers)
  CH_hits <- CH_hits[which(!(CH_hits$exSize == 0 & CH_hits$CH_hit == 0)),]
  
  ## make size.bins list
  size.bins <- list()
  bin.start <- 0
  if (bin.end == F){
    bin.end <- max(CH_hits$exSize)
  }
  while (bin.start < bin.end){
    bin <- list(bin.start, bin.start + bin.size) # list(start, end)
    size.bins <- c(size.bins, list(bin))
    bin.start <- bin.start + bin.size
  }
  
  ## find cut-offs for each SLIDING size bin
  cutoffs.out <- data.frame()
  
  for (i in 1:length(size.bins)){ # find start/end size bin cut-offs
    if (bin.size == 2500){
      alt.start.i <- window.bin.n/2 + 0.5
      alt.end.i <- length(size.bins) - alt.start.i + 1
      if (i <= alt.start.i){
        window.start.i <- 1
        window.end.i <- window.bin.n
      }
      if (i >= alt.end.i){
        window.start.i <- length(size.bins) - window.bin.n + 1 ## does not work across all n !
        window.end.i <- length(size.bins)
      }
      if ((i > alt.start.i)&(i < alt.end.i)){
        window.start.i <- i - (window.bin.n - 1)/2
        window.end.i <- i + (window.bin.n - 1)/2
      }
    }
    # if (bin.size == 500){
    #   if (i <= beg.i){ # For the first bin, quantile(c((n,n+1,n+2,n+3,n+4,n+5,n+6),c(.5,.8,.9))
    #     # and increasing until you reach bin until you reach the main formula.
    #     window.start.i <- 1
    #     window.end.i <- window.bin.n
    #   }
    #   if (i >= (length(size.bins)-(beg.i-1))){ # For the last bin, quantile(c((n-6,n-5,n-4,n-3,n-2,n-1,n),c(.5,.8,.9))
    #     window.start.i <- length(size.bins) - window.bin.n
    #     window.end.i <- length(size.bins)
    #   }
    #   if ((i > beg.i) & (i < (length(size.bins)-(beg.i-1)))){ # quantile(c((n-3,n-2,n-1,n,n+1,n+2,n+3),c(.5,.8,.9))
    #     window.start.i <- i-beg.i
    #     window.end.i <- i+beg.i+1
    #   }
    # }
    start.exsize <- size.bins[[i]][[1]] # non-sliding start exsize of size bin 
    end.exsize <- size.bins[[i]][[2]] # non-sliding end exsize of size bin 
    sliding.window.start.exsize <- size.bins[[window.start.i]][[1]]
    sliding.window.end.exsize <- size.bins[[window.end.i]][[2]]
    
    ## filter CH hits for only those within sliding window exsize
    CH_hits_in_sliding_bin <- CH_hits[which(CH_hits$exSize >= sliding.window.start.exsize &
                                      CH_hits$exSize < sliding.window.end.exsize),] 
    
    ## add size bin start & end columns
    cutoffs.df <- data.frame(size.bin.start=start.exsize)
    cutoffs.df$size.bin.end <- end.exsize
    
    ## add cutoff columns
    cutoffs.df$median <- median(CH_hits_in_sliding_bin$CH_hit)
    cutoffs.df$mean <- mean(CH_hits_in_sliding_bin$CH_hit)
    cutoffs.df$percentile_0.8 <- quantile(CH_hits_in_sliding_bin$CH_hit, probs = 0.8)
    cutoffs.df$percentile_0.9 <- quantile(CH_hits_in_sliding_bin$CH_hit, probs = 0.9)
    cutoffs.df$bin.count <- nrow(CH_hits_in_sliding_bin)
    
    cutoffs.out <- rbind(cutoffs.out, cutoffs.df)
  }
  return(cutoffs.out)
}

## SPARK bin size 2.5 kb cutoffs table
SPARK.2.5kb.n7.window.cutoffs.df <- Get_Sliding_CutOffs_Table(SPARK.CH.hits, 
                                                               bin.end=30000, bin.size = 2500,
                                                               window.bin.n = 7)
SPARK.2.5kb.n5.window.cutoffs.df <- Get_Sliding_CutOffs_Table(SPARK.CH.hits, 
                                                                  bin.end=30000, bin.size = 2500,
                                                                  window.bin.n = 5)

SPARK.2.5kb.n7.pRec.window.cutoffs.df <- Get_Sliding_CutOffs_Table(SPARK.CH.hits.pRec, 
                                                                  bin.end=30000, bin.size = 2500,
                                                                  window.bin.n = 7)

# write.table(SPARK.2.5kb.n7.window.cutoffs.df, "./data/SPARK.2.5kb.n7.window.cutoffs.df.tsv",
#             sep="\t", row.names=F, quote=F, col.names=T)
# write.table(SPARK.2.5kb.n7.pRec.window.cutoffs.df, "./SPARK.2.5kb.n7.pRec.window.cutoffs.df.tsv",
#             sep="\t", row.names=F, quote=F, col.names=T)


########################################################################################################################################################################################
#### Plot CUT-OFFs for each size bin - sliding window approach ####

Get_Sliding_Window_Cutoff_Plots <- function(window.cutoffs.df, CH.hits.df, file.name){
  ### Plots median, mean, 80th, and 90th percentiles for each 500 bp size bin 
  ###   using sliding window approach (takes average of 3 size bins above and below)

  ## Take average of 20-30kb size bins as single bin
  window.cutoffs.df.20kb.less <- window.cutoffs.df[which(window.cutoffs.df$size.bin.start <= 20000),]
  window.cutoffs.df.20kb.more <- window.cutoffs.df[which(window.cutoffs.df$size.bin.start >= 20000),]
  
  window.cutoffs.df.20kb.less[41,]$size.bin.end <- 30000
  window.cutoffs.df.20kb.less[41,]$median <- mean(window.cutoffs.df.20kb.more$median)
  window.cutoffs.df.20kb.less[41,]$mean <- mean(window.cutoffs.df.20kb.more$mean)
  window.cutoffs.df.20kb.less[41,]$percentile_0.8 <- mean(window.cutoffs.df.20kb.more$percentile_0.8)
  window.cutoffs.df.20kb.less[41,]$percentile_0.9 <- mean(window.cutoffs.df.20kb.more$percentile_0.9)
  window.cutoffs.df.20kb.less[41,]$bin.count <- sum(window.cutoffs.df.20kb.more$bin.count)
                                         
  sliding.window.df <- window.cutoffs.df.20kb.less
  # write.table(sliding.window.df, "../MSSNG.SSC.sliding.window.df.tsv", sep="\t", row.names=F, quote=F, col.names=T)
  
  ## Melt sliding window data frame
  sliding.window.df.melt <- melt(sliding.window.df, 
                                 id.vars = c("size.bin.start","size.bin.end", "bin.count"))

  ## Keep only probands in CH hits data frame
  CH.hits.df.probands <- CH.hits.df[which(CH.hits.df$Relation == "Proband"),]
  CH.hits.df.probands <- CH.hits.df.probands[which(CH.hits.df.probands$exSize <= 30000),]  # remove > 30kb
  CH.hits.df.probands <- CH.hits.df.probands[which(!(CH.hits.df.probands$exSize == 0 & CH.hits.df.probands$CH_hit == 0)),] # remove (0,0)
  
  ## Plot sliding window averages linear model fit - spline
  # png(file="all.cutoff.png", width=900, height=700)
  pdf(file=sprintf("./figures/SPARK.all.cutoff_%s.pdf", file.name), width=9, height=7)
  
  plot(x=sliding.window.df$size.bin.start, y = sliding.window.df$percentile_0.9, col = 1,
       main = sprintf("SPARK Sliding Window Cut-offs (bin size = %s)",file.name), 
       xlab = "Size Bin Start (bp)", ylab = "Control No. SNV",
       cex.main = 1.5, font.main= 1, cex.lab = 1,
       ylim = c(0,7))
  lines(spline(sliding.window.df$size.bin.start, sliding.window.df$percentile_0.9), col = 1)
  points(x=sliding.window.df$size.bin.start, y = sliding.window.df$percentile_0.8, col = 2)
  lines(spline(sliding.window.df$size.bin.start, sliding.window.df$percentile_0.8), col = 2)
  points(x=sliding.window.df$size.bin.start, y = sliding.window.df$median, col = 3)
  lines(spline(sliding.window.df$size.bin.start, sliding.window.df$median), col = 3)
  points(x=sliding.window.df$size.bin.start, y = sliding.window.df$mean, col = 4)
  lines(spline(sliding.window.df$size.bin.start, sliding.window.df$mean), col = 4)
  grid (NULL,NULL, lty = "dotted", col = "grey") 
  legend(1,7, legend=c("90th percentile", "80th percentile", "median", "mean"),
         col= 1:4, lty = 1, cex=1.2, text.font=50,  title="Cut-off")

  dev.off()
}

Get_Sliding_Window_Cutoff_Plots(SPARK.2.5kb.n7.window.cutoffs.df, SPARK.CH.hits, "2.5kb window.n=7")
Get_Sliding_Window_Cutoff_Plots(SPARK.2.5kb.n5.window.cutoffs.df, SPARK.CH.hits, "2.5kb window.n=5")


