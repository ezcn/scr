require(ggplot2) 
#library(readr)
#library(ggrepel)
require(dplyr)
require(RColorBrewer)


gg.manhattan <- function(df, threshold, hlight, col, ylims, title){
#https://github.com/pcgoddard/Burchardlab_Tutorials/wiki/GGplot2-Manhattan-Plot-Function
#Format GWAS results as you would for qqman: SNP CHR BP P (tab-delim)
  # format df
  df.tmp <- df %>% 
    
    # Compute chromosome size
    group_by(CHR) %>% 
    summarise(chr_len=max(BP)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(df, ., by=c("CHR"="CHR")) %>%
    
    # Add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate( BPcum=BP+tot) %>%
    
    # Add highlight and annotation information
    mutate( is_highlight=ifelse(SNP %in% hlight, "yes", "no")) %>%
    mutate( is_annotate=ifelse(Beta < threshold, "yes", "no"))
  
  # get chromosome center positions for x-axis
  axisdf <- df.tmp %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
  ggplot(df.tmp, aes(x=BPcum, y=Beta)) + ########y=-log10(P))) +
    # Show all points
    geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=2) +
    #scale_color_manual(values = rep(col, 22 )) +

    # custom X axis:
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0), limits = ylims) + # expand=c(0,0)removes space between plot area and x axis 
    
    # add plot and axis titles
    ggtitle(paste0(title)) +
    labs(x = "Chromosome") +
    
    # add genome-wide sig and sugg lines
   # geom_hline(yintercept = -log10(sig)) +
    #geom_hline(yintercept = -log10(sugg), linetype="dashed") +
    
    # Add highlighted points
    #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
    
    # Add label using ggrepel to avoid overlapping
  ######## geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(SNP), alpha=0.7), size=5, force=1.3) +
    
    # Custom the theme:
    theme_bw(base_size = 22) +
    theme( 
      plot.title = element_text(hjust = 0.5),
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
}





genome <- data.frame(SNP=character(),
                  CHR=integer(), 
                  BP=integer(), 
		  Beta=numeric(),
                  stringsAsFactors=FALSE)

for (ch in seq(1,22)){ 
    folder="/lustre/scratch115/projects/ukbb500k_t151/selectionBalancing/dataOut/"
    filename=paste(folder, "ukb_impv3_chr", ch, "_04.betaout", sep="") 
    #print (filename)  
    myd=read.table(filename, header=T )
    myd$CHR=ch
    myd$SNP=paste(ch, myd$Position, sep="_")
    mydf<- data.frame(cbind(SNP=myd$SNP, CHR=myd$CHR, BP=myd$Position, Beta=myd$Beta1 ) )
    #print (length(mydf$SNP))
    genome=rbind(genome, mydf) 
    #print (length(genome$CHR))
}

gz1 <- gzfile("betascanGenome.gz", "w")
write.table(genome, file=gz1, quote = FALSE, sep = "\t", row.names = FALSE)


#hlight=c()
#threshold=0
#ylims=c(min(genome$Beta),max(genome$Beta))
#title="BetaScan test"
#gg.manhattan(myd, threshold, hlight, col, ylims, title)
 
