library(ggplot2)

datapop2L <- read.table("all_pops_2L_header_PoolSNP_SNAPE.stats.txt", h=F)
datapop2R <- read.table("all_pops_2R_header_PoolSNP_SNAPE.stats.txt", h=F)
datapop3L <- read.table("all_pops_3L_header_PoolSNP_SNAPE.stats.txt", h=F)
datapop3R <- read.table("all_pops_3R_header_PoolSNP_SNAPE.stats.txt", h=F)
datapop4 <- read.table("all_pops_4_header_PoolSNP_SNAPE.stats.txt", h=F)
datapopY <- read.table("all_pops_Y_header_PoolSNP_SNAPE.stats.txt", h=F)
datapopX <- read.table("all_pops_X_header_PoolSNP_SNAPE.stats.txt", h=F)

datachrom <- datapop2L
# Common polymorphic positions
datachrom$V5 <- datapop2L$V5 + datapop2R$V5 + datapop3L$V5 + datapop3R$V5 + datapop4$V5 + datapopY$V5 + datapopX$V5
# Polymorphic positions with same frequency
datachrom$V6 <- datapop2L$V6 + datapop2R$V6 + datapop3L$V6 + datapop3R$V6 + datapop4$V6 + datapopY$V6 + datapopX$V6
# Frequency difference less than 0.05
datachrom$V7 <- datapop2L$V7 + datapop2R$V7 + datapop3L$V7 + datapop3R$V7 + datapop4$V7 + datapopY$V7 + datapopX$V7
# Frequency difference less than 0.1
datachrom$V8 <- datapop2L$V8 + datapop2R$V8 + datapop3L$V8 + datapop3R$V8 + datapop4$V8 + datapopY$V8 + datapopX$V8
# Frequency difference higher than 0.1
datachrom$V9 <- datapop2L$V9 + datapop2R$V9 + datapop3L$V9 + datapop3R$V9 + datapop4$V9 + datapopY$V9 + datapopX$V9
# Polymorphic position in PoolSNP and missing data in SNAPE
datachrom$V13 <- datapop2L$V13 + datapop2R$V13 + datapop3L$V13 + datapop3R$V13 + datapop4$V13 + datapopY$V13 + datapopX$V13
# Polymorphic position in PoolSNP and monomorphic position in SNAPE
datachrom$V14 <- datapop2L$V14 + datapop2R$V14 + datapop3L$V14 + datapop3R$V14 + datapop4$V14 + datapopY$V14 + datapopX$V14
# Polymorphic position in SNAPE and missing data in PoolSNP
datachrom$V15 <- datapop2L$V15 + datapop2R$V15 + datapop3L$V15 + datapop3R$V15 + datapop4$V15 + datapopY$V15 + datapopX$V15
# Polymorphic position in SNAPE and monomorphic position in PoolSNP
datachrom$V16 <- datapop2L$V16 + datapop2R$V16 + datapop3L$V16 + datapop3R$V16 + datapop4$V16 + datapopY$V16 + datapopX$V16
# Monomorphic position in PoolSNP and missing data in SNAPE
datachrom$V17 <- datapop2L$V17 + datapop2R$V17 + datapop3L$V17 + datapop3R$V17 + datapop4$V17 + datapopY$V17 + datapopX$V17
# Monomorphic position in SNAPE and missing data in PoolSNP
datachrom$V18 <- datapop2L$V18 + datapop2R$V18 + datapop3L$V18 + datapop3R$V18 + datapop4$V18 + datapopY$V18 + datapopX$V18
# Total differences between both callers
datachrom$V19 <- datapop2L$V19 + datapop2R$V19 + datapop3L$V19 + datapop3R$V19 + datapop4$V19 + datapopY$V19 + datapopX$V19
# Total positions analyzed
datachrom$V20 <- datapop2L$V20 + datapop2R$V20 + datapop3L$V20 + datapop3R$V20 + datapop4$V20 + datapopY$V20 + datapopX$V20


datanames <- read.table("names_246pops.txt")
datachrom$V21 <- datanames[,1]

# Delete 22 problematic samples
datachrom <- datachrom[-c(1,2,5,28,73,75,80,81,102,114,140,151,152,167,168,185,186,191,210,232), ]

write.table(datachrom, file = "caller_comparison_226pop.txt", sep = "\t")



###################################################################
### PLOTTING
###################################################################


### FIRST PLOT
datafreq2chrom <- data.frame(
  type = c(rep("Same Freq", length(datachrom$V1)), 
           rep("Freq diff =< 0.05", length(datachrom$V1)),
           rep("Freq diff =< 0.1", length(datachrom$V1)),
           rep("Freq diff > 0.1", length(datachrom$V1))),
  value = c((datachrom$V6/datachrom$V5), 
            ((datachrom$V7)/datachrom$V5),
            ((datachrom$V8)/datachrom$V5),
            ((datachrom$V9)/datachrom$V5))
)

datafreq2chrom$names <-  datachrom$V21

freq_plot_chrom <- ggplot(data=datafreq2chrom, aes(x=names, y=value, fill=type, order=type)) + 
  geom_bar(stat="identity") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  #geom_text(aes(label=round(value, digits=2)), vjust=1.6, size=1.2, col = "red") + 
  labs(title="PoolSNP and SNAPE polymorphic calls",x="Population", y = "Percentage")

### SECOND PLOT

diff_plot_chrom <- ggplot(data=datachrom, aes(x=V21, y=(V19/V20))) + 
  geom_bar(stat="identity") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) + 
  #geom_text(aes(label=round(value, digits=2)), vjust=1.6, size=1.2, col = "red") + 
  ylim(0,1) +
  labs(title="PoolSNP and SNAPE different calls",x="Population", y = "Percentage")

### THIRD PLOT
datadiffchrom <- data.frame(
  type = c(rep("Polymorphic in PoolSNP and missing data in SNAPE", length(datachrom$V1)), 
           rep("Polymorphic in PoolSNP and monomorphic in SNAPE", length(datachrom$V1)),
           rep("Polymorphic in SNAPE and monomorphic in PoolSNP", length(datachrom$V1)),
           rep("Polymorphic in SNAPE and missing data in PoolSNP", length(datachrom$V1)),
           rep("Monomorphic in PoolSNP and missing data in PoolSNP", length(datachrom$V1)),
           rep("Monomorphic in SNAPE and missing data in PoolSNP", length(datachrom$V1))),
  value = c((datachrom$V13/datachrom$V19), 
            ((datachrom$V14)/datachrom$V19),
            ((datachrom$V15)/datachrom$V19),
            ((datachrom$V16)/datachrom$V19),
            ((datachrom$V17)/datachrom$V19),
            ((datachrom$V18)/datachrom$V19))
)

datadiffchrom$names <-  datachrom$V21

diff_type_chrom_plot <- ggplot(data=datadiffchrom, aes(x=names, y=value, fill=type, order=type)) + 
  geom_bar(stat="identity") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  #geom_text(aes(label=round(value, digits=2)), vjust=1.6, size=1.2, col = "red") + 
  labs(title="PoolSNP and SNAPE different calls differences",x="Population", y = "Percentage")


### ALL PLOTS TOGETHER
library(cowplot)
plot_grid(freq_plot_chrom + theme(legend.position = "bottom", legend.title = element_blank()), diff_plot_chrom, diff_type_chrom_plot + theme(legend.position ="bottom", legend.title = element_blank()), labels = c('A', 'B', 'C'), label_size = 12, nrow=3)

