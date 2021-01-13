# plot SupFigure S5 - Margot
library(ggplot2)
library("gridExtra")

dat<-read.table(file="SupFigureS5_data.txt",header=T,fill =T,sep="\t")

plot1<-ggplot(dat, aes(x=CHROM, y=Pi)) + geom_boxplot(, show.legend = FALSE) + labs(y = "Nucleotide diversity") + theme_bw() + theme(axis.title.x=element_blank()) + facet_grid(~ DATASET) + scale_x_discrete(limits=c("X","2L","2R","3L","3R"))
plot2<-ggplot(dat, aes(x=CHROM, y=Watterson)) + geom_boxplot(, show.legend = FALSE) + labs(y = "Theta Watterson") + theme_bw()  + theme(axis.title.x=element_blank()) + facet_grid(~ DATASET) + scale_x_discrete(limits=c("X","2L","2R","3L","3R"))
plot3<-ggplot(dat, aes(x=CHROM, y=Tajima_D)) + geom_boxplot(, show.legend = FALSE) +  labs(y = "Tajima's D") +  theme_bw() + theme(axis.title.x=element_blank()) + facet_grid(~ DATASET) + scale_y_continuous(breaks=seq(-1, 2, 0.5)) + scale_x_discrete(limits=c("X","2L","2R","3L","3R")) 
plot4<-ggplot(dat, aes(x=CHROM, y=pn_ps)) + geom_boxplot(, show.legend = FALSE) + labs(y =expression('p'[N]*'/p'[S])) +  theme_bw() + theme(axis.title.x=element_blank()) + labs(x = "Chromosome") + facet_grid(~ DATASET) + scale_x_discrete(limits=c("X","2L","2R","3L","3R"))

# make final SupFigure S5
library("gridExtra")
grid.arrange(plot1, plot2, plot3, plot4, 
        ncol=1, nrow=4, heights=c(1,1,1,1))




