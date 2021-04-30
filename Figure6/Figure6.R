# plot Figure6 - Margot

library(ggplot2)
library("gridExtra")
library(scales)

dat<-read.table(file="/Users/martinkapun/Documents/GitHub/data-paper/Figure6/Figure6_data.txt",header=T,fill =T,sep="\t")

dat_PoolSNP<-dat[dat$FILE=="PoolSNP",]
dat_SNAPE<-dat[dat$FILE=="SNAPE",]

# remove 20 pops to be excluded from SNAPE 
Meta<-read.table(file="~/Documents/GitHub/data-paper/Figure2/data/classify_pops.txt",header=T)
exclude <- filter(Meta,Status=="Exclude")
dat_SNAPE.sub <-filter(dat_SNAPE, !(POP %in% exclude$POP))

dat2 <-data.frame(rbind(dat_PoolSNP,dat_SNAPE.sub))
dat2 <-filter(dat2,Continent %in% c("Africa","North America","Europe"))

# PoolSNP plots

plot1<-ggplot(dat2, aes(x=Continent, y=Pi)) +
  geom_boxplot(, show.legend = FALSE) +
  labs(y = expression(paste("Nucleotide diversity (",italic(pi),")")), tag = "A") +
  theme_bw() +
  theme(plot.tag = element_text(size = 22, angle = 0),
        axis.title.y = element_text(size = 18, angle = 90),
        axis.title.x=element_blank(),
        axis.text.x = element_text(size = 12, angle = 00),
        axis.text.y = element_text(size = 12, angle = 00)) +
  facet_grid(.~FILE,space="free") +
  scale_x_discrete(limits=c("Africa","North America","Europe"),
                   labels = wrap_format(10)) +
  scale_y_continuous(labels=function(x){
    sprintf("%.3f", x)},
    limits=c(0.0037, 0.007),
    breaks=seq(0.004, 0.007, 0.001))  

plot2<-ggplot(dat2, aes(x=Continent, y=Watterson)) +
  geom_boxplot(, show.legend = FALSE) +
  labs(y = expression(paste("Watterson's ", italic(theta))), tag = "C") +
  theme_bw()  +
  theme(plot.tag = element_text(size = 22, angle = 0),
        axis.title.x=element_blank(),
        axis.title.y = element_text(size = 18, angle = 90),
        axis.text.x = element_text(size = 12, angle = 00),
        axis.text.y = element_text(size = 12, angle = 0)) +
  facet_grid(~ FILE, 
             scales="free") +
  scale_x_discrete(limits=c("Africa","North America","Europe"),
                   labels = wrap_format(10)) +
  scale_y_continuous(labels=function(x){
    sprintf("%.3f", x)},
    limits=c(0.0026, 0.0068),
    breaks=seq(0.003, 0.006, 0.001))  

plot3<-ggplot(dat2, aes(x=Continent, y=Tajima_D)) +
  geom_boxplot(, show.legend = FALSE) + 
  labs(y = expression(paste("Tajima's ",italic(D))), tag = "B") + 
  theme_bw() +
  theme(plot.tag = element_text(size = 22, angle = 0),
        axis.title.x=element_blank(),
        axis.title.y = element_text(size = 18, angle = 90),
        axis.text.x = element_text(size = 12, angle = 00),
        axis.text.y = element_text(size = 12, angle = 00)) +
  facet_wrap(~FILE, scales = "free_x") +
  scale_x_discrete(limits=c("Africa","North America","Europe"),
                   labels = wrap_format(10)) +
  scale_y_continuous(labels=function(x){
    sprintf("%.3f", x)},
    limits=c(-0.285, 1.65),
    breaks=seq(0, 1.5, 0.5)) 

pdf("~/Documents/GitHub/data-paper/Figure6/figure/Figure6.pdf",widt=10,height=6)
# make final Figure 6
grid.arrange(plot1, 
             plot3, 
             plot2, 
             layout_matrix = rbind(c(1,2),
                                   c(1,3)))
dev.off()
        

