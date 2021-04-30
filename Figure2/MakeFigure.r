library(tidyverse)
library(gridExtra)
library(scales)
library(ggpubr)


## Set global Alpha value
ALPHA=0.15

## set line color
linecol="green"

## get expected values for DGRP data
null.l=read.table("~/Documents/GitHub/data-paper/Figure2/data/null_subset.pnps",header=F)
colnames(null.l)<-c("x","Chrom","NS","SS","y")
NL<-filter(null.l,Chrom == "genomewide")

NL.0 <- filter(NL,x==0)$y

## at first read SNAPE pnps data and filter genomewide values
DATA.snape=read.table("~/Documents/GitHub/data-paper/Figure2/data/SNAPE_full.pnps.gz",
                      header=T,
                      stringsAsFactors = F)
summary(DATA.snape)

DATA.snape.group <- DATA.snape %>%
  filter(Chrom =="genomewide") %>%
  group_by(MAF,Chrom) %>%
  summarise(
    n=n(),
    pnps.m = mean(pNpS),
    pnps.sd=sd(pNpS),
    pnps.se =sd(pNpS)/sqrt(sum(n()))
  )
#DATA.snape.group

## keep pnps without MAF filtering only
DATA.snape.group.MAF0 <- DATA.snape %>%
  filter(MAF == 0 & Chrom =="genomewide")
#DATA.snape.group.MAF0

# read CSV of private SNPs
DATA.ps=read.csv("~/Documents/GitHub/data-paper/Figure2/data/SNAPE_full.ps",
                 header=T,
                 stringsAsFactors = F,
                 sep = "\t")

DATA=merge(DATA.ps,DATA.snape.group.MAF0, by.x="POP",by.y="POP")

DATA$private=log10(DATA$N)

## calculated Mean/SD and threshold based on Mean+2SD for pNpS data
Mean.pNpS=mean(DATA$pNpS)
SD.pNpS=sd(DATA$pNpS)
th.pNpS=Mean.pNpS+1.96*SD.pNpS

##classify
DATA$TH.pNpS <-DATA$pNpS
DATA$TH.pNpS[DATA$TH.pNpS<th.pNpS]<-NA
DATA$TH.pNpS[DATA$TH.pNpS>=th.pNpS]<-"Exclude"
DATA$TH.pNpS[is.na(DATA$TH.pNpS)]<-"Keep"

## calculated Mean/SD and threshold based on Mean+1.96SD for private SNP data
Mean.private=mean(DATA$private)
SD.private=sd(DATA$private)
th.private=Mean.private+1.96*SD.private

#classify
DATA$TH.private <-DATA$private
DATA$TH.private[DATA$TH.private<th.private]<-NA
DATA$TH.private[DATA$TH.private>=th.private]<-"Exclude"
DATA$TH.private[is.na(DATA$TH.private)]<-"Keep"

# make new final classification, where any pop will be excluded that is excluded based on either measurement (pNpS and/or private SNPs)
DATA$Status <- rep("Keep",length(DATA$TH.private))
DATA$Status[DATA$TH.pNpS=="Exclude"]<-"Exclude"
DATA$Status[DATA$TH.private=="Exclude"]<-"Exclude"

write.table(DATA,"~/Documents/GitHub/data-paper/Figure2/data/classify_pops.txt",quote = F,row.names = F)

DATA.RAL<-DATA %>%
  filter(POP=="NC_ra_03_n")

Classify.plot<-ggplot(DATA,aes(x=pNpS,y=private,col=Status))+
  geom_point(size=10)+
  xlab(expression(italic("p")["N"]/italic("p")["S"]))+
  ylab("No. of private SNPs")+
  theme_bw()+
  geom_hline(yintercept=th.private,lty=2,col="black",lwd=1)+
  geom_vline(xintercept=th.pNpS,lty=2,col="black",lwd=1)+
  geom_vline(xintercept=NL.0,lty=2,col=linecol,lwd=1)+
  geom_point(data=DATA.RAL,
             aes(x=pNpS,y=private),col=linecol,shape=18,size=10)+
  theme(axis.title.y = element_text(size = 20, angle = 90),
        axis.title.x = element_text(size = 20, angle = 0),
        axis.text=element_text(size=18),
        strip.text =element_text(size=20),
        legend.position = c(0.88, 0.88),
        legend.title = element_text(color = "black", size = 20),
        legend.text=element_text(size=16))+
  scale_fill_manual(values=c(alpha(c("red"),
                                   alpha=ALPHA),
                             alpha(c("black"),
                                   alpha=ALPHA)))+
  scale_colour_manual(values=c(alpha(c("red"),
                                     alpha=ALPHA),
                               alpha(c("black"),
                                     alpha=ALPHA)))

# pdf("~/Desktop/classify.pdf",width=10,height=8)
# Classify.plot
# dev.off()

LIST<-DATA$POP[DATA$Status=="Exclude"]
STATUSLIST=data.frame(POP=DATA$POP,Status=DATA$Status)

DATA.poolsnp=read.table("~/Documents/GitHub/data-paper/Figure2/data/PoolSNP.pnps.gz",
                        header=T,
                        na.strings = "NA")

## keep pnps without MAF filtering only
DATA.poolsnp.group.MAF0 <- DATA.poolsnp %>%
  filter(MAF == 0.01 & MAC==10 & MAFc ==0 & Chrom =="genomewide")
#DATA.snape.group.MAF0

DATA.group <- merge(DATA.poolsnp.group.MAF0,DATA.snape.group.MAF0,by.x="POP",by.y="POP")
plot(DATA.group$pNpS.x,DATA.group$pNpS.y)

DATA.poolsnp.RAL<-DATA.poolsnp %>%
  filter(POP=="NC_ra_03_n")

DATA.poolsnp<-merge(DATA.poolsnp,STATUSLIST,by.x="POP",by.y="POP")

DATA.poolsnp.group <- DATA.poolsnp %>%
  filter(Chrom =="genomewide" & MAF==0.01) %>%
  group_by(MAC) %>%
  summarise(
    n=n(),
    pnps.m = mean(pNpS),
    pnps.sd=sd(pNpS),
    pnps.se =sd(pNpS)/sqrt(sum(n()))
  )

DATA.poolsnp.group$NAME<-rep("PoolSNP",
                             nrow(DATA.poolsnp.group))

DATA.poolsnp.plot <- ggplot(DATA.poolsnp.group,
                            aes(x=MAC,y=pnps.m))+
  geom_jitter(data=filter(DATA.poolsnp,Chrom=="genomewide" & MAF==0.01),
              aes(x=MAC,y=pNpS,col=Status))+
  geom_ribbon(aes(ymin = pnps.m-pnps.sd, ymax = pnps.m+pnps.sd),
              colour = NA,
              fill=alpha(c("blue"),alpha=0.2))+
  geom_line(lwd=1.5,
            col="blue")+
  geom_point(data=filter(DATA.poolsnp.RAL,Chrom=="genomewide" & MAF==0.01),
             aes(x=MAC,y=pNpS),col=linecol,shape=18,size=3)+
  xlab("Minor allele count (MAC)")+
  ylab(expression(italic("p")["N"]/italic("p")["S"]))+
  geom_abline(slope = 0,intercept=NL.0,lty=2,col=linecol,lwd=1)+
  theme_bw()+
  theme(axis.title.y = element_text(size = 20, angle = 90),
        axis.title.x = element_text(size = 20, angle = 0),
        axis.text=element_text(size=18),
        strip.text =element_text(size=20),
        legend.title = element_blank(),
        legend.position = "none")+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.05),
                     limits=c(0,2.5))+
  scale_fill_manual(values=c(alpha(c("red"),
                                   alpha=ALPHA),
                             alpha(c("black"),
                                   alpha=ALPHA)))+
  scale_colour_manual(values=c(alpha(c("red"),
                                     alpha=ALPHA),
                               alpha(c("black"),
                                     alpha=ALPHA)))+
  facet_wrap(vars(NAME))
DATA.poolsnp.plot


DATA.snape<-merge(DATA.snape,STATUSLIST,by.x="POP",by.y="POP")
DATA.snape.group$NAME<-rep("SNAPE",nrow(DATA.snape.group))

DATA.snape.RAL<-DATA.snape %>%
  filter(POP=="NC_ra_03_n" & Chrom =="genomewide")


DATA.snape.plot <- ggplot(DATA.snape.group,
                          aes(x=MAF,y=pnps.m))+
  geom_jitter(data=filter(DATA.snape,Chrom=="genomewide"),
              aes(x=MAF,y=pNpS,col=Status))+
  geom_ribbon(aes(ymin = pnps.m-pnps.sd, ymax = pnps.m+pnps.sd),
              colour = NA,
              fill=alpha(c("blue"),alpha=0.2))+
  geom_point(data=DATA.snape.RAL,
             aes(x=MAF,y=pNpS),col=linecol,shape=18,size=3)+
  geom_line(lwd=1.5,
            col="blue")+
  xlab("Minor allele frequency (MAF)")+
  ylab(expression(italic("p")["N"]/italic("p")["S"]))+
  ylim(0,2.5)+
  geom_line(data=NL,aes(x=x,y=y),lty=2,col=linecol,lwd=1)+
  theme_bw()+
  theme(axis.title.y = element_text(size = 20, angle = 90),
        axis.title.x = element_text(size = 20, angle = 0),
        axis.text=element_text(size=18),
        strip.text =element_text(size=20),
        legend.title = element_blank(),
        legend.position = "none")+
  scale_fill_manual(values=c(alpha(c("red"),
                                   alpha=ALPHA),
                             alpha(c("black"),
                                   alpha=ALPHA)))+
  scale_colour_manual(values=c(alpha(c("red"),
                                     alpha=ALPHA),
                               alpha(c("black"),
                                     alpha=ALPHA)))+
  facet_wrap(vars(NAME))
DATA.snape.plot


FP<-ggarrange(Classify.plot,
              ggarrange(DATA.poolsnp.plot,
                        DATA.snape.plot,
                        labels = c("B","C"),
                        font.label = list(size = 28, face = "bold"),
                        nrow=2),
              ncol=2,
              labels="A",
              common.legend = T,
              legend = "bottom",
              font.label = list(size = 28, face = "bold"))

ggsave("~/Documents/GitHub/data-paper/Figure2/Figure2.pdf",
       FP,
       device="pdf",
       width=15,
       height=9)
