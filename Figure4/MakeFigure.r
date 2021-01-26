library(tidyverse)
library(gridExtra)
library(scales)

DATA.poolsnp=read.table("~/GitHub/data-paper/Figure4/data/PoolSNP.pnps.gz",
                        header=T,
                        na.strings = "NA")
#summary(DATA.poolsnp)
LIST=as.character(read.table("~/GitHub/data-paper/Figure4/data/exlcude.txt",
                             header=F)[,1])

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
              aes(x=MAC,y=pNpS),
              col=rgb(0,0,0,0.1))+
  geom_ribbon(aes(ymin = pnps.m-pnps.sd, ymax = pnps.m+pnps.sd),
              colour = NA,
              fill=alpha(c("blue"),alpha=0.2))+
  geom_line(lwd=1.5,
            col="blue")+
  xlab("Minor allele count (MAC)")+
  ylab("pN/pS")+
  #geom_abline(slope = 0,intercept=0.4454857483957618,lty=2,col="red",lwd=1.5)+
  theme_bw()+
  theme(axis.title.y = element_text(size = 26, angle = 90),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text=element_text(size=18),
        strip.text =element_text(size=20))+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.05),
                     limits=c(0,2.5))+
  facet_wrap(vars(NAME))
DATA.poolsnp.plot

DATA.poolsnp.group.f <- DATA.poolsnp %>%
  filter(!POP %in% LIST & Chrom =="genomewide" & MAF==0.01) %>%
  group_by(MAC) %>%
  summarise(
    n=n(),
    pnps.m = mean(pNpS),
    pnps.sd=sd(pNpS),
    pnps.se =sd(pNpS)/sqrt(sum(n()))
  )

DATA.poolsnp.group.f$NAME<-rep("PoolSNP",nrow(DATA.poolsnp.group.f))

DATA.poolsnp.plot.f <- ggplot(DATA.poolsnp.group.f,
                              aes(x=MAC,y=pnps.m))+
  geom_jitter(data=filter(DATA.poolsnp,Chrom=="genomewide" & MAF==0.01 & !POP %in% LIST),
              aes(x=MAC,y=pNpS),
              col=rgb(0,0,0,0.1))+
  geom_ribbon(aes(ymin = pnps.m-pnps.sd, ymax = pnps.m+pnps.sd),
              colour = NA,
              fill=alpha(c("blue"),alpha=0.2))+
  geom_line(lwd=1.5,
            col="blue")+
  xlab("Minor allele count (MAC)")+
  ylab("pN/pS")+
  #geom_abline(slope = 0,intercept=0.4454857483957618,lty=2,col="red",lwd=1.5)+
  theme_bw()+
  theme(axis.title.y = element_text(size = 26, angle = 90),
        axis.title.x = element_text(size = 26, angle = 00),
        axis.text=element_text(size=18))+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.05),
                     limits=c(0.2,0.6))
DATA.poolsnp.plot.f

DATA.snape=read.table("~/GitHub/data-paper/Figure4/data/SNAPE.pnps.gz",
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
DATA.snape.group

DATA.snape.group$NAME<-rep("SNAPE",nrow(DATA.snape.group))

null.l=read.table("~/GitHub/data-paper/Figure4/data/null.pnps",header=F)
colnames(null.l)<-c("x","Chrom","NS","SS","y")

NL<-filter(null.l,Chrom == "genomewide")

DATA.snape.plot <- ggplot(DATA.snape.group,
                          aes(x=MAF,y=pnps.m))+
  geom_jitter(data=filter(DATA.snape,Chrom=="genomewide"),
              aes(x=MAF,y=pNpS),
              col=rgb(0,0,0,0.1))+
  geom_ribbon(aes(ymin = pnps.m-pnps.sd, ymax = pnps.m+pnps.sd),
              colour = NA,
              fill=alpha(c("blue"),alpha=0.2))+
  geom_line(lwd=1.5,
            col="blue")+
  xlab("Minor allele frequency (MAF)")+
  ylab("pN/pS")+
  ylim(0,2.5)+
  #geom_line(DATA.poolsnp=NL,aes(x=x,y=y),lty=2,col="red",lwd=1.5)+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text =element_text(size=20))+
  facet_wrap(vars(NAME))
DATA.snape.plot

DATA.snape.group.f <- DATA.snape %>%
  filter(!POP %in% LIST) %>%
  filter(Chrom =="genomewide") %>%
  group_by(MAF,Chrom) %>%
  summarise(
    n=n(),
    pnps.m = mean(pNpS),
    pnps.sd=sd(pNpS),
    pnps.se =sd(pNpS)/sqrt(sum(n()))
  )

DATA.snape.group.f$NAME<-rep("SNAPE",
                             nrow(DATA.snape.group.f))

DATA.snape.plot.f <- ggplot(DATA.snape.group.f,
                            aes(x=MAF,y=pnps.m))+
  geom_jitter(data=filter(DATA.snape,Chrom=="genomewide"& !POP %in% LIST),
              aes(x=MAF,y=pNpS),
              col=rgb(0,0,0,0.1))+
  geom_ribbon(aes(ymin = pnps.m-pnps.sd, ymax = pnps.m+pnps.sd),
              colour = NA,
              fill=alpha(c("blue"),alpha=0.2))+
  geom_line(lwd=1.5,
            col="blue")+
  xlab("Minor allele frequency (MAF)")+
  ylab("pN/pS")+
  ylim(0.2,0.6)+
  #geom_line(DATA.poolsnp=NL,aes(x=x,y=y),lty=2,col="red",lwd=1.5)+
  theme_bw()+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 26, angle = 00),
        axis.text=element_text(size=18))
DATA.snape.plot.f

pdf("~/GitHub/data-paper/Figure4/Figure4.pdf",width=15,height=9)

grid.arrange(DATA.poolsnp.plot,
             DATA.snape.plot,
             DATA.poolsnp.plot.f,
             DATA.snape.plot.f,
             layout_matrix = rbind(c(rep(1,15),rep(2,13)),c(rep(1,15),rep(2,13)),c(rep(1,15),rep(2,13)),c(rep(3,15),rep(4,13)),c(rep(3,15),rep(4,13)))) ## by eyeballing
dev.off()
