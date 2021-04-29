### libraries
  library(data.table)
  library(stringr)
  library(ggplot2)
  library(cowplot)
  library(patchwork)

### datasets
  ### Maria's list
    mList <- c("AT_Mau_14_01","AT_Mau_14_02","AT_See_14_44",
    "CH_Cha_14_42","CH_Cha_14_43","DE_Mun_14_32",
    "DK_Kar_14_39","DK_Kar_14_41","ES_Gim_14_34",
    "ES_Gim_14_35","FI_Aka_14_36","FI_Aka_14_37",
    "FI_Ves_14_38","FL_ho_10_fall","FR_Got_14_08",
    "FR_Vil_14_05","FR_Vil_14_07","ME_bo_09_fall.r1",
    "ME_bo_09_fall.r2","ON_su_15_spring","PA_li_09_spring",
    "PA_li_10_spring","PA_li_11_fall","PA_li_15_spring","PT_Rec_14_33",
    "RU_Val_14_51","TR_Yes_14_03","TR_Yes_14_04","UA_Che_14_47",
    "UA_Che_14_48","UA_Kha_14_45","UA_Kha_14_46","UA_Kyi_14_49","UA_Uma_14_50","UK_She_14_09","UK_Sou_14_10")

  ### bad samps list as implemented
    good.samps <- fread("/Users/alanbergland/Documents/GitHub/data-paper/SupplementalTable1/identify_problematic_samples/good.samps.txt", header=F)
    all.samps <- fread("/Users/alanbergland/Documents/GitHub/data-paper/SupplementalTable1/identify_problematic_samples/all.samps.txt", header=F)
    bad.samps <- setdiff(all.samps$V1, good.samps$V1)

    setdiff(mList, bad.samps)

  ## for private snps
    ### scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/priv_ag.Rdata ~/.

    load("~/priv_ag.Rdata")
    priv.ag[,method:="snape"]

  ## for pn/ps
    DATA.poolsnp=read.table("~/Documents/GitHub/data-paper/Figure2/data/PoolSNP.pnps.gz",
                            header=T,
                            na.strings = "NA")
    DATA.poolsnp$method="poolsnp"

    DATA.snape=read.table("~/Documents/GitHub/data-paper/Figure2/data/SNAPE.pnps.gz",
                          header=T,
                          stringsAsFactors = F)
    DATA.snape$method="snape"

    pnps <- rbind(as.data.table(DATA.poolsnp), as.data.table(DATA.snape), fill=T)[Chrom=="genomewide"][MAC==10 | is.na(MAC)][MAF==0.05]

    ### samps
      samps <- fread("/Users/alanbergland/Documents/GitHub/data-paper/SupplementalTable1/samps.csv")

### read depth info
  rd <- fread("/Users/alanbergland/Documents/GitHub/DEST_freeze1/populationInfo/sequencingStats/rd.csv")

### tajima's D
  ss <- fread("/Users/alanbergland/Documents/GitHub/data-paper/Figure6/Figure6_data.txt")
  ss <- ss[FILE=="PoolSNP"]

### a bit of manipulation to get the private SNPs extracted
  setkey(pnps, POP, method)
  setkey(priv.ag, POP, method)
  m <- merge(pnps, priv.ag)
  m <- merge(m, samps, by.x="POP", by.y="sampleId")
  setkey(m, POP)
  m[J(bad.samps), bad:=T]
  m[is.na(bad), bad:=F]
  m <- merge(m, rd[auto==T], by.x="POP", by.y="sampleId")
  m <- merge(m, ss)

### plot

  p1 <- ggplot() +
  geom_point(data=m, aes(x=log10(N), y=pNpS, color=set), size=2, alpha=.8) +
  geom_point(data=m[J(bad.samps)], aes(x=log10(N), y=pNpS), shape=5, size=2.5, stroke=1.25) +
  ylab("pN/pS") +
  xlab("log10(Number of private SNPs)\n(estimate)") +
  theme_bw()

  p2 <- ggplot() +
  geom_point(data=m, aes(x=log10(N), y=Tajima_D, color=set), size=2, alpha=.8) +
  geom_point(data=m[J(bad.samps)], aes(x=log10(N), y=Tajima_D), shape=5, size=2.5, stroke=1.25) +
  ylab("Tajima_D (estimated from PoolSNP)") +
  xlab("log10(Number of private SNPs)\n(estimate)") +
  theme_bw()

  p3 <- ggplot() +
  geom_point(data=m, aes(x=pNpS, y=Tajima_D, color=set), size=2, alpha=.8) +
  geom_point(data=m[J(bad.samps)], aes(x=pNpS, y=Tajima_D), shape=5, size=2.5, stroke=1.25) +
  ylab("Tajima_D (estimated from PoolSNP)") +
  xlab("pn/ps") +
  theme_bw()

  p4 <- ggplot() +
  geom_point(data=m, aes(x=mu.25, y=Tajima_D, color=set), size=2, alpha=.8) +
  geom_point(data=m[J(bad.samps)], aes(x=mu.25, y=Tajima_D), shape=5, size=2.5, stroke=1.25) +
  ylab("Tajima_D (estimated from PoolSNP)") +
  xlab("rd") +
  theme_bw()


  ggsave(p1 + p2 + p3 + p4, file="~/priv_pnps.png", width=10, height=7)
