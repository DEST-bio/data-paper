### libraries
  library(data.table)
  library(stringr)
  library(ggplot2)
  library(cowplot)
  library(patchwork)
### datasets
  ### orignal bad 36
    good.samps <- fread("/Users/alanbergland/Documents/GitHub/data-paper/SupplementalTable1/identify_problematic_samples/good.samps.txt", header=F)
    all.samps <- fread("/Users/alanbergland/Documents/GitHub/data-paper/SupplementalTable1/identify_problematic_samples/all.samps.txt", header=F)
    bad.samps <- setdiff(all.samps$V1, good.samps$V1)



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

    pnps <- rbind(as.data.table(DATA.poolsnp), as.data.table(DATA.snape), fill=T)[Chrom=="genomewide"][MAC==10 | is.na(MAC)][MAF==0.01]

    ### samps
      samps <- fread("/Users/alanbergland/Documents/GitHub/data-paper/SupplementalTable1/samps.csv")

### a bit of manipulation to get the private SNPs extracted
  setkey(pnps, POP, method)
  setkey(priv.ag, POP, method)
  m <- merge(pnps, priv.ag)
  m <- merge(m, samps, by.x="POP", by.y="sampleId")
  setkey(m, POP)

### plot
  ggplot() +
  geom_point(data=m, aes(x=N, y=pNpS, color=set), size=2) +
  geom_point(data=m[J(bad.samps)], aes(x=N, y=pNpS), shape=5, size=2.5) +
  ylab("pN/pS") +
  xlab("Number of private SNPs (estimate)") +
  theme_bw()
