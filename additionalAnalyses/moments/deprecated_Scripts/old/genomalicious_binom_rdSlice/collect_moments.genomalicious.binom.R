
### module load gcc/7.1.0 openmpi/3.1.4 R/3.6.3; R


### library
  library(data.table)
  library(foreach)
  library(sp)

### load data
  fs <- list.files("/project/berglandlab/moments/moments_output_genomalicious_binom", full.names=T)

### read
  o <- foreach(i=fs, .errorhandling="remove")%do%{
    #i<-fs[254]
    tmp <- fread(i)
  }
  o <- rbindlist(o)
  setnames(o, "-2LL_model", "LL")
  o <- o[Pair_name!="Pair_name"]
  o.ag <- o[,list(divergence_time=divergence_time[which.min(na.omit(AIC))],
                  theta=theta[which.min(na.omit(AIC))],
                  pop1_size=pop1_size[which.min(na.omit(AIC))],
                  pop2_size=pop2_size[which.min(na.omit(AIC))]),
            list(Pair_name)]

  o.ag[,caller:=tstrsplit(Pair_name, "\\.")[[2]]]
  o.ag[,pair:=paste(tstrsplit(Pair_name, "\\.")[[3]], tstrsplit(Pair_name, "\\.")[[4]], sep=".")]

### tack in some meta data

  ### load samps
    setwd("/scratch/aob2x/")
    samps <- fread("DEST_freeze1/populationInfo/samps_10Nov2020.csv")

  ### Get E/W cluster IDs
    clusters <- fread("DEST_freeze1/populationInfo/Cluster_Assingment/DEST_Sample_clusters.txt")
    samps <- merge(samps, clusters[,c("sampleId", "Continental_clusters"), with=F], by="sampleId")

  ### make full pairwise table
    pairs <- CJ(samps[Continental_clusters=="1.Europe_W" & set=="DrosEU"]$sampleId,
                samps[Continental_clusters=="3.Europe_E" & set=="DrosEU"]$sampleId)
    pairs[,popset:="EW"]

    pairs_E <- CJ(samps[Continental_clusters=="3.Europe_E" & set=="DrosEU"]$sampleId,
                samps[Continental_clusters=="3.Europe_E" & set=="DrosEU"]$sampleId)
    pairs_E[,popset:="EE"]

    pairs_W <- CJ(samps[Continental_clusters=="1.Europe_W" & set=="DrosEU"]$sampleId,
                samps[Continental_clusters=="1.Europe_W" & set=="DrosEU"]$sampleId)
    pairs_W[,popset:="WW"]

    pairs <- rbindlist(list(pairs, pairs_E, pairs_W))

  ### attach location information
    setnames(pairs, "V1", "sampleId")
    pairs <- merge(samps, pairs, by="sampleId")[,c("sampleId", "V2", "lat", "long", "popset"), with=F]
    setnames(pairs, c("sampleId", "V2", "lat", "long"), c("V1", "sampleId", "lat.V1", "long.V1"))

    pairs <- merge(samps, pairs, by="sampleId")[,c("V1", "sampleId", "lat.V1", "long.V1", "lat", "long", "popset"), with=F]
    setnames(pairs, c("sampleId", "lat", "long"), c("V2", "lat.V2", "long.V2"))

    pairs[,dist:=spDists(x=as.matrix(pairs[,c("long.V1", "lat.V1"), with=F]),
                          y=as.matrix(pairs[,c("long.V2", "lat.V2"), with=F]), diagonal=T)]

    pairs[,id:=as.character(1:dim(pairs)[1])]
    pairs[,pair:=paste(V1, V2, sep=".")]

### merge
  o_wide <- dcast(o.ag, pair ~ caller, value.var=c("divergence_time", "theta", "pop1_size", "pop2_size"))

  o_wide <- merge(o_wide, pairs, by="pair")
  o_wide
  o_wide[,dt:=as.numeric(as.character(divergence_time_PoolSNP))]

  summary(lm(divergence_time_PoolSNP~as.numeric(as.character(dist)), o_wide))
  summary(lm(divergence_time_PoolSNP~popset, o_wide))

  cor.test(o_wide$theta_PoolSNP , o_wide$theta_SNAPE)
  cor.test(o_wide$pop1_size_PoolSNP , o_wide$pop1_size_SNAPE)
  cor.test(o_wide$pop2_size_PoolSNP, o_wide$pop2_size_SNAPE)

  cor.test(o_wide$divergence_time_PoolSNP,  o_wide$divergence_time_SNAPE)
)
  save(o, o.ag, o_wide, file="~/moments_o_genomalicious_binomial.Rdata")


  # scp aob2x@rivanna.hpc.virginia.edu:~/moments_o.Rdata ~/.
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(patchwork)

  load("~/moments_o.Rdata")
  o_wide[,divergence_time:=as.numeric(as.character(divergence_time_PoolSNP))]
  o.ag[,theta:=theta/83960116]
  oo <- melt(o.ag, id.vars=c("pair", "caller", "Pair_name"))
  oow <- dcast(oo, variable+pair~caller, value.var="value")


  ps2 <- oow[variable=="pop2_size" & SNAPE > 1e6]
  ps1 <- oow[variable=="pop1_size" & SNAPE > 1e6]
  theta <- oow[variable=="theta" & PoolSNP > 6e5]
  dt <- oow[variable=="divergence_time" & PoolSNP > 5e4]

  flags <- rbindlist(list(ps2, ps1, theta, dt))

tab <- table(unlist(tstrsplit(flags$pair, "\\.")))
tab <- data.table(sampleId=names(tab), N=as.numeric(tab))
tab[,x:=rank(-N, ties="random")]

flags[,use:=apply(foreach(pi=tab[N>4]$sampleId, .combine="cbind")%do%grepl(pi, flags$pair),
                  1, any)]
setkey(oow, pair)
setkey(flags, pair)
oow[J(flags),outlier:=pair]

oow[,use:=apply(foreach(pi=tab[N>5]$sampleId, .combine="cbind")%do%ifelse(grepl(pi, oow$pair), pi, ""),
                  1, function(x) paste(x[x!=""], collapse="."))]
oow[use=="", use:="others"]

options(digits = 2, scipen = -2)
p1 <-
ggplot(data=oow, aes(x=(SNAPE), y=(PoolSNP), color=use)) +
geom_point() +
facet_wrap(~variable, scales="free") + geom_abline(slope=1, intercept=0) +
#geom_text_repel(data=flags[use==T], aes(label=pair), size=2, force=2) +
theme_bw() +
theme(legend.position="bottom")

p2 <- ggplot(data=tab, aes(x=x, y=N, fill=(N>5))) + geom_bar(stat="identity") +
  geom_text(aes(y=N+5, label=sampleId), angle=0, size=2) +
  ylim(0,25) + theme_bw() + coord_flip() +
  theme(legend.position="none") +
  xlab("Rank")

layout <- "
AAAAABB
AAAAABB
AAAAABB
"


mega <- p1 + p2 +
plot_layout(design = layout) +
plot_annotation(tag_levels = 'A')


ggsave(mega, height=5, w=7, file="~/moments.pdf")
