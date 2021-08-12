
### module load gcc/7.1.0 openmpi/3.1.4 R/3.6.3; R


\
### load data
  fs <- list.files("/project/berglandlab/moments/moments_output", full.names=T)

### read
  o <- foreach(i=fs, .errorhandling="remove")%do%{
    #i<-fs[254]
    tmp <- fread(i)
  }
  o <- rbindlist(o)
  setnames(o, "-2LL_model", "LL")
  o <- o[Pair_name!="Pair_name"]
  o <- na.omit(o)
  o.ag <- o[,list(divergence_time=divergence_time[which.min(AIC)],
                  theta=theta[which.min(AIC)],
                  pop1_size=pop1_size[which.min(AIC)],
                  pop2_size=pop2_size[which.min(AIC)]),
            list(Pair_name)]

  o.ag[,caller:=tstrsplit(Pair_name, "\\.")[[3]]]
  o.ag[,pair:=paste(tstrsplit(Pair_name, "\\.")[[1]], tstrsplit(Pair_name, "\\.")[[2]], sep=".")]

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

  cor.test(o_wide$theta_PoolSNP , o_wide$theta_SNAPE)
  cor.test(o_wide$pop1_size_PoolSNP , o_wide$pop1_size_SNAPE)
  cor.test(o_wide$pop2_size_PoolSNP, o_wide$pop2_size_SNAPE)

  cor.test(o_wide$divergence_time_PoolSNP,  o_wide$divergence_time_SNAPE)

  save(o, o.ag, o_wide, file="~/moments_o.Rdata")
