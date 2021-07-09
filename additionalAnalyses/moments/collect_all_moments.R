### module load gcc/7.1.0 openmpi/3.1.4 R/3.6.3; R

### libraries
  library(data.table)
  library(SeqArray)
  library(foreach)
  library(sp)
  library(doMC)
  registerDoMC(8)

### get files
  fs <- list.files("/scratch/aob2x/moments_general/output", full.names=T)

  o <- foreach(i=fs, .errorhandling="remove")%dopar%{
    print(which(i==fs))
    #i<-fs[254]
    tmp <- fread(i)
  }
  o <- rbindlist(o)
  setnames(o, "-2LL_model", "LL")
  o <- o[Pair_name!="Pair_name"]
  o <- na.omit(o)
  o[,AIC:=as.numeric(AIC)]
  o.ag <- o[,list(divergence_time=as.numeric(divergence_time[which.min((AIC))]),
                  theta=as.numeric(theta[which.min((AIC))])/as.numeric(L[which.min((AIC))]),
                  pop1_size=as.numeric(pop1_size[which.min((AIC))]),
                  pop2_size=as.numeric(pop2_size[which.min((AIC))]),
                  pop1_mig=as.numeric(mig_pop1[which.min((AIC))]),
                  pop2_mig=as.numeric(mig_pop2[which.min((AIC))]),
                  .N),
            list(Pair_name)]

  o.ag[,SNP_caller:=tstrsplit(Pair_name, "\\.")[[1]]]
  o.ag[,SFS_method:=tstrsplit(Pair_name, "\\.")[[2]]]
  o.ag[,RD_filter:=tstrsplit(Pair_name, "\\.")[[3]]]
  o.ag[,pair:=paste(tstrsplit(Pair_name, "\\.")[[4]], tstrsplit(Pair_name, "\\.")[[5]], sep=".")]

### tack in some meta data

  ### load samps
    setwd("/scratch/aob2x/")
    samps <- fread("DEST_freeze1/populationInfo/samps_10Nov2020.csv")

  ### Get cluster combos IDs
    clusters <- fread("DEST_freeze1/populationInfo/Cluster_Assingment/DEST_Sample_clusters.txt")
    samps <- merge(samps, clusters[,c("sampleId", "Continental_clusters"), with=F], by="sampleId")

    o.ag[,sampleId:=tstrsplit(Pair_name, "\\.")[[4]]]
    o.ag <- merge(o.ag, clusters[,c("sampleId", "Continental_clusters"), with=F], by="sampleId")
    o.ag[,sampleId:=tstrsplit(Pair_name, "\\.")[[5]]]
    o.ag <- merge(o.ag, clusters[,c("sampleId", "Continental_clusters"), with=F], by="sampleId")

    o.ag[Continental_clusters.x==Continental_clusters.y & Continental_clusters.x=="3.Europe_E", popset:="EE"]
    o.ag[Continental_clusters.x==Continental_clusters.y & Continental_clusters.x=="1.Europe_W", popset:="WW"]
    o.ag[Continental_clusters.x!=Continental_clusters.y, popset:="EW"]

    o.ag[,V1:=tstrsplit(pair, "\\.")[[1]]]
    o.ag[,V2:=tstrsplit(pair, "\\.")[[2]]]

    o.ag[,locale1:=paste(tstrsplit(V1, "_")[[1]], tstrsplit(V1, "_")[[2]], sep="_")]
    o.ag[,locale2:=paste(tstrsplit(V2, "_")[[1]], tstrsplit(V2, "_")[[2]], sep="_")]
    o.ag[,sameLocale:=locale1==locale2]

    save(o, o.ag, file="~/moments_out.Rdata")


    o.ag[]

    summary(lm(log10(divergence_time)~popset, o.ag[SFS_method=="binom"][RD_filter=="all"][SNP_caller=="PoolSNP"]))
    summary(lm(log10(divergence_time)~sameLocale*popset, o.ag[SFS_method=="binom"][RD_filter=="all"][SNP_caller=="PoolSNP"]))











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


  ### pull in samps
    o.ag[,sampleId:=tstrsplit(pair, "\\.")[1]]
    o.ag <- merge(o.ag, samps[,c("sampleId", "status"), with=F], by="sampleId")
    o.ag[,sampleId:=tstrsplit(pair, "\\.")[2]]
    o.ag <- merge(o.ag, samps[,c("sampleId", "status"), with=F], by="sampleId")




tmp <- o[,list(rMin=rank(rAIC, ties.method="min"), rMax=rank(rAIC, ties.method="max"), .N), list(Pair_name)]
tmp[,frac:=(rMax-rMin)/N]
tmp[rMin==1][N>50]
