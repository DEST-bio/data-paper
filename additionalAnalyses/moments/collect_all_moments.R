### module load gcc/7.1.0 openmpi/3.1.4 R/3.6.3
### R

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
  o.ag <- o[,list(divergence_time=divergence_time[which.min((AIC))],
                  theta=theta[which.min((AIC))]/L[which.min((AIC))],
                  pop1_size=pop1_size[which.min((AIC))],
                  pop2_size=pop2_size[which.min((AIC))],
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


### pull in samps
  o.ag[,sampleId:=tstrsplit(pair, "\\.")[1]]
  o.ag <- merge(o.ag, samps[,c("sampleId", "status"), with=F], by="sampleId")
  o.ag[,sampleId:=tstrsplit(pair, "\\.")[2]]
  o.ag <- merge(o.ag, samps[,c("sampleId", "status"), with=F], by="sampleId")

  save(o, o.ag, file="~/moments_out.Rdata")






tmp <- o[,list(rMin=rank(rAIC, ties.method="min"), rMax=rank(rAIC, ties.method="max"), .N), list(Pair_name)]
tmp[,frac:=(rMax-rMin)/N]
tmp[rMin==1][N>50]
