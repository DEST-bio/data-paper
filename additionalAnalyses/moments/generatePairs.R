### module load gcc/7.1.0 openmpi/3.1.4 R/3.6.3
### R

### libraries
  library(data.table)
  library(SeqArray)
  library(foreach)
  library(bedr)
  library(sp)

### load samps
  setwd("/scratch/aob2x/")
  samps <- fread("DEST_freeze1/populationInfo/samps_10Nov2020.csv")

### some basic sample filtering
  samps <- samps[status=="Keep"]
  samps <- samps[propSimNorm<=0.01]

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

### subsample to ~1000 pairs evenly across the distance distribution
    setkey(1234)
    table(round(pairs$dist/20)*20)


    pairs.sample <- pairs[,list(id=rep(sample(id, 200, replace=F))),
                            list(dist.bin=round(dist/20)*20)]
    job_groups <- expand.grid(data_source=c("PoolSNP", "SNAPE"),
                                                   sfs_method=c("counts", "binom"),
                                                   rd_filter=c("all", "median_1sd"))

    pairs.sample <- pairs.sample[,list(data_source=job_groups$data_source,
                                      sfs_method=job_groups$sfs_method,
                                      rd_filter=job_groups$rd_filter,
                                        dist.bin=dist.bin),
                                  list(id)]



    pairs.sample <- merge(pairs.sample, pairs, by="id")

    table(pairs.sample$popset)

### write file
  write.csv(pairs.sample, "/scratch/aob2x/data-paper/additionalAnalyses/moments/pairs.csv", row.names=F)
