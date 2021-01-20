library(data.table)


### set working directory
#setwd("/scratch/aob2x/dest")
setwd("/Users/alanbergland/Documents/GitHub/")

### load meta data
samps <- fread("./DEST_freeze1/populationInfo/samps.csv")

### load sequencing summary stats
pcr <- fread("./DEST_freeze1/populationInfo/sequencingStats/pcr.csv")
rd <- fread("./DEST_freeze1/populationInfo/sequencingStats/rd.csv")[auto==T]
simulans <- fread("./DEST_freeze1/populationInfo/sequencingStats/simulans.csv")[auto==T]

setkey(rd, sampleId)
setkey(pcr, sampleId)
setkey(simulans, sampleId)
mps <- merge(rd, pcr)
mps <- merge(mps, simulans)

### quick fix for Ukrainian samples
  mps[sampleId=="UA_Pir_14_26", sampleId:="UA_Pyr_14_26"]
  mps[sampleId=="UA_Pir_15_21", sampleId:="UA_Pyr_15_21"]
  mps[sampleId=="UA_Pyr_16_48", sampleId:="UA_Pir_16_48"]

mps <- merge(mps, samps[,c("sampleId", "nFlies")], all.y=T, by="sampleId")
setnames(mps, "mu.25", "AveReadDepth.25")

### a few small fixes

  mps[,effRD.25:=(AveReadDepth.25 * 2*nFlies) / (AveReadDepth.25 + 2*nFlies)]

### rank x-axis to mean read depth
  mps <- mps[auto.x==T]
  mps[,x:=rank(AveReadDepth.25, ties.method="first")]
  mps[,x.id:=factor(sampleId, levels=mps$sampleId[mps$x])]

  mps[,propMissing:=nmissing.25/10000]

  setnames(mps, "AveReadDepth.25", "AveReadDepth")
  setnames(mps, "auto.x", "auto")
  setnames(mps, "effRD.25", "effRD")
  setnames(mps, "nmissing.25", "nmissing")

  write.csv(mps[,-c("x", "x.id", "auto.y", "q", "nmissing", "nFlies", "auto")],
            file="data-paper/SupplementalTable2/SupplementalTable2.csv",
            quote=F, row.names=F)
