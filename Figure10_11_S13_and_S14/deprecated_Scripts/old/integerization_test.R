# module load gcc/7.1.0 openmpi/3.1.4 R/3.6.3; R

### libraries
  library(data.table)
  library(SeqArray)
  library(foreach)
  library(sp)
  library(doMC)
  registerDoMC(40)

### load samps
  samps <- fread("/scratch/aob2x/DEST_freeze1/populationInfo/samps_10Nov2020.csv")
  samps <- samps[set!="dgn"]

### get SNP tables
  ### PoolSNP
    genofile.poolsnp <- seqOpen(paste("/project/berglandlab/DEST/gds/dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.gds", sep=""), allow.duplicate=TRUE)
    snps.dt.poolsnp <- data.table(chr=seqGetData(genofile.poolsnp, "chromosome"),
                          pos=seqGetData(genofile.poolsnp, "position"),
                          variant.id=seqGetData(genofile.poolsnp, "variant.id"),
                          nAlleles=seqNumAllele(genofile.poolsnp),
                          missing=seqMissing(genofile.poolsnp, .progress=T))
    snps.dt.poolsnp <- snps.dt.poolsnp[nAlleles==2]

  ### SNAPE
    genofile.snape <- seqOpen(paste("/project/berglandlab/DEST/gds/dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.gds", sep=""), allow.duplicate=FALSE)
    snps.dt.snape <- data.table(chr=seqGetData(genofile.snape, "chromosome"),
                          pos=seqGetData(genofile.snape, "position"),
                          variant.id=seqGetData(genofile.snape, "variant.id"),
                          nAlleles=seqNumAllele(genofile.snape),
                          missing=seqMissing(genofile.snape, .progress=T))
    snps.dt.snape <- snps.dt.snape[nAlleles==2]


### function
  getConvertedSites <- function(pop="AT_Mau_14_01", type="PoolSNP") {
    if(type=="PoolSNP") {
      ### get polymorphism data
        message("get poly")
        setkey(snps.dt.poolsnp, chr)
        seqSetFilter(genofile.poolsnp, sample.id=pop,
                      variant.id=snps.dt.poolsnp[J(c("2L", "2R", "3L", "3R"))]$variant.id)

        ### get allele frequency data
        ad <- seqGetData(genofile.poolsnp, "annotation/format/AD")
        dp <- seqGetData(genofile.poolsnp, "annotation/format/DP")
        dat <- ad$data/dp$data

    } else if(type=="SNAPE") {
      message("get poly")
      setkey(snps.dt.snape, chr)
      seqSetFilter(genofile.snape, sample.id=pop,
                    variant.id=snps.dt.snape[J(c("2L", "2R", "3L", "3R"))]$variant.id)

      ### get allele frequency data
      ad <- seqGetData(genofile.snape, "annotation/format/AD")
      dp <- seqGetData(genofile.snape, "annotation/format/DP")
      dat <- ad$data/dp$data

    }

      dim(dat)
      rownames(dat) <- pop
      dat <- t(dat)
      dat <- na.omit(dat)

    message("integerize")
    ### load average effective read depth
      dep <- fread("/scratch/aob2x/DEST_freeze1/populationInfo/sequencingStats/rd.csv")[auto==T]
      setkey(dep, sampleId)
      setkey(samps, sampleId)

      neff <- merge(dep[J(colnames(dat))], samps[J(colnames(dat))], by="sampleId")[,c("sampleId", "nFlies", "mu.25"), with=F]
      neff[,ne:=round((2*nFlies*mu.25)/(2*nFlies+mu.25))]
      neff[,ne:=floor(ne/2)*2]


    ### integerize
      dat.i <- data.table(af=dat[,1], int=round(dat[,1]*neff$ne[1]))
      dat.i <- dat.i[af>0]

    ### return
      data.table(pop=pop, type=type, fracLost=mean(dat.i$int==0), nPoly=dim(dat.i)[1])

  }

### iterate
  o <- foreach(pop.i=samps[1:40]$sampleId, .combine="rbind")%dopar%{
    foreach(type.i=c("PoolSNP", "SNAPE"), .combine="rbind")%do%{
      getConvertedSites(type=type.i, pop=pop.i)
    }
  }

save(o, file="~/nLost.Rdata")

### scp aob2x@rivanna.hpc.virginia.edu:~/nLost.Rdata ~/.

library(data.table)
library(ggplot2)
load("~/nLost.Rdata")

ggplot(data=o, aes(x=nPoly, y=fracLost, color=type)) + geom_point()























library(ggplot2)
library(data.table)
library(foreach)

### load samps
  samps <- fread("/Users/alanbergland/Documents/GitHub/DEST_freeze1/populationInfo/samps_10Nov2020.csv")
  dep <- fread("/Users/alanbergland/Documents/GitHub/DEST_freeze1/populationInfo/sequencingStats/rd.csv")[auto==T]

### calculate denimonator
  neff <- merge(dep, samps, by="sampleId")[,c("sampleId", "nFlies", "mu.25", "set"), with=F]
  neff[,ne:=round((2*nFlies*mu.25)/(2*nFlies+mu.25))]
  neff[,ne:=floor(ne/2)*2]
  neff <- neff[set!="dgn"]

###
  afs <- seq(from=.001, to=.5, by=.0002)

  o <- foreach(i=1:dim(neff)[1], .combine="rbind")%do%{
    data.table(af=afs, sampleId=neff[i]$sampleId, fs=round(neff[i]$ne*afs), set=neff[i]$set)
  }

  ggplot(data=o, aes(x=af, y=fs, group=sampleId)) + geom_line() + facet_wrap(~set)


  o.ag <- o[,list(prop=mean(fs==0)), list(af, set)]

  ggplot(data=o.ag, aes(x=af, y=prop)) + geom_point() + facet_wrap(~set)
