# module load gcc/7.1.0 openmpi/3.1.4 R/3.6.3; R

args = commandArgs(trailingOnly=TRUE)
job=as.numeric(args[1])
message(job)
#job<-2

### libraries
  library(data.table)
  library(SeqArray)
  library(foreach)
  library(sp)
  library(doMC)
  registerDoMC(2)

### load in pairs file
  #pairs <- fread("/project/berglandlab/moments/pairs.csv")
  pairs <- fread("/scratch/aob2x/pairs.csv")
  head(pairs)
  pairs[job]

### load samps
  samps <- fread("/scratch/aob2x/DEST_freeze1/populationInfo/samps_10Nov2020.csv")

### open GDS file & make SNP table
  if (pairs[job]$type=="PoolSNP") {
    #q(save="no")
    message("PoolSNP")
    genofile <- seqOpen(paste("/project/berglandlab/DEST/gds/dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.gds", sep=""))

    ### run once
      #message("making snp table")
      #snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
      #                      pos=seqGetData(genofile, "position"),
      #                      variant.id=seqGetData(genofile, "variant.id"),
      #                      nAlleles=seqNumAllele(genofile),
      #                      missing=seqMissing(genofile, .progress=T))
#
      ### choose number of alleles
      # snps.dt <- snps.dt[nAlleles==2]
      # save(snps.dt, file="/project/berglandlab/moments/PoolSNP.snp.dt.Rdata")
      load(file="/project/berglandlab/moments/PoolSNP.snp.dt.Rdata")


  } else if (pairs[job]$type=="SNAPE") {
    message("SNAPE")
    genofile <- openfn.gds(paste("/project/berglandlab/DEST/gds/dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.gds", sep=""))
    #genofile <- seqOpen(paste("/project/berglandlab/DEST/gds/dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.gds", sep=""))

    message("making snp table")
    snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                          pos=seqGetData(genofile, "position"),
                          variant.id=seqGetData(genofile, "variant.id"),
                          nAlleles=seqNumAllele(genofile),
                          missing=seqMissing(genofile, .progress=T))

    ## choose number of alleles
     snps.dt <- snps.dt[nAlleles==2]
     #save(snps.dt, file="/project/berglandlab/moments/SNAPE.snp.dt.Rdata")

  }


### get polymorphism data
  message("get poly")
  setkey(snps.dt, chr)
  seqSetFilter(genofile, sample.id=as.character(pairs[job, c("V1", "V2"), with=F]),
                variant.id=snps.dt[J(c("2L", "2R", "3L", "3R"))]$variant.id)

  ### get allele frequency data
  ad <- seqGetData(genofile, "annotation/format/AD")
  dp <- seqGetData(genofile, "annotation/format/DP")

  dat <- ad$data/dp$data
  dim(dat)
  rownames(dat) <- seqGetData(genofile, "sample.id")
  dat <- t(dat)
  dat <- na.omit(dat)

### fold
  #f.hat <- dat[,1]/2 + dat[,2]/2
  #
  #dat[f.hat>.5, 1] <- 1 - dat[f.hat>.5, 1]
  #dat[f.hat>.5, 2] <- 1 - dat[f.hat>.5, 2]

### turn into integers
  message("integerize")
  ### load average effective read depth
    dep <- fread("/scratch/aob2x/DEST_freeze1/populationInfo/sequencingStats/rd.csv")[auto==T]
    setkey(dep, sampleId)
    setkey(samps, sampleId)
    neff <- merge(dep[J(colnames(dat))], samps[J(colnames(dat))], by="sampleId")[,c("sampleId", "nFlies", "mu.25"), with=F]
    neff[,ne:=round((2*nFlies*mu.25)/(2*nFlies+mu.25))]
    neff[,ne:=floor(ne/2)*2]



    neff <- merge(dep, samps, by="sampleId")[,c("sampleId", "nFlies", "mu.25"), with=F]


  ### integerize
    dat[,1] <- round(dat[,1]*neff$ne[1])
    dat[,2] <- round(dat[,2]*neff$ne[2])

### turn in to 2D-SFS
  sfs <- foreach(i=0:(neff$ne[1]), .combine="rbind")%do%{
    foreach(j=0:(neff$ne[2]), .combine="rbind")%dopar%{
      #i<-
      print(paste(i, j, sep=" / "))
      data.table(N=sum(dat[,1]==i & dat[,2]==j), i=i, j=j)
    }
  }

#ggplot(data=sfs, aes(x=i, y=j, fill=log10(N))) + geom_tile()

### make filter
  sfs.filter <- as.numeric(sfs$N==0)
  sfs.filter[1] <- 1

### write output
  colnames(dat)
  setwd(paste("/project/berglandlab/moments/", pairs[job]$type, sep=""))
  getwd()
  fileConn <- file(paste(paste(c(colnames(dat), "unfolded"), collapse="."), ".", pairs[job]$type, ".fs", sep=""))
  message(paste(paste(c(colnames(dat), "unfolded"), collapse="."), ".", pairs[job]$type, ".fs", sep=""))
  writeLines(c(paste(c(neff$ne[1]+1, neff$ne[2]+1, "unfolded", dQuote(  colnames(dat), F)), collapse=" "),
               paste(sfs$N, collapse=" "),
               paste(sfs.filter, collapse=" ")), fileConn)
  close(fileConn)
