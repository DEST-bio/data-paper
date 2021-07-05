# module load gcc/7.1.0 openmpi/3.1.4 R/3.6.3; R

args = commandArgs(trailingOnly=TRUE)
job=as.numeric(args[1])
message(job)
#job<-1

### libraries
  library(data.table)
  library(SeqArray)
  library(foreach)
  library(sp)
  library(doMC)
  registerDoMC(2)
  library(genomalicious)

### load in pairs file
  pairs <- fread("/project/berglandlab/moments/pairs.csv")

  head(pairs)
  pairs[job]

### load samps
  samps <- fread("/scratch/aob2x/DEST_freeze1/populationInfo/samps_10Nov2020.csv")

### open GDS file & make SNP table
  if (pairs[job]$type=="PoolSNP") {
    #q(save="no")
    message("PoolSNP")
    genofile <- seqOpen(paste("/project/berglandlab/DEST/gds/dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.gds", sep=""))

  } else if (pairs[job]$type=="SNAPE") {
    message("SNAPE")
    #genofile <- openfn.gds(paste("/project/berglandlab/DEST/gds/dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.gds", sep=""))
    genofile <- seqOpen(paste("/project/berglandlab/DEST/gds/dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.gds", sep=""))

    #save(snps.dt, file="/project/berglandlab/moments/SNAPE.snp.dt.Rdata")

  }

  message("making snp table")
  snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                        pos=seqGetData(genofile, "position"),
                        variant.id=seqGetData(genofile, "variant.id"),
                        nAlleles=seqNumAllele(genofile),
                        missing=seqMissing(genofile, .progress=T))

  ## choose number of alleles
   snps.dt <- snps.dt[nAlleles==2]

### get polymorphism data
  message("get poly")
  setkey(snps.dt, chr)
  seqSetFilter(genofile, sample.id=as.character(pairs[job, c("V1", "V2"), with=F]),
                variant.id=snps.dt[J(c("2L", "2R", "3L", "3R"))]$variant.id)


### does file already exist?
  fn <- paste("/scratch/aob2x/moments_binom/",
            seqGetData(genofile, "sample.id")[1],
            ".",
            seqGetData(genofile, "sample.id")[2],
            ".",
            pairs[job]$type,
            ".delim", sep="")

  if(file.exists(fn)) q(save="no")

### if not, do the rest

  ### get allele frequency data
  ad <- seqGetData(genofile, "annotation/format/AD")
  dp <- seqGetData(genofile, "annotation/format/DP")

  dat <- ad$data/dp$data
  dim(dat)
  rownames(dat) <- seqGetData(genofile, "sample.id")
  dat <- t(dat)
  #dat <- na.omit(dat)

### downsample to Neff
  ### load average effective read depth
    dep <- fread("/scratch/aob2x/DEST_freeze1/populationInfo/sequencingStats/rd.csv")[auto==T]
    setkey(dep, sampleId)
    setkey(samps, sampleId)
    neff <- merge(dep[J(colnames(dat))], samps[J(colnames(dat))], by="sampleId")[,c("sampleId", "nFlies", "mu.25"), with=F]
    neff[,ne:=round((2*nFlies*mu.25)/(2*nFlies+mu.25))]
    neff[,ne:=floor(ne/2)] ### note, this is the number of diploid individuals...


### convert to format for genomalicious
  dat <- as.data.table(dat)
  dat[,locus:=seqGetData(genofile, "variant.id")]
  dat[,ref:=seqGetData(genofile, "$ref")]
  dat[,alt:=seqGetData(genofile, "$alt")]

  datl <- melt(dat, id.vars=c("locus", "ref", "alt"))
  setnames(datl, c("variable", "value"), c("POOL", "FREQ"))

  datl <- merge(datl, neff[,c("sampleId", "ne")], by.x="POOL", by.y="sampleId")
  datl[,FREQ:=1-FREQ]
  datl <- na.omit(datl)

  dadi <- dadi_inputs_pools(
    datl,
    poolCol = "POOL",
    locusCol = "locus",
    refCol = "ref",
    altCol = "alt",
    freqCol = "FREQ",
    indsCol = "ne",
    poolSub = NULL,
    methodSFS = "probs"
  )
  dadi <- na.omit(dadi)

  fn <- paste("/scratch/aob2x/moments_binom/",
            seqGetData(genofile, "sample.id")[1],
            ".",
            seqGetData(genofile, "sample.id")[2],
            ".",
            pairs[job]$type,
            ".delim", sep="")

  write.table(dadi,
              file=fn,
              sep="\t", quote=F, row.names=F)

### write meta-data file
  meta <- data.table(id=paste(seqGetData(genofile, "sample.id"), collapse="|"),
                     file=fn,
                     L=83960116,
                     pop1=seqGetData(genofile, "sample.id")[1],
                     pop2=seqGetData(genofile, "sample.id")[2],
                     projection1=neff[sampleId==seqGetData(genofile, "sample.id")[1]]$ne*2,
                     projection2=neff[sampleId==seqGetData(genofile, "sample.id")[2]]$ne*2)

    meta.fn <- paste("/scratch/aob2x/moments_binom/",
              seqGetData(genofile, "sample.id")[1],
              ".",
              seqGetData(genofile, "sample.id")[2],
              ".",
              pairs[job]$type,
              ".meta", sep="")

    write.table(meta,
                file=meta.fn,
                sep="\t", quote=F, row.names=F, col.names=F)
