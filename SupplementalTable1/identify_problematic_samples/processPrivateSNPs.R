### module load intel/18.0 intelmpi/18.0 R/3.6.3; R

### libraries
  library(data.table)

### load
  dat <- fread("/scratch/aob2x/SNAPE.AllSamps.0.001.delim")
  samps <- fread("/scratch/aob2x/dest/data-paper/SupplementalTable1/identify_problematic_samples/all.samps.txt", header=F)

### merge
  dat[,i:=as.numeric(gsub(";", "", V7))]
  samps[,i:=1:dim(samps)[1]]
  setnames(samps, "V1", "POP")

  priv <- merge(dat, samps, by="i")

  priv.ag <- priv[,list(.N), POP]

### save
  write.csv(priv.ag, file="/scratch/aob2x/number_privateSNPs_nofilter.csv", quote=F, row.names=F)

  scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/number_privateSNPs_nofilter.csv ~/.
