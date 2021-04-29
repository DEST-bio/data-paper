### module load intel/18.0 intelmpi/18.0 R/3.6.3 bcftools; R

### libraries
  library(data.table)
  library(magrittr)

### load
  dat <- fread("/project/berglandlab/SNAPE.AllSamps.0.001.delim")
  #samps <- fread("/scratch/aob2x/dest/data-paper/SupplementalTable1/identify_problematic_samples/all.samps.txt", header=F)

###
  samps <- system("bcftools view -h /project/berglandlab/DEST/vcf/dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf.gz", intern=T) %>%
  last() %>%
  tstrsplit(., "\t") %>%
  unlist()
  samps <- samps[-c(1:9)]
  samps <- data.table(POP=samps)

### merge
  dat[,i:=as.numeric(gsub(";", "", V7))]
  dat.ag <- dat[,.N,i]

  samps[,i:=1:dim(samps)[1]]

  priv.ag <- merge(dat.ag, samps, by="i")


### save
  write.csv(priv.ag, file="/scratch/aob2x/number_privateSNPs_nofilter.csv", quote=F, row.names=F)

  scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/number_privateSNPs_nofilter.csv ~/.
