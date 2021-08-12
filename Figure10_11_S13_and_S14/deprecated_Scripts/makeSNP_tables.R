#module load gcc/7.1.0 openmpi/3.1.4 R/3.6.3; R

### libraries
  library(data.table)
  library(SeqArray)


### PoolSNP
  message("PoolSNP")
  genofile <- seqOpen(paste("/project/berglandlab/DEST/gds/dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.gds", sep=""))
  message("making snp table")
  snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                        pos=seqGetData(genofile, "position"),
                        variant.id=seqGetData(genofile, "variant.id"),
                        nAlleles=seqNumAllele(genofile),
                        missing=seqMissing(genofile, .progress=T))

  ## choose number of alleles
   snps.dt <- snps.dt[nAlleles==2]

  ### save
    save(snps.dt, file="/scratch/aob2x/moments_general/snp_table_PoolSNP.Rdata")

  seqClose(genofile)

### SNAPE
  genofile <- seqOpen(paste("/project/berglandlab/DEST/gds/dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.gds", sep=""))

  message("making snp table")
  snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                        pos=seqGetData(genofile, "position"),
                        variant.id=seqGetData(genofile, "variant.id"),
                        nAlleles=seqNumAllele(genofile),
                        missing=seqMissing(genofile, .progress=T))

  ## choose number of alleles
   snps.dt <- snps.dt[nAlleles==2]


 ### save
   save(snps.dt, file="/scratch/aob2x/moments_general/snp_table_SNAPE.Rdata")

 seqClose(genofile)
