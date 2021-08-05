### module load gcc/7.1.0 bedops/2.4.1 bedtools/2.26.0 tabix  openmpi/3.1.4 R/3.6.3
### R

### libraries
  library(data.table)

### load samps
  setwd("/scratch/aob2x/")
  samps <- fread("DEST_freeze1/populationInfo/samps_10Nov2020.csv")


### get data
  getL <- function(pop1) {
    #pop1 <- "AT_Mau_14_01"
    bed1 <- fread(paste("/project/berglandlab/DEST/dest_mapped/pipeline_output/", pop1, "/", pop1, ".bed.gz", sep=""))
    setkey(bed1, V1)
    bed1 <- bed1[J(c("2L", "2R", "3L", "3R"))]

    sum(bed1$V3-bed1$V2)

  }
  samps[set=="DrosEU",nMasked:=foreach(p=samps[set=="DrosEU"]$sampleId, .combine="c")%dopar%getL(pop1=p)]

  median(samps$nMasked, na.rm=T) ## == 25030090

### get autosomal genome-size
  #head -n5 /scratch/aob2x/dest/referenceGenome/r6/holo_dmel_6.12.dict | sed '1d'

  #@SQ     SN:2L   LN:23513712     M5:b6a98b7c676bdaa11ec9521ed15aff2b     UR:file:/mnt/spicy_2/dest/pipeline/reference/holo_dmel_6.12.fa
  #@SQ     SN:2R   LN:25286936     M5:2ecce4ca1c0106bb7c63a78b93ab49ba     UR:file:/mnt/spicy_2/dest/pipeline/reference/holo_dmel_6.12.fa
  #@SQ     SN:3L   LN:28110227     M5:3c3ea06b22af8cc59809dbf8d154791e     UR:file:/mnt/spicy_2/dest/pipeline/reference/holo_dmel_6.12.fa
  #@SQ     SN:3R   LN:32079331     M5:420540d26d86fbf14185d2f2d68af9c4     UR:file:/mnt/spicy_2/dest/pipeline/reference/holo_dmel_6.12.fa

  (23513712+25286936+28110227+32079331) - 25030090  ### == 83960116
