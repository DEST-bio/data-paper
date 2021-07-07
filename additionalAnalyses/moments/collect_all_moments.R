### module load gcc/7.1.0 openmpi/3.1.4 R/3.6.3
### R

### libraries
  library(data.table)
  library(SeqArray)
  library(foreach)
  library(bedr)
  library(sp)


### counts & full data-set
  fs.counts.full <- list.files("/project/berglandlab/moments/moments_output_genomalicious", full.names=T)

### binom & full data-set
  fs.binom.full <- list.files("/project/berglandlab/moments/moments_output_genomalicious_binom", full.names=T)

### counts & rdSlice
  fs.counts.full <- list.files("/project/berglandlab/moments/moments_output_genomalicious_rdSlice", full.names=T)

### binom & rdSlice
  fs.binom.full <- list.files("/project/berglandlab/moments/moments_output_genomalicious_binom_rdSlice", full.names=T)
