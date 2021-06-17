### module load gcc/7.1.0 openmpi/3.1.4 R/3.6.3
### R


### library
  library(data.table)
  library(foreach)

### load data
  fs <- list.files("/project/berglandlab/moments/moments_output", full.names=T)

### read
  o <- foreach(i=fs)%do%{
    #i<-fs[2]
    tmp <- fread(i)
  }
  o <- rbindlist(o)
  setnames(o, "-2LL_model", "LL")
  o[,list(divergence_time=divergence_time[which.max(LL)], var_dt=var(divergence_time)), list(Pair_name)]

  save(o, file="~/moments_o.Rdata")


  # scp aob2x@rivanna.hpc.virginia.edu:~/moments_o.Rdata ~/.
