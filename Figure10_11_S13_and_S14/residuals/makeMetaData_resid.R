### module load intel/18.0 intelmpi/18.0 R/3.6.3; R

### libraries
  library(data.table)
  library(foreach)

### load data
  load("/project/berglandlab/moments_jcbn_keric/Model_test/o.best.model_search.Rdata")
  bm <- as.data.table(finding_the_best_model)
  bm[inference_method=="bound", inference_method:="widebounds"]
  bm[Demo_model=="SM", model_type:="split_mig"]
  bm[Demo_model=="IM", model_type:="IMbg"]
  bm[,index:=c(1:dim(bm)[1])]

### get pool sizes from metadata files
  neff <- foreach(i=bm$index, .combine="rbind")%do%{
    message(i)
    fn <- gsub("delim", "meta", bm[i]$fs_file)
    meta <- fread(fn)
    data.table(index=i, pop_n1=meta$V6, pop_n2=meta$V7)
  }
  bm <- merge(bm, neff, by="index")

### format table
  resid.meta <- data.table(
    fs_file= bm$fs_file,
    Pair_name=bm$Pair_name,
    param=bm$inference_method,
    model_type=bm$model_type,
    model_sym=bm$migration_model,
    pop_name1=bm$pop1,
    pop_name2=bm$pop2,
    pool_n1=bm$pop_n1,
    pool_n2=bm$pop_n2,
    nu1=bm$nu1,
    nu2=bm$nu2,
    m12=bm$m12,
    m21=bm$m21,
    Ts=bm$Ts,
    theta=bm$theta,
    nu1B=bm$nu1B,
    nu2B=bm$nu2B,
    nu1F=bm$nu1F,
    nu2F=bm$nu2F)

### export
  write.table(resid.meta, sep="\t", row.names=F, quote=F, file="/project/berglandlab/moments/resid_meta.delim")
