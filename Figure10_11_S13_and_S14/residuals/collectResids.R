### module load intel/18.0 intelmpi/18.0 R/3.6.3; R

### libraries
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(10)

### load metadata file
  resid.meta <- fread("/project/berglandlab/moments/resid_meta.delim")

### load resid files
  fl <- list.files("/project/berglandlab/moments/resids/", "delim", full=T)

  allResid <- foreach(fl.i=fl)%dopar%{
    message(paste(which(fl.i==fl), length(fl), sep=" / "))
    #fl.i<-fl[1]
    resid.tmp <- fread(fl.i)
    resid.tmp <- resid.tmp[,-1][-1,]

    resid.tmp.l <- data.table(resid=expand.grid(as.matrix(resid.tmp))[,1],
                       n1=rep(1:dim(resid.tmp)[1], dim(resid.tmp)[2]),
                       n2=rep(1:dim(resid.tmp)[2], each=dim(resid.tmp)[1]))

    fl.is <- last(tstrsplit(fl.i, "/"))
    fl.is <- tstrsplit(fl.is, "\\.")

    resid.tmp.l[,Pair_name:=paste(unlist(fl.is[1:5]), collapse=".")]
    resid.tmp.l[,param:=fl.is[[6]]]
    resid.tmp.l[,model_type:=fl.is[[7]]]
    resid.tmp.l[,model_sym:=fl.is[[8]]]
    resid.tmp.l[,max_n1:=max(n1)]
    resid.tmp.l[,max_n2:=max(n2)]
    resid.tmp.l[,pop_name1:=fl.is[[3]]]
    resid.tmp.l[,pop_name2:=fl.is[[4]]]
    resid.tmp.l
  }
  allResid <- rbindlist(allResid)

### save
  save(allResid, file="~/allResid.Rdata")
