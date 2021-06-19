library(ggplot2)
library(data.table)
library(foreach)

### load samps
  samps <- fread("/Users/alanbergland/Documents/GitHub/DEST_freeze1/populationInfo/samps_10Nov2020.csv")
  dep <- fread("/Users/alanbergland/Documents/GitHub/DEST_freeze1/populationInfo/sequencingStats/rd.csv")[auto==T]

### calculate denimonator
  neff <- merge(dep, samps, by="sampleId")[,c("sampleId", "nFlies", "mu.25", "set"), with=F]
  neff[,ne:=round((2*nFlies*mu.25)/(2*nFlies+mu.25))]
  neff[,ne:=floor(ne/2)*2]
  neff <- neff[set!="dgn"]

###
  afs <- seq(from=.001, to=.5, by=.0002)

  o <- foreach(i=1:dim(neff)[1], .combine="rbind")%do%{
    data.table(af=afs, sampleId=neff[i]$sampleId, fs=round(neff[i]$ne*afs), set=neff[i]$set)
  }

  ggplot(data=o, aes(x=af, y=fs, group=sampleId)) + geom_line() + facet_wrap(~set)


  o.ag <- o[,list(prop=mean(fs==0)), list(af, set)]

  ggplot(data=o.ag, aes(x=af, y=prop)) + geom_point() + facet_wrap(~set)
