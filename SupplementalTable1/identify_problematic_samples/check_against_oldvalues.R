library(googlesheets4)
library(data.table)

load("~/priv_ag.Rdata")

old_priv <- as.data.frame(read_sheet(ss="https://docs.google.com/spreadsheets/d/1PJEQE6Pn8H_vwA0cgn8e4nECrWrSGzKj6kgg7b8a3Ms/edit#gid=0",
                            sheet="Hoja 2"))
old_priv$old_N<-unlist(old_priv[,2])
old_priv$old_N<-as.numeric(as.character(old_priv$old_N))
old_priv <- as.data.table(old_priv)
setnames(old_priv, "Private SNPs in VCF", "old_priv")
m <- merge(priv.ag, old_priv, by.x="POP", by.y="Sample")

plot(N~old_N, m)
plot(I(sort(m$N))~sort(m$N))
