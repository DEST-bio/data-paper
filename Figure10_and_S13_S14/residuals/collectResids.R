library(data.table)
library(ggplot2)
library(lattice)
dat <- as.matrix(fread("/Users/alanbergland/PoolSNP.binom.AT_See_16_1.AT_Mau_15_50.EW_widebounds_split_mig_asymm_resids.delim", skip=1))
datl <- data.table(resid=expand.grid(dat[,-1])[,1],
                   n1=rep(1:dim(dat[,-1])[1], dim(dat[,-1])[2]),
                   n2=rep(1:dim(dat[,-1])[2], each=dim(dat[,-1])[1]))
datl <- na.omit(datl)

ggplot(datl, aes(resid)) + geom_histogram()

ggplot(datl, aes(x=n1, y=n2, fill=resid)) + geom_tile()
histogram(datl[resid>-60 & resid<60]$resid)
levelplot(dat[,-1])

obs<-as.matrix(fread("/Users/alanbergland/PoolSNP.binom.AT_See_16_1.AT_Mau_15_50.EW_widebounds_split_mig_asymm_resids.data", skip=1))
levelplot(log10(obs[,-1]))


model<-as.matrix(fread("/Users/alanbergland/PoolSNP.binom.AT_See_16_1.AT_Mau_15_50.EW_widebounds_split_mig_asymm_resids.model", skip=1))
levelplot(log10(model[,-1]*417843.212514047))


417843.212514047
