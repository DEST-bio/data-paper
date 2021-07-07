### scp aob2x@rivanna.hpc.virginia.edu:~/moments_out.Rdata ~/.

### libraries
  library(ggplot2)
  library(data.table)

### load data
  load("~/moments_out.Rdata")

### rescale theta
  o.ag[,theta:=theta/83960116]

### re-arrange the dimensions of hte data
  oo <- melt(o.ag[,-"sampleId",with=F][N>25], id.vars=c("pair", "SNP_caller", "SFS_method", "RD_filter", "Pair_name", "status.x", "status.y"))
  oow <- dcast(oo, variable+pair+SFS_method+RD_filter~SNP_caller, value.var="value")


### divergence time plot
  ggplot(data=oow, aes(x=SNAPE, y=PoolSNP)) + geom_point() +
  facet_wrap(variable~SFS_method+RD_filter, scale="free", ncol=4) +
  geom_abline(slope=1, intercept=0) +
  theme_bw() +
  theme(legend.position="bottom")
  
