### scp aob2x@rivanna.hpc.virginia.edu:~/moments_out.Rdata ~/.

### libraries
  library(ggplot2)
  library(data.table)
  library(cowplot)
  library(patchwork)

### load data
  load("~/moments_out.Rdata")

### re-arrange the dimensions of hte data
  oo <- melt(o.ag[,-c("sampleId", "status.x", "status.y"),with=F][N>25], id.vars=c("pair", "SNP_caller", "SFS_method", "RD_filter", "Pair_name"))
  oow <- dcast(oo, variable+pair+SFS_method+RD_filter~SNP_caller, value.var="value")



### full set
### divergence time plot
  dt.plot <-
  ggplot(data=oow[variable%in%c("divergence_time")], aes(x=SNAPE, y=PoolSNP)) + geom_point() +
  facet_grid(SFS_method~RD_filter) +
  geom_abline(slope=1, intercept=0) +
  theme_bw() +
  theme(legend.position="bottom") +
  ggtitle("Divergence Time (years)")

  theta.plot <-
  ggplot(data=oow[variable%in%c("theta")], aes(x=SNAPE, y=PoolSNP)) + geom_point() +
  facet_grid(SFS_method~RD_filter) +
  geom_abline(slope=1, intercept=0) +
  theme_bw() +
  theme(legend.position="bottom") +
  ggtitle("theta")

### mega-plot

  mega.plot <- dt.plot + theta.plot

  ggsave(mega.plot, file="~/moments_full.pdf", height=5, w=10)


### best combination
  ### divergence time plot
    dt.plot <-
    ggplot(data=oow[variable%in%c("divergence_time")][SFS_method=="binom"][RD_filter=="all"], aes(x=SNAPE, y=PoolSNP)) + geom_point() +
    facet_grid(SFS_method~RD_filter) +
    geom_abline(slope=1, intercept=0) +
    theme_bw() +
    theme(legend.position="bottom") +
    ggtitle("Divergence Time (years)")

    theta.plot <-
    ggplot(data=oow[variable%in%c("theta")][SFS_method=="binom"][RD_filter=="all"], aes(x=SNAPE, y=PoolSNP)) + geom_point() +
    facet_grid(SFS_method~RD_filter) +
    geom_abline(slope=1, intercept=0) +
    theme_bw() +
    theme(legend.position="bottom") +
    ggtitle("theta")

  ### mega-plot

    mega.plot <- dt.plot + theta.plot

    ggsave(mega.plot, file="~/moments_best.pdf", height=5, w=10)


### Are the clusters of theta & dt the same sets?
  dt.dt <- oow[variable%in%c("divergence_time")][SFS_method=="binom"]
  theta.dt <- oow[variable%in%c("theta")][SFS_method=="binom"]

  dt.dt[,cluster:=PoolSNP<SNAPE]
  theta.dt[,cluster:=PoolSNP<SNAPE]

  m <- merge(dt.dt, theta.dt, by="pair")
  table(m$cluster.x, m$cluster.y)



### rank AIC
  o.rank <- o[,list(rank=rank(AIC, ties="first"),
                    deltaAIC=AIC-min(AIC, na.rm=T),
                    divergence_time, theta=theta/L, .N),
              list(Pair_name)]


  o.rank[,SNP_caller:=tstrsplit(Pair_name, "\\.")[[1]]]
  o.rank[,SFS_method:=tstrsplit(Pair_name, "\\.")[[2]]]
  o.rank[,RD_filter:=tstrsplit(Pair_name, "\\.")[[3]]]
  o.rank[,pair:=paste(tstrsplit(Pair_name, "\\.")[[4]], tstrsplit(Pair_name, "\\.")[[5]], sep=".")]

  om <- merge(o.rank, theta.dt[SFS_method=="binom"][RD_filter=="all"], by="pair")


  rank_plot <-
  ggplot(data=om[grepl("binom.all", Pair_name)][!is.na(cluster)],
        aes(x=rank, y=divergence_time, group=Pair_name)) +
  geom_line() +
  facet_grid(SNP_caller~cluster) +
  ylim(0, 20000)
  ggsave(rank_plot, file="~/rank_plot.pdf")


  delta_plot <-
  ggplot(data=om[grepl("binom.all", Pair_name)][!is.na(cluster)],
        aes(x=deltaAIC, y=divergence_time, group=Pair_name)) +
  geom_line() +
  geom_point() +
  facet_grid(SNP_caller~cluster) +
  ylim(0, 10000) +
  xlim(0, 500)
  ggsave(delta_plot, file="~/delta_plot.pdf")


  delta_rank_plot <-
  ggplot(data=om[grepl("binom.all", Pair_name)][!is.na(cluster)],
        aes(x=log10(deltaAIC+1), y=rank, group=Pair_name)) +
  geom_line() +
  facet_grid(SNP_caller~cluster)

  ggsave(delta_rank_plot, file="~/delta_rank_plot.pdf")
[

  pp <-
  ggplot(data=om[grepl("binom.all", Pair_name)][!is.na(cluster)][order(-deltaAIC)][rank<=5],
        aes(x=divergence_time, y=theta, color=log10(deltaAIC+1))) +
  geom_point() +
  facet_wrap(SNP_caller~cluster) +
  ylim(0.003, 0.01) +
  xlim(0, 10000)

  ggsave(pp, file="~/pp.pdf")


  bp <-
