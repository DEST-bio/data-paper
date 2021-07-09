### scp aob2x@rivanna.hpc.virginia.edu:~/moments_out.Rdata ~/.

### libraries
  library(ggplot2)
  library(data.table)
  library(cowplot)
  library(patchwork)

### load data
  load("~/moments_out.Rdata")

### re-arrange the dimensions of hte data
  oo <- melt(o.ag[,-c("sampleId", "V1", "V2", "locale1", "locale2", "Continental_clusters.x", "Continental_clusters.y"),with=F][N>25],
              id.vars=c("pair", "SNP_caller", "SFS_method", "RD_filter", "Pair_name", "popset", "sameLocale"))
  oow <- dcast(oo, variable+pair+SFS_method+RD_filter+popset+sameLocale~SNP_caller, value.var="value")

  oow[SFS_method=="binom", SFS_method:="probs"]

### some basic stats
  o.ag[SFS_method=="binom"][RD_filter=="all"][SNP_caller=="PoolSNP"][,
      list(mu=10^mean(log10(divergence_time)),
            sd=10^sd(log10(divergence_time))),
      list(popset)]

  t.test(log10(o.ag[SFS_method=="binom"][RD_filter=="all"][SNP_caller=="PoolSNP"][popset=="EW"]$divergence_time),
         log10(o.ag[SFS_method=="binom"][RD_filter=="all"][SNP_caller=="PoolSNP"][popset=="EE"]$divergence_time))

  t.test(log10(o.ag[SFS_method=="binom"][RD_filter=="all"][SNP_caller=="PoolSNP"][popset=="EW"]$divergence_time),
        log10(o.ag[SFS_method=="binom"][RD_filter=="all"][SNP_caller=="PoolSNP"][popset=="WW"]$divergence_time))


### full set
### divergence time plot
  dt.plot <-
  ggplot(data=oow[variable%in%c("divergence_time")][RD_filter=="all"], aes(x=log10(1+SNAPE), y=log10(1+PoolSNP), color=popset)) +
  geom_point() +
  facet_grid(SFS_method~.) +
  geom_abline(slope=1, intercept=0) +
  theme_bw() +
  theme(legend.position="bottom") +
  labs(x=expression(log[10](Divergence~Time[SNAPE]+1)),
       y=expression(log[10](Divergence~Time[PoolSNP]+1)))

  theta.plot <-
  ggplot(data=oow[variable%in%c("theta")][RD_filter=="all"], aes(x=SNAPE, y=PoolSNP, color=popset)) +
  geom_point() +
  facet_grid(SFS_method~.) +
  geom_abline(slope=1, intercept=0) +
  theme_bw() +
  theme(legend.position="bottom") +
  labs(x=expression(Theta[SNAPE]),
       y=expression(Theta[PoolSNP]))


### mega-plot

  mega.plot <- dt.plot + theta.plot +
  theme(legend.position="bottom") + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A')


  ggsave(mega.plot, file="~/moments_full.pdf", height=5, w=10)

### best combination
  ### divergence time plot
    dt.plot <-
    ggplot(data=oow[variable%in%c("divergence_time")][SFS_method=="binom"][RD_filter=="all"],
            aes(x=log10(1+SNAPE), y=log10(1+PoolSNP), color=popset)) +
    geom_point() +
    facet_grid(SFS_method~RD_filter) +
    geom_abline(slope=1, intercept=0) +
    theme_bw() +
    theme(legend.position="bottom") +
    ggtitle("Divergence Time (years)") +
    ylim(0,4.2) + xlim(0,4.2)

    theta.plot <-
    ggplot(data=oow[variable%in%c("theta")][SFS_method=="binom"][RD_filter=="all"], aes(x=SNAPE, y=PoolSNP, color=popset)) +
    geom_point() +
    facet_grid(SFS_method~RD_filter) +
    geom_abline(slope=1, intercept=0) +
    theme_bw() +
    theme(legend.position="bottom") +
    ggtitle("theta") +
    ylim(0.0035,0.0065) + xlim(0.0035,0.0065)

  ### mega-plot

    mega.plot <- dt.plot + theta.plot

    ggsave(mega.plot, file="~/moments_best.pdf", height=5, w=10)


###
  o[]

### Are the clusters of theta & dt the same sets? Yes
  dt.dt <- oow[variable%in%c("divergence_time")][SFS_method=="binom"]
  theta.dt <- oow[variable%in%c("theta")][SFS_method=="binom"]

  dt.dt[,cluster:=PoolSNP<SNAPE]
  theta.dt[,cluster:=PoolSNP<SNAPE]

  m <- merge(dt.dt, theta.dt, by="pair")
  table(m$cluster.x, m$cluster.y)


### time split differences between the clusters
  oow[,popset:=factor(popset, levels=c("WW", "EW", "EE"))]

  poolSNP_divtime.boxplot <- ggplot(data=oow[variable%in%c("divergence_time")][SFS_method=="binom"][RD_filter=="all"],
        aes(x=popset, y=log10(PoolSNP+1))) +
  geom_boxplot() +
  ylab("log10(Divergence Time, years)") +
  geom_signif(comparisons = list(c("EE", "EW")), map_signif_level = F) +
  geom_signif(comparisons = list(c("WW", "EW")), map_signif_level = F, step_increase=1)

  ggsave(poolSNP_divtime.boxplot, file="~/box_plot.pdf")

  summary(lm((PoolSNP)~popset, oow[variable%in%c("divergence_time")][SFS_method=="binom"][RD_filter=="all"]))













### rank AIC
  o.rank <- o[,list(rank=rank(AIC, ties="first"),
                    deltaAIC=AIC-min(AIC, na.rm=T),
                    divergence_time, theta=theta, .N),
              list(Pair_name)]


  o.rank[,SNP_caller:=tstrsplit(Pair_name, "\\.")[[1]]]
  o.rank[,SFS_method:=tstrsplit(Pair_name, "\\.")[[2]]]
  o.rank[,RD_filter:=tstrsplit(Pair_name, "\\.")[[3]]]
  o.rank[,pair:=paste(tstrsplit(Pair_name, "\\.")[[4]], tstrsplit(Pair_name, "\\.")[[5]], sep=".")]
  o.rank[,divergence_time:=as.numeric(divergence_time)]
  o.rank[,theta:=as.numeric(theta)]
  o.rank[,deltaAIC_bin:=floor(deltaAIC/20)*20]
  #o.rank[deltaAIC_bin>1000,deltaAIC_bin:=5000]

  o.rank <- o.rank[SFS_method=="binom"][RD_filter=="all"]

  pairs_zero <- o.rank[rank==1][divergence_time<350][SNP_caller=="SNAPE"]$Pair_name
  setkey(o.rank, "Pair_name")
  o.rank[,zero:=F]
  o.rank[J(pairs_zero), zero:=T]

  rank.plot <-
  ggplot(data=o.rank[deltaAIC<1000][order(deltaAIC_bin)],
          aes(x=rank, y=log10(divergence_time), group=pair)) +
  geom_line(color="grey") +
  geom_point(aes(color=as.factor(deltaAIC_bin))) +
  facet_grid(zero~SNP_caller) +
  theme_bw()


  ggsave(rank.plot, file="~/rank_plot.png", height=5, w=10)


  o.rank.ag <- o.rank[,list(pr=mean(divergence_time<=350), deltaAIC_mean=mean(deltaAIC)),
                        list(SNP_caller, zero, rank)]

  pr.plot <- ggplot(o.rank.ag[rank<=10], aes(x=rank, y=pr)) + geom_line()+ facet_grid(zero~SNP_caller)
  dAIC.plot <- ggplot(o.rank.ag[rank<=10], aes(x=rank, y=deltaAIC_mean)) + geom_line()+ facet_grid(zero~SNP_caller)

  ggsave(pr.plot + dAIC.plot, file="~/rank_plot_ag.png", h=5, w=10)




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
