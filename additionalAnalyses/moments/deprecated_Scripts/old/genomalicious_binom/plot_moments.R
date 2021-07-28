
# scp aob2x@rivanna.hpc.virginia.edu:~/moments_o_genomalicious_binomial.Rdata ~/.

  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(patchwork)

  load("~/moments_o_genomalicious_binomial.Rdata")
  o_wide[,divergence_time:=as.numeric(as.character(divergence_time_PoolSNP))]
  o.ag[,theta:=theta/83960116]
  oo <- melt(o.ag, id.vars=c("pair", "caller", "Pair_name"))
  oow <- dcast(oo, variable+pair~caller, value.var="value")


  ps2 <- oow[variable=="pop2_size" & SNAPE > 1e6]
  ps1 <- oow[variable=="pop1_size" & SNAPE > 1e6]
  theta <- oow[variable=="theta" & PoolSNP > 6e5]
  dt <- oow[variable=="divergence_time" & PoolSNP > 5e4]

  flags <- rbindlist(list(ps2, ps1, theta, dt))

tab <- table(unlist(tstrsplit(flags$pair, "\\.")))
tab <- data.table(sampleId=names(tab), N=as.numeric(tab))
tab[,x:=rank(-N, ties="random")]

flags[,use:=apply(foreach(pi=tab[N>4]$sampleId, .combine="cbind")%do%grepl(pi, flags$pair),
                  1, any)]
setkey(oow, pair)
setkey(flags, pair)
oow[J(flags),outlier:=pair]

oow[,use:=apply(foreach(pi=tab[N>5]$sampleId, .combine="cbind")%do%ifelse(grepl(pi, oow$pair), pi, ""),
                  1, function(x) paste(x[x!=""], collapse="."))]
oow[use=="", use:="others"]

options(digits = 2, scipen = -2)
p1 <-
ggplot(data=oow, aes(x=(SNAPE), y=(PoolSNP))) +
geom_point() +
facet_wrap(~variable, scales="free") + geom_abline(slope=1, intercept=0) +
#geom_text_repel(data=flags[use==T], aes(label=pair), size=2, force=2) +
theme_bw() +
theme(legend.position="bottom")

p2 <- ggplot(data=tab, aes(x=x, y=N, fill=(N>5))) + geom_bar(stat="identity") +
  geom_text(aes(y=N+5, label=sampleId), angle=0, size=2) +
  ylim(0,25) + theme_bw() + coord_flip() +
  theme(legend.position="none") +
  xlab("Rank")

layout <- "
AAAAABB
AAAAABB
AAAAABB
"


mega <- p1 + p2 +
plot_layout(design = layout) +
plot_annotation(tag_levels = 'A')


ggsave(mega, height=5, w=7, file="~/moments.pdf")




ggplot(oo) +
geom_boxplot(aes(x=caller, y=value)) +
facet_wrap(~variable, scales="free")

geom_line(aes(x=caller, y=value, group=pair)) +
