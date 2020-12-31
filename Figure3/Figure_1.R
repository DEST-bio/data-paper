library(ggplot2)
library(cowplot)
library(data.table)

setwd("/Users/alanbergland/Documents/GitHub/")


### map
  ### made run this script to generate input data: /DEST/populationInfo/makePlots.R

  ### save figure objects :
  # save(world, world.plot, samps.ag.ag, file="DEST/figures/figure_1/worldMap.Rdata")

  load(file="DEST/figures/figure_1/worldMap.Rdata")
  samps.ag.ag[set=="dgn", set:="DGN"]
  	world.plot <- 	ggplot() +
  				geom_polygon(data = world,
  							aes(x=long, y = lat, group = group), fill="lightgrey") +
  				geom_point(data = samps.ag.ag,
  							aes(x=long, y=lat, size=I((n-1)/2 + 4), color=set), alpha=.5) +
  				xlab("Longitude") + ylab("Latitude") + scale_fill_manual(values="black") +
          theme_cowplot()


### timeline
  ### made here: DEST/figures/figure_1/worldMap.Rdata
  #		save(samps, samps.ag, multi_sample, file="DEST/figures/figure_1/timeline.Rdata")
  load(file="DEST/figures/figure_1/timeline.Rdata")

      samps[,season:=factor(season, levels=c("spring", "fall", "frost"))]

  		multi_sample <- ggplot() +
  		geom_line(data= samps[J(samps.ag[maxDelta>10])], aes(x=as.Date(yday, origin = as.Date("2018-01-01")), y=lat, group=locality, linetype=continent)) +
  		geom_point(data=samps[J(samps.ag[maxDelta>10])], aes(x=as.Date(yday, origin = as.Date("2018-01-01")), y=lat, group=locality, color=season)) +
  		facet_grid(.~year) +
      theme_cowplot() +
  		theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="right", legend.direction="vertical") +
  		scale_x_date(date_labels = "%b", limits = as.Date(c(110,355), origin = as.Date("2018-01-01"))) +
  		xlab("Collection Date") + ylab("Latitude")

### summary stats
  ### DEST summary stat stuff from: DEST/utils/averageReadDepth_GDS.R
  #save(summaryStat.plot, mpsl, file="DEST/figures/figure_1/DEST_summary_plot_data.Rdata")
  load(file="DEST/figures/figure_1/DEST_summary_plot_data.Rdata")

    summaryStat.plot <- ggplot(data=mpsl[(variable=="propMissing" & value<.2) | (variable!="propMissing")], aes(x=xf, y=value, color=continent, fill=continent)) +
            geom_point(pch=21, alpha=.5, size=2) +
            facet_grid(variable~set, scales="free", space="free_x", switch="y") +
            theme_cowplot() +
            theme(axis.text.x = element_blank(), panel.spacing = unit(1, "lines"),
                  strip.placement = "outside", strip.background = element_blank(),
                  legend.position="bottom") +
            xlab("") + ylab("")


### mega plot
  AB_plot <- plot_grid(world.plot, multi_sample, nrow=2, align="v", axis="l" )
  ABC_plot <- plot_grid(AB_plot, summaryStat.plot, ncol=2, rel_widths=c(2,1.2))
  ggsave(ABC_plot, file="/Users/alanbergland/Documents/work/Proposals:Grants/2020.06 NSF-Career/DEST_summary.pdf", h=8, w=14)
