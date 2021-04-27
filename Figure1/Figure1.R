### Make basic plots for DEST samples
### Final data has columns: sampleId, country, city, collectionDate, lat, long, season, nFlies, locality, type (inbred/pooled), continent
### Alan Bergland, Oct 3, 2018; updated Feb 2020


### ijob -c1 -p standard -A berglandlab
### module load gcc/7.1.0  openmpi/3.1.4 R/3.6.0; R


### libraries
	library(data.table)
	library(gdata)
	library(cowplot)
	library(data.table)
	library(foreach)
	library(ggplot2)
	library(ggmap)
	library(maps)
	library(mapdata)
	library(ggthemes)
	library(patchwork)

### set working directory
	#setwd("/scratch/aob2x/dest")
	setwd("/Users/alanbergland/Documents/GitHub/")

### load meta data
  samps <- fread("./DEST_freeze1/populationInfo/samps.csv")

### load sequencing summary stats
	pcr <- fread("./DEST_freeze1/populationInfo/sequencingStats/pcr.csv")
	rd <- fread("./DEST_freeze1/populationInfo/sequencingStats/rd.csv")[auto==T]
	simulans <- fread("./DEST_freeze1/populationInfo/sequencingStats/simulans.csv")[auto==T]

	setkey(rd, sampleId)
	setkey(pcr, sampleId)
	setkey(simulans, sampleId)
	mps <- merge(rd, pcr)
	mps <- merge(mps, simulans)
	### quick fix for Ukrainian samples
		mps[sampleId=="UA_Pir_14_26", sampleId:="UA_Pyr_14_26"]
		mps[sampleId=="UA_Pir_15_21", sampleId:="UA_Pyr_15_21"]
		mps[sampleId=="UA_Pyr_16_48", sampleId:="UA_Pir_16_48"]

	mps <- merge(mps, samps, all.y=T, by="sampleId")
	setnames(mps, "mu.25", "AveReadDepth.25")

	### a few small fixes
	  mps[continent=="North_America", continent:="NorthAmerica"]

	  mps[,effRD.25:=(AveReadDepth.25 * 2*nFlies) / (AveReadDepth.25 + 2*nFlies)]

	### rank x-axis to mean read depth
	  mps <- mps[auto.x==T]
	  mps[,x:=rank(AveReadDepth.25, ties.method="first")]
	  mps[,x.id:=factor(sampleId, levels=mps$sampleId[mps$x])]

	  mps[,propMissing:=nmissing.25/10000]

		setnames(mps, "AveReadDepth.25", "AveReadDepth")
		setnames(mps, "auto.x", "auto")
		setnames(mps, "effRD.25", "effRD")

	### wide to long
	  mpsl <- melt(mps,
	              id.vars=c("x", "sampleId", "continent", "auto", "set"),
	              measure.vars=c("AveReadDepth", "propMissing", "nFlies", "effRD", "pcrDup", "propSimNorm"))

	  mpsl[,xf:=as.factor(x)]






### time plot
	### find sites with multiple time points
		samps.ag <- samps[,list(nSamps=length(locality),
							nSpring=sum(season=="spring"),
							nFall=sum(season=="fall"),
							nTime=length(unique(collectionDate)),
							maxDelta=max(yday) - min(yday),
							lat=mean(lat),
							long=mean(long)),
					list(locality, year, continent, set)]

		setkey(samps.ag, locality, year)
		setkey(samps, locality, year)


### world map plot

	world <- as.data.table(map_data("world"))

	samps.ag.ag <- samps.ag[,list(n=sum(nTime), lat=mean(lat), long=mean(long)), list(locality, set)]
	samps.ag.ag[,set:=factor(set, levels=c("DrosEU", "DrosRTEC", "dgn"))]

### make maps

	min.lat.eu <- 35
	max.lat.eu <- 55
	min.long.eu <- -10
	max.long.eu <- 37
	# [long>=min.long.eu & long<= max.long.eu][lat>=min.lat.eu & lat<=max.lat.eu]
	#[longitude>=min.long.eu & longitude<= max.long.eu][latitude>=min.lat.eu & latitude<=max.lat.eu]


	### summary plot
	  summaryStat.plot <- ggplot(data=mpsl, aes(x=xf, y=value, color=continent, fill=continent)) +
	  geom_point(pch=21, alpha=.95, size=.75) +
	  facet_grid(variable~set, scales="free", space="free_x") +
	  theme_few() + scale_color_tableau() +
	  theme(axis.text.x = element_blank(),
					panel.spacing = unit(.25, "lines"),
					axis.ticks.x = element_blank(),
					legend.position="bottom",
				  legend.text=element_text(size=8),
					axis.text.y = element_text(size=8)) +
	  xlab("Population") + ylab("")

	### plot multi-sample populations

		multi_sample <- ggplot() +
		geom_line(data= samps[J(samps.ag[maxDelta>10])], aes(x=as.Date(yday, origin = as.Date("2018-01-01")), y=lat, group=locality, linetype=continent)) +
		geom_point(data=samps[J(samps.ag[maxDelta>10])], aes(x=as.Date(yday, origin = as.Date("2018-01-01")), y=lat, group=locality)) +
		facet_grid(.~year) +
		scale_x_date(date_labels = "%b", limits = as.Date(c(110,355), origin = as.Date("2018-01-01"))) +
		xlab("Collection Date") + ylab("Latitude") +
		theme_bw() + scale_colour_colorblind() +
		theme(axis.text = element_text(angle = 45, hjust = 1, size=8),
					legend.text=element_text(size=8))


	 world.plot <- 	ggplot() +
	 geom_polygon(data = world,
	 			aes(x=long, y = lat, group = group), fill="lightgrey") +
	 geom_point(data = samps.ag.ag,
	 			aes(x=long, y=lat, size=I((n-1)/5 + 1), color=set), alpha=.75) +
	 xlab("Longitude") + ylab("Latitude") + scale_fill_manual(values="black") +
	 theme_bw() + scale_colour_colorblind() +
	 theme(axis.text = element_text(angle = 45, hjust = 1, size=8),
			   legend.text=element_text(size=8))



	#ggsave(world.plot, file="./DEST/populationInfo/worldPlot.pdf", height=4, width=6)



### big plot
layout <- "
AAAAACC
AAAAACC
AAAAACC
BBBBBCC
BBBBBCC
BBBBBCC"

bigplot <- world.plot + multi_sample + summaryStat.plot +
plot_annotation(tag_levels = 'A', theme=theme(plot.tag = element_text(size = 19))) +
plot_layout(design = layout)



ggsave(bigplot, file="~/Figure3.pdf", h=7, w=9)
