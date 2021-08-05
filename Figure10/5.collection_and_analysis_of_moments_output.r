module load intel/18.0 intelmpi/18.0
module load goolf/7.1.0_3.1.4
module load gdal proj R/4.0.0
R

### libraries
  library(data.table)
  library(SeqArray)
  library(foreach)
  library(sp)
  library(doMC)
  registerDoMC(8)
  library(tidyverse)
  library(magrittr)
  library(gmodels)
  library(FactoMineR)
  library(factoextra)
  library(ggridges)
  library(viridis)
  library(hrbrthemes)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(ggdist)
  library(patchwork)
  
### get files

	paths <- c(
	"/project/berglandlab/moments_jcbn_keric/binombound_thetacount_outputs/output/",
	"/project/berglandlab/moments_jcbn_keric/alan/moments_BoundsCount/output/",
	"/project/berglandlab/moments_jcbn_keric/jcbn_runs/"
	)

o.l <- list()
for(j in 1:length(paths)){

	file_path = paths[j]
	
  ### run loop internal
    fs <- list.files(file_path, full.names=T)
    fsl <- gsub(file_path, "", fs)
    
    fsl %<>%
    	data.frame(files=.) %>% 
    	separate(files, into=c(
    	"caller",
    	"SFSmethod",
    	"pop1",
    	"pop2",
    	"Demo_cluster",
    	"inference_method",
    	"file"
    	), sep ="\\.")
    
    fsl$caller = gsub("/", "",  fsl$caller)
	fsl$Demo_cluster = gsub("_output", "",  fsl$Demo_cluster)
  
    o = list()
    for(i in 1:length(fs)){ #Open i
    
	  ##Read the file
      tmp <- fread(fs[i])
	  ##Add metadata	  	  
      tmp[,caller:=fsl[i,"caller"]]
      tmp[,SFSmethod:=fsl[i,"SFSmethod"]]
      tmp[,pop1:=fsl[i,"pop1"]]
      tmp[,pop2:=fsl[i,"pop2"]]
      tmp[,Demo_cluster:=fsl[i,"Demo_cluster"]]
      tmp[,inference_method:=fsl[i,"inference_method"]]
      
      if(fsl[i,"inference_method"] == "theta"){
      names(tmp)[10] = "theta"
      }
      
	  ##Evaluate models      		
      		##Add model fir
      		min_aic = min(tmp$AIC)
      		tmp %<>%
      			mutate(AIC_label = 
      				 ifelse(.$AIC == min_aic, 
      										"Best",
      										"NotBest"))
    o[[i]] <- tmp
    } # Close i
    o.l[[j]] <- rbindlist(o, fill=T)
} # close j

o.l.all <- rbindlist(o.l, fill=T)


#rm(list = ls())
save(o.l.all, file = "DataFrom4.1.and4.2.Rdata")

####
load("./DataFrom4.1.and4.2.Rdata")

##### Bring data from first run
load("/project/berglandlab/moments_jcbn_keric/moments_all_backconvert.RData")

df[,c(
"Pair_name",
"L",
"pop1_size",
"pop2_size",
"divergence_time",
"mig_pop1",
"mig_pop2",
"theta",
"LL",
"AIC",
"Nref",
"Ts",
"nu1",
"nu2",
"mij",
"m12"
)] -> run1

#Update names
names(run1)[9] = "-2LL_model"

#Expand metadata
run1 %<>%
	separate(Pair_name, remove = F,
			 into = c("caller","SFSmethod","etc1","pop2","pop1"),
			 sep = "\\.")
			 
#Remove etcs			 
run1 = run1[,-4]

run1 %<>%
	mutate(inference_method = "wide-bounds")

#find best models

run1 %>%
	group_by(Pair_name) %>%
	summarize(MinAIC=min(AIC)) ->
	min_aic_run1

left_join(run1, min_aic_run1) %>% 
	mutate(AIC_label = ifelse(.$AIC == .$MinAIC, "Best", "NotBest")) ->
	run1_aic

## Set 2
o.l.all[,c(
"Pair_name",
"L",
"pop1_size",
"pop2_size",
"divergence_time",
"mig_pop1",
"mig_pop2",
"theta",
"nu1",
"nu2",
"Ts",
"m12",
"-2LL_model",
"AIC",
"caller",
"SFSmethod",
"pop1",
"pop2",
"Demo_cluster",
"inference_method",
"AIC_label")] -> run2_3_aic

run2_3_aic$inference_method = gsub("bound", "shallow-bounds", run2_3_aic$inference_method)
run2_3_aic$inference_method = gsub("theta", "theta-prior",  run2_3_aic$inference_method)


##### Merging datasets
names(run1_aic) 
names(run2_3_aic)

run1_aic$mig_pop1 = as.numeric(run1_aic$mig_pop1)
run1_aic$mig_pop2 = as.numeric(run1_aic$mig_pop2)

run2_3_aic$mig_pop1 = as.numeric(run2_3_aic$mig_pop1) 
run2_3_aic$mig_pop2 = as.numeric(run2_3_aic$mig_pop2)

run1_aic$`-2LL_model` = as.numeric(run1_aic$`-2LL_model` )

run2_3_aic$`-2LL_model` = as.numeric(run2_3_aic$`-2LL_model`) 

list_allsets = list()
list_allsets[[1]] = run1_aic
list_allsets[[2]] = run2_3_aic

merged_datasets <- rbindlist(list_allsets, fill=T)


#Add metadata
meatadata <- "/scratch/yey2sn/moments/data-paper/additionalAnalyses/moments/pairs_all.txt"
meatadata_df <- fread(meatadata)

names(meatadata_df)[3:5] = c("caller","pop1","pop2") 

left_join(merged_datasets, meatadata_df) -> merged_datasets

save(merged_datasets, file = "AllDataMerged_FromBounds_and_Theta.Rdata")


################################
#### Start here:
################################

load("./AllDataMerged_FromBounds_and_Theta.Rdata")
###### 
#####
# Select only the best model
o.best = merged_datasets %>%
		.[which(.$AIC_label == "Best"),] %>%
		mutate(theta_est = theta/L )
		
################################
#### Explore theta
################################

	o.best %>%
		.[which(.$inference_method == "theta-prior"),] %>%
		.[,c(
		"pop1",
		"pop2",
		"caller",
		"SFSmethod",
		"inference_method",
		"theta_est"
		)] %>%
	dcast(pop1+pop2+SFSmethod+inference_method ~ caller, 
	value.var = "theta_est",
	mean ) ->
	theta_by_caller 

# Chart
ggplot(theta_by_caller, aes(x=x) ) +
  # Top
  geom_histogram( aes(x = PoolSNP, y = ..density..), fill="#69b3a2" ) +
  geom_label( aes(x=0.0040, y=2000, label="PoolSNP"), color="#69b3a2") +
  # Bottom
  geom_histogram( aes(x = SNAPE, y = -..density..), fill= "#404080") +
  geom_label( aes(x=0.0040, y=-2000, label="SNAPE"), color="#404080") +
	facet_grid(SFSmethod~inference_method) +
		geom_vline(xintercept = (0.005),
				linetype = "dashed",
				color = "red",
				alpha = 0.5) +
				xlim(0.0044,0.006) +
				xlab(expression(theta))  ->
	theta_plot_Caller

ggsave(theta_plot_Caller,
		file="theta_plot_Caller.pdf",
		 width = 4,
	     height = 3)

################################
#### Explore Ts
################################

	o.best %>%
		.[which(.$inference_method == "theta-prior"),] %>%
		.[which(.$caller == "PoolSNP"),] %>%
		.[,c(
		"pop1",
		"pop2",
		"SFSmethod",
		"inference_method",
		"popset",
		"Ts"
		)] %>%
	dcast(pop1+pop2+popset ~ SFSmethod, 
	value.var = "Ts",
	mean ) %>%
	mutate(stat = "Ts") ->
	Ts_by_caller 

Ts_by_caller$binom %>% sd
Ts_by_caller$counts %>% sd


	o.best %>%
		.[which(.$inference_method == "theta-prior"),] %>%
		.[which(.$caller == "PoolSNP"),] %>%
	ggplot( aes(x=SFSmethod, y=-log10(Ts), fill=popset)) +
	 ggdist::stat_halfeye(
    adjust = 1.5, 
    width = .6, 
    .width = 0, 
    justification = -.2, 
    point_colour = NA
  ) + 
  geom_point(
    ## draw horizontal lines instead of points
    shape = 95,
    size = 4,
    alpha = .2
  ) + 
  ylab(expression(-log[10](T[s]))) +
  xlab("AF discretization method") +
  theme_classic() +
  theme(legend.pos = "none") +
  scale_fill_brewer(palette = "Pastel1") ->
	Ts_plot_Caller

ggsave(Ts_plot_Caller,
		file="Ts_plot_Caller.pdf",
		 width = 4,
	     height = 3)
  
##################################
# Explore Divergence
##################################
	o.best %>%
	.[which(.$inference_method == "theta-prior"),] %>%
	.[which(.$caller == "PoolSNP"),] %>%
	.[which(.$SFSmethod == "binom"),] %>%
	ggplot(aes(
	x=popset,
	y=log10(divergence_time),
	fill = as.factor(popset)
		)) +
  ggdist::stat_halfeye(
    adjust = 1.5, 
    width = .6, 
    .width = 0, 
    justification = -.3, 
    point_colour = NA) + 
  geom_point(
    size = 0.9,
    alpha = .1,
    position = position_jitter(
      seed = 1, width = .1)
  )  +
   geom_boxplot(
    size = 0.2,
    width = .25, 
    outlier.shape = NA
  ) + 
  theme_bw() + 
  theme(legend.position = "none") +
  coord_flip() +
  ylab(expression(Log[10]("Divergence Time"))) +
  xlab("Demographic Cluster") +
  scale_fill_brewer(palette = "Pastel1") ->
	Divergence_plot

ggsave(Divergence_plot,
		file="Divergence_plot.pdf",
		 width = 5,
	     height = 4)


#Table
 	o.best %>%
	.[which(.$inference_method == "theta-prior"),] %>%
	.[which(.$caller == "PoolSNP"),] %>%
	.[which(.$SFSmethod == "binom"),] %>%
 	group_by(popset) %>%
 	summarize(Mean_dt = ci(divergence_time)[1],
 			  Low_dt = ci(divergence_time)[2],
 			  High_dt = ci(divergence_time)[3],
 			  Median_dt = median(divergence_time),
 	)
  
##################################
# Explore distance vs divergence
##################################

	o.best %>%
	.[which(.$inference_method == "theta-prior"),] %>%
	.[which(.$caller == "PoolSNP"),] %>%
	.[which(.$SFSmethod == "binom"),] %>%
	ggplot(aes(
	x=dist,
	y=log10(divergence_time),
	fill = popset
		)) +
	geom_point(alpha = 0.3,
				shape = 21) +
	geom_smooth(color = "black",
				method = "lm") +
	theme_bw() +
	theme(legend.position = "none") +
	facet_wrap(~popset, scales = "free_x") +
	ylab(expression(Log[10]("Divergence Time"))) +
	xlab("Distance (Km)") +
	scale_fill_brewer(palette = "Pastel1")  ->
	DT_dist_plot

ggsave(DT_dist_plot,
		file="DT_dist_plot.pdf",
		 width = 9,
	     height = 2.4)

	o.best %>%
	.[which(.$inference_method == "theta-prior"),] %>%
	.[which(.$caller == "PoolSNP"),] %>%
	.[which(.$SFSmethod == "binom"),] -> divergence_dat

cor.test(divergence_dat$dist[which(divergence_dat$popset == "EE")],
	     log10(divergence_dat$divergence_time[which(divergence_dat$popset == "EE")])
	     )

cor.test(divergence_dat$dist[which(divergence_dat$popset == "WW")],
	     log10(divergence_dat$divergence_time[which(divergence_dat$popset == "WW")])
	     )

cor.test(divergence_dat$dist[which(divergence_dat$popset == "EW")],
	     log10(divergence_dat$divergence_time[which(divergence_dat$popset == "EW")])
	     )


##################################
# Compound plots <---------------------------
##################################

ggsave(((theta_plot_Caller+Ts_plot_Caller+Divergence_plot)/DT_dist_plot),
		file="compound_figure.pdf",
		 width = 8,
	     height = 4)

##################################
# Explore world-plots vs divergence
##################################
world <- ne_countries(scale = "medium", returnclass = "sf")

	o.best %>%
	.[which(.$inference_method == "theta-prior"),] %>%
	.[which(.$caller == "PoolSNP"),] %>%
	.[which(.$SFSmethod == "binom"),] %>%
	.[which(.$dist > 1   ),] -> demographic_estimation


demographic_estimation$divergence_time %>% quantile
#        0%          25%          50%          75%         100% 
#    2.71326    555.82104    710.36148    938.03919 312874.39189 
#   0%       25%       50%       75%      100% 
#0.4334914 2.7449350 2.8514794 2.9722210 5.4953700 

demographic_estimation %<>%
	mutate(divergence_bin = ifelse(.$divergence_time < 500 , "(a) >500 y",
							ifelse(.$divergence_time < 1000, "(b) 500-1000 y", 
							ifelse(.$divergence_time > 5000, "(d) >5000 y", 
							"(c) 1000-5000 y"))))


ggplot(data = world) +
  geom_sf(fill= "antiquewhite") +
  coord_sf(xlim =  c(-12, 41.00), 
  			ylim = c(32.00, 63.00), 
  			expand = FALSE) + 
  			theme(panel.grid.major = element_line(color = gray(.5), 
  			linetype = "dashed", size = 0.5), 
  			panel.background = element_rect(fill = "aliceblue")) + 
  			geom_segment(data=demographic_estimation,
             aes(x = long.V1 , 
             	 y = lat.V1, 
             	 xend = long.V2 , 
             	 yend = lat.V2,
             	 color = divergence_bin), 
             alpha = 0.9,                  
             size = 0.6,
             lineend = "round") +
             theme(legend.position = "none") +
             ylab("Latitude") +
             xlab("Longitude") +
			facet_grid(divergence_bin~popset) ->
			Divergence_map

ggsave(Divergence_map,
		file="Divergence_map.net.pdf",
		 width = 10,
	     height = 9)
ggsave(Divergence_map,
		file="Divergence_map.net.png",
		 width = 10,
	     height = 9)

##################################
# Explore distance vs migration
##################################

	o.best %>%
	.[which(.$inference_method == "theta-prior"),] %>%
	.[which(.$caller == "PoolSNP"),] %>%
	.[which(.$SFSmethod == "binom"),] %>%
	.[,c("dist.bin","mig_pop1","mig_pop2","popset")] %>%
	melt(id=c("dist.bin","popset")) %>% 
	ggplot(aes(
	x=as.factor(dist.bin),
	y=value,
	linetype = variable,
	color = popset
		)) +
	geom_boxplot() +
	#geom_smooth(
	#			method = "lm") +
	theme_bw() +
	theme(legend.position = "none") +
	#facet_wrap(~popset, scales = "free_x") +
	#ylab("mig_pop1") +
	xlab("Distance (Km)") +
	scale_color_brewer(palette = "Pastel1")  ->
	Mighration_dist_plot

ggsave(Mighration_dist_plot,
		file="Mighration_dist_plot.pdf",
		 width = 3,
	     height = 3)
