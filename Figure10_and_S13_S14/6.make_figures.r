### This code builds the panels for figures 10 and 11 

#########################
### figure 10
#########################

#### Figure 10A

# frist load the data for this ... this was made using script 5.
#load("/project/berglandlab/moments_jcbn_keric/R_data/AllDataMerged_FromBounds_and_Theta.Rdata")

#merged_datasets %>%
#		.[which(.$AIC_label == "Best"),] %>%
#		mutate(theta_est = theta/L ) ->
#		o.best.caller.test
#
#save(o.best.caller.test, file = "data.for.10A.Rdata")

load("./data.for.10A.Rdata")

	o.best.caller.test %>%
		.[which(.$inference_method == "wide-bounds"),] %>%
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
				xlim(0.0044,0.006) +
				xlab(expression(theta))  ->
	theta_plot_Caller

ggsave(theta_plot_Caller,
		file="theta_plot_Caller.pdf",
		 width = 4,
	     height = 3)


#### Figure 10B

#load("/project/berglandlab/moments_jcbn_keric/R_data/AllDataMerged_FromBounds_and_Theta#.Rdata")
#
#merged_datasets %>%
#		.[which(.$AIC_label == "Best"),] %>%
#		mutate(theta_est = theta/L ) ->
#		o.best.caller.test
#
#	o.best.caller.test %>%
#		.[which(.$inference_method == "wide-bounds"),] %>%
#		.[which(.$caller == "PoolSNP"),] %>%
#		.[,c(
#		"pop1",
#		"pop2",
#		"SFSmethod",
#		"inference_method",
#		"popset",
#		"mij"
#		)] %>%
#	dcast(pop1+pop2+popset ~ SFSmethod, 
#	value.var = "mij",
#	mean ) %>%
#	mutate(stat = "mij") ->
#	Mij_by_caller 

#save(Mij_by_caller, file = "data.for.10B.Rdata")

load("./data.for.10B.Rdata")


Mij_by_caller$binom %>% sd
Mij_by_caller$counts %>% sd

	o.best.caller.test %>%
		.[which(.$inference_method == "wide-bounds"),] %>%
		.[which(.$caller == "PoolSNP"),] %>%
		.[complete.cases(.$popset)] %>%
	ggplot( aes(x=SFSmethod, y=(mij), fill=popset)) +
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
  ylab(expression((m[ij]))) +
  xlab("AF discretization method") +
  theme_classic() +
  theme(legend.pos = "none") +
  scale_fill_brewer(palette = "Pastel1") ->
	Mij_plot_Caller

ggsave(Mij_plot_Caller,
		file="TMij_plot_Caller.pdf",
		 width = 4,
	     height = 3)

#### Figure 10C

# Load final runs for model testing
#load("/project/berglandlab/moments_jcbn_keric/final_runs/AllData_IM_SM_Mig_models.Rdata")
#
#o.best = o.l.all_mapped %>%
#		.[which(.$AIC_label == "Best"),] %>%
#		mutate(theta_est = theta/L )
#
#save(o.best, file = "data.for.10C.and.all11.Rdata")

load("./data.for.10C.and.all11.Rdata")

### How many times was a model conisdered best?
o.best %>%
  group_by(Pair_name) %>%
  summarise(AIC = min(AIC) ) -> min_aic_models


left_join(min_aic_models, o.best) %>% 
  .$full_model %>%
  table() %>%
  prop.table() %>% 
  data.frame(Model = `.`) %>%
  ggplot(
    aes(
      x=`Model..`,
      y=Model.Freq
    )
  ) + 
  geom_bar(stat = "identity") +
  theme_classic() ->
  best_models

ggsave(best_models,
       file = "best_models.pdf")


#### Figure 10D

### Load final runs for model testing
#load("/project/berglandlab/moments_jcbn_keric/final_runs/AllData_IM_SM_Mig_models.Rdata")
#
#o.best = o.l.all_mapped %>%
#		.[which(.$AIC_label == "Best"),] %>%
#		mutate(theta_est = theta/L )
#
#### Do model AIC improve with more iterations?
#### This makes ---> Figure 10 panels D
#### 
#o.l.all_mapped %>% 
#  as.data.frame() %>% 
#  mutate(Sample_size_bin = 
#           ifelse( .$Within_Sim_id <= 9, 
#            print("1-9"),
#            signif(Within_Sim_id, digits=1))) ->
#  o.l.all_mapped_for_bin_plotting
 
#save(o.l.all_mapped_for_bin_plotting, file = "data.for.10D.Rdata")

load("./data.for.10D.Rdata")

        
o.l.all_mapped_for_bin_plotting %>%
  group_by(Pair_name) %>%
  summarise(Min_AC_glob = min(AIC)) -> global_min_AIC



left_join(o.l.all_mapped_for_bin_plotting, global_min_AIC) %>% 
  mutate(full_model = paste(Demo_model, 
                            migration_model, 
                            sep = "_") ) %>%
  mutate(Delta_AIC = as.numeric(AIC)-as.numeric(Min_AC_glob)) %>%
  ggplot(
    aes(
      log10(Delta_AIC+1),
      fill=formal_model_name
    )
  ) + geom_histogram(
    aes(y = ((..count..)/sum(..count..))*100),
    #position = "dodge"
    position = "fill"
    ) +
  ylab(
    expression(paste(
     "Proportion of models at ",
     delta, "AIC bin", sep = ""))) +
  xlab(expression(paste(
    Log[10], delta, "(AIC)", sep = ""))) +
  theme(legend.position = "bottom") +
  labs(fill = "Model") +
  facet_wrap(~as.factor(Sample_size_bin)) ->
  delta_aic

ggsave(delta_aic,
       file = "panel_delta_aic.pdf",
       #width = 4,
       #height = 4
       )

#########################
#########################
#########################
#########################

## Figure 11A

#####
######
###### Plot the divergence time
######
load("./data.for.10C.and.all11.Rdata")

o.best %>%
  .[which(.$Run_Type == "final_runs"),] %>% 
  group_by(Demo_cluster) %>%
  summarise(DT_m = mean(divergence_time),
            DT_low = ci(divergence_time)[2],
            DT_high = ci(divergence_time)[3],
            Median_Dt = median(divergence_time),
            )


o.best %>%
  .[which(.$Run_Type == "final_runs"),] %>% 
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

#ggsave(Divergence_plot,
#       file="Divergence_plot.pdf",
#       width = 5,
#       height = 4)
#

## Figure 11B
load("./data.for.10C.and.all11.Rdata")

#### Divergence as a function of distance
o.best %>%
  .[which(.$Run_Type == "final_runs"),] %>% 
  ggplot(aes(
    y=log10(divergence_time),
    x=(dist),
    fill=Demo_cluster,
    color = Demo_cluster
  )) + 
  #geom_density2d() +
  geom_point(shape = 21,
             alpha = 0.3,
             color = "black") +
  geom_smooth(method = "lm",
              aes(linetype =Demo_cluster),
              color = "black") +
  #facet_wrap(~Demo_cluster, 
  #           scales = "free",
  #           ncol =1) +
  ylab(expression(Log[10]("Divergence Time"))) +
  xlab("Distance x100 Km") +
  theme_bw() +
  theme(legend.position = "none") +
  scale_color_brewer(palette = "Pastel1") ->
  divergence_dist_plot

### Summary stats
cor.test(o.best$dist[which(o.best$full_model == "SM_asymm" & o.best$Demo_cluster == "EE" )],
         log10(o.best$divergence_time[which(o.best$full_model == "SM_asymm" & o.best$Demo_cluster == "EE" )])
)

 cor.test(o.best$dist[which(o.best$full_model == "SM_asymm" & o.best$Demo_cluster == "WW" )],
          log10(o.best$divergence_time[which(o.best$full_model == "SM_asymm" & o.best$Demo_cluster == "WW" )])
 )
 
 cor.test(o.best$dist[which(o.best$full_model == "SM_asymm" & o.best$Demo_cluster == "EW" )],
          log10(o.best$divergence_time[which(o.best$full_model == "SM_asymm" & o.best$Demo_cluster == "EW" )])
 )


## Figure 11C
 load("./data.for.10C.and.all11.Rdata")
 
#### Migration as a function of distance

#Import cluster metatdta
cluster_metadata <- fread("/project/berglandlab/moments_jcbn_keric/final_runs/DEST_freeze1/populationInfo/Cluster_Assingment/DEST_Sample_clusters.txt")

cluster_metadata %>%
  .[,c("sampleId","Continental_clusters")] -> base_met

pop1_met = base_met
names(pop1_met) = c("pop1","Pop1_cluster")

pop2_met = base_met
names(pop2_met) = c("pop2", "Pop2_cluster")

o.best %>%
  .[which(.$Run_Type == "final_runs"),] %>% 
  left_join(.,pop1_met ) %>% 
  left_join(pop2_met) %>% 
  mutate(M_1_to_2 = paste(Pop1_cluster, Pop2_cluster, sep = "_to_")) %>%
  mutate(M_2_to_1 = paste(Pop2_cluster, Pop1_cluster, sep = "_to_")) ->
  o.best.migrationCalibrated

## Get absolute rates of migration
o.best.migrationCalibrated %>% 
  .[which(.$Demo_cluster == "EW"),] %>%
  group_by(M_1_to_2) %>%
  summarise(m12_mean = ci(m12)[1],
            m12_l = ci(m12)[2],
            m12_h = ci(m12)[3],
            migrants12_mean =ci(mig_pop1)[1],
            migrants12_l =ci(mig_pop1)[2],
            migrants12_h =ci(mig_pop1)[3]
            )

o.best.migrationCalibrated %>% 
  .[which(.$Demo_cluster == "EW"),] %>%
  group_by(M_2_to_1) %>%
  summarise(m21_mean = ci(m21)[1],
            m21_l = ci(m21)[2],
            m21_h = ci(m21)[3],
            migrants21_mean =ci(mig_pop2)[1],
            migrants21_l =ci(mig_pop2)[2],
            migrants21_h =ci(mig_pop2)[3]
  )


o.best.migrationCalibrated %>%
  .[which(.$Demo_cluster == "EW"),] %>%
  .[which(.$Run_Type == "final_runs"),] %>% 
  #.[,c("m12", "M_1_to_2", "Demo_cluster", "dist.bin")] %>% 
  group_by(Demo_cluster, dist.bin, M_1_to_2) %>% 
  summarise(Mij_mean = ci(m12)[1],
            Mij_low = ci(m12)[2],
            Mij_hihg = ci(m12)[3],
            Migrants_mean = ci(mig_pop1)[1],
              Migrants_low = ci(mig_pop1)[2],
              Migrants_high = ci(mig_pop1)[3]
            ) -> mean_m12
names(mean_m12)[3] = "direction"  
          
o.best.migrationCalibrated %>%
  .[which(.$Demo_cluster == "EW"),] %>%
  .[which(.$Run_Type == "final_runs"),] %>% 
  #.[,c("m21", "M_2_to_1", "Demo_cluster", "dist.bin")] %>% 
  group_by(Demo_cluster, dist.bin, M_2_to_1) %>% 
  summarise(Mij_mean = ci(m21)[1],
            Mij_low = ci(m21)[2],
            Mij_hihg = ci(m21)[3],
            Migrants_mean = ci(mig_pop2)[1],
            Migrants_low = ci(mig_pop2)[2],
            Migrants_high = ci(mig_pop2)[3]
            ) -> mean_m21
names(mean_m21)[3] = "direction"  
  
 rbind(mean_m12, mean_m21) -> migration_stats

 
 migration_stats$direction = gsub(".Europe", "",  migration_stats$direction)
 migration_stats$direction = gsub("1_", "",  migration_stats$direction)
 migration_stats$direction = gsub("3_", "",  migration_stats$direction)
 migration_stats$direction = gsub("_to_", " to ",  migration_stats$direction)
 
 migration_stats 
 
 migration_stats %>%
  ggplot(aes(
    y = Mij_mean,
    ymin =Mij_low,
    ymax =Mij_hihg,
    x=dist.bin
  )) + 
  geom_errorbar(width = 2) +
  geom_line() +
  geom_point(size = 3, 
             shape = 21,
             fill = "white") +
  theme_bw() +
  facet_grid(~direction) +
  ylab(expression(M[i%->%j])) +
  xlab("Distance bin (Km x 100)") ->
  mij_dist_plot


 ######
 ######
 ######
##################################
##Compound plot
##################################

ggsave(Divergence_plot + divergence_dist_plot + mij_dist_plot,
       file = "Compound_DT_plot.pdf",
       width = 8,
       height = 2.5)

#######
#######
#######
#######

#### Extra Stuff:

##################################
## Get Ne for each cluster
##################################

 cluster_metadata <- fread("/project/berglandlab/moments_jcbn_keric/final_runs/DEST_freeze1/populationInfo/Cluster_Assingment/DEST_Sample_clusters.txt")
 
 cluster_metadata %>%
   .[,c("sampleId","Continental_clusters")] -> base_met
 
  o.best %>% 
   .[which(.$Run_Type == "final_runs"),] %>% 
   .[,c("pop1", "pop1_size")] -> pop1
 names(pop1) = c("sampleId", "Ne")
 
 o.best %>% 
   .[which(.$Run_Type == "final_runs"),] %>% 
   .[,c("pop2", "pop2_size")] -> pop2
 names(pop2) = c("sampleId", "Ne")
 
 rbind(pop1, pop2 ) %>%
   left_join(base_met) %>%
   group_by(Continental_clusters) %>%
    summarise(Ne_mean = ci(Ne)[1],
             Ne_l = ci(Ne)[2],
             Ne_h = ci(Ne)[3]
   )
  
##################################
# Explore pi differences east vs west
##################################

pi_data <- fread("/project/berglandlab/moments_jcbn_keric/final_runs/data-paper/Figure6/Figure6_data.txt")

cluster_metadata <- fread("/project/berglandlab/moments_jcbn_keric/final_runs/DEST_freeze1/populationInfo/Cluster_Assingment/DEST_Sample_clusters.txt")

cluster_metadata %>%
  .[,c("sampleId","Continental_clusters")] -> base_met


names(pi_data)[1] = "sampleId"

left_join(pi_data, base_met) ->
  pi_data_annot

pi_data_annot %>%
  group_by(Continental_clusters) %>%
  summarise(Mean = mean(Pi, na.rm = T))

t.test(
  pi_data_annot$Pi[which(pi_data_annot$Continental_clusters == "3.Europe_E")],
  pi_data_annot$Pi[which(pi_data_annot$Continental_clusters == "1.Europe_W")],
)
