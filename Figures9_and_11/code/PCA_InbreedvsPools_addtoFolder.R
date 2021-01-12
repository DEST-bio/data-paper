#module load gcc/7.1.0  
#module load openmpi/3.1.4
#module load gdal
#module load proj
#module load goolf R/4.0.0
#R

library(adegenet)
library(tidyverse)
library(magrittr)
library(reshape2)
library(vroom)
library(reshape2)
library(FactoMineR)
library(factoextra)
library(vcfR)
library(patchwork)
library(reshape2)
library(zoo)
library(matrixStats)
library(data.table)
library(SeqArray)
library(LEA)
library(gdsfmt)
library(SNPRelate)
library(patchwork)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggExtra)
library(vcfR)

load("/Users/jcbnunez/Downloads/DEST_DGN_AllSNPs_Metadata.Rdata")

DEST_DGN_metadata$country = gsub("USA","United States", DEST_DGN_metadata$country)

dat_filt_maf_LD500_naimp[which(DEST_DGN_metadata$Continental_clusters == "North_America" & DEST_DGN_metadata$set %in% c("dgn","DrosRTEC") ),] %>%  PCA(scale.unit = F, graph = F, ncp = 50) -> LD500_naimp_PCA_50PCs_object

LD500_naimp_PCA_50PCs_object$ind$coord %>% as.data.frame() %>% mutate(sampleId = rownames(.)) %>% left_join(., DEST_DGN_metadata) -> PCA_coords_metadata

PC1_PVE=round(LD500_naimp_PCA_50PCs_object$eig[1,2],2)
PC2_PVE=round(LD500_naimp_PCA_50PCs_object$eig[2,2],2)
PC3_PVE=round(LD500_naimp_PCA_50PCs_object$eig[3,2],2)


#Graph
label.df_PCA <- PCA_coords_metadata %>% 
  separate(locality, into =c("state_code", "city_abrv")) %>%
  group_by(state_code) %>% 
  summarize(Dim.1 = mean(Dim.1), Dim.2 = mean(Dim.2), Dim.3 = mean(Dim.3), Dim.4 = mean(Dim.4) ) 

ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_point(data =PCA_coords_metadata ,  
             aes(x=Dim.1, 
                 y=Dim.2, 
                 size = as.factor(type), 
                 color = as.factor(type), 
                 shape = as.factor(type)), 
                 alpha = 0.8)  +
      theme_classic() + 
      theme(legend.position = "bottom") + 
      scale_shape_manual(values = 16:17) +
      scale_color_manual(values = c("firebrick3","gold3")) +
      xlab(paste(PC1_PVE,"%")) + ylab(paste(PC2_PVE,"%")) +
      scale_size_manual(values = c(5,2)) +  
      ggrepel::geom_label_repel(data = label.df_PCA, size = 1.7 ,  alpha = 0.7, aes(x=Dim.1, y=Dim.2, ,label = state_code)) -> SET1


ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_point(data =PCA_coords_metadata,  
             aes(x=Dim.1, 
                 y=Dim.3, 
                 size = as.factor(type), 
                 color = as.factor(type), 
                 shape = as.factor(type)), 
             alpha = 0.8)  +
  theme_classic() + 
  theme(legend.position = "none") + 
  scale_shape_manual(values = 16:17) +
  scale_color_manual(values = c("firebrick3","gold3")) +
  xlab(paste(PC1_PVE,"%")) + ylab(paste(PC2_PVE,"%")) +
  scale_size_manual(values = c(5,2))  + 
  ggrepel::geom_label_repel(data = label.df_PCA, size = 1.7 ,  alpha = 0.7, aes(x=Dim.1, y=Dim.3, ,label = state_code)) -> SET2


ggplot() + geom_vline(xintercept = 0, linetype = "dashed") + geom_hline(yintercept = 0, linetype = "dashed") + geom_point(data =PCA_coords_metadata ,  aes(x=Dim.1, y=Dim.2, fill = as.factor(country), shape = as.factor(Continental_clusters)),size = 3,  alpha = 0.8)  +  ggrepel::geom_label_repel(data = label.df_PCA, size = 1.7 ,  alpha = 0.7, aes(x=Dim.1, y=Dim.2, ,label = state_code))  + theme_classic() + theme(legend.position = "none") + xlab(paste(PC1_PVE,"%")) + ylab(paste(PC2_PVE,"%")) + scale_shape_manual(values = c(21,22,23,24)) -> DEST_Woldwide_PCA12

ggplot() + geom_vline(xintercept = 0, linetype = "dashed") + geom_hline(yintercept = 0, linetype = "dashed") + geom_point(data =PCA_coords_metadata ,  aes(x=Dim.1, y=Dim.3, fill = as.factor(country), shape = as.factor(Continental_clusters)),size = 3,  alpha = 0.8)  +  ggrepel::geom_label_repel(data = label.df_PCA, size = 1.7 ,  alpha = 0.7, aes(x=Dim.1, y=Dim.3, ,label = state_code))  + theme_classic() + theme(legend.position = "none") + xlab(paste(paste(PC1_PVE,"%"))) + ylab(paste(PC3_PVE,"%")) + scale_shape_manual(values = c(21,22,23,24))  -> DEST_Woldwide_PCA13

#PLOT:
(DEST_Woldwide_PCA12+DEST_Woldwide_PCA13)/(SET1+SET2)

