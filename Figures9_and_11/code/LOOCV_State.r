#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(adegenet)
library(tidyverse)

load("./data_to_train_State_DAPC.Rdata")
#dat_filt_maf_LD500_naimp_forStateTrain
#dat_filt_maf_LD500_naimp_Labels_StatesTrain
load("./AIM_SNPs.Rdata")


dat_filt_maf_LD500_naimp_forStateTrain[,which(colnames(dat_filt_maf_LD500_naimp_forStateTrain) %in% AIMS_Subset$SNPid)] -> datSNPs_AIMset

dim(datSNPs_AIMset)[1] == dim(dat_filt_maf_LD500_naimp_Labels_StatesTrain)[1]

xval_tmp <- xvalDapc(
	datSNPs_AIMset[-as.numeric(args[1]),],
	 dat_filt_maf_LD500_naimp_Labels_StatesTrain$State[-as.numeric(args[1])],
 	 n.pca.max = 50,
 	  training.set = 0.9,
 	  result = "groupMean",
 	   center = TRUE, 
 	   scale = FALSE,
 	   n.pca = NULL, 
 	   n.rep = 30, 
 	   xval.plot = FALSE)


model=xval_tmp$DAPC
LOOCV=as.data.frame(t(datSNPs_AIMset[as.numeric(args[1]),]))

pred.sup <- predict.dapc(model, 
                         newdata=LOOCV
                         )

pred.sup$posterior %>% t() %>% as.data.frame -> pred.sup.posterior

top_posterior=slice_max( pred.sup.posterior, V1 ) 

LABEL=dat_filt_maf_LD500_naimp_Labels_StatesTrain$State[as.numeric(args[1])]

label_posterior=pred.sup.posterior$V1[which(rownames(pred.sup.posterior) == LABEL )]

write.table( cbind(rownames(top_posterior),
                   top_posterior,
                   LABEL,
                   label_posterior,
                   as.numeric(args[1]), 
                   rownames(datSNPs_AIMset)[as.numeric(args[1])],
                   model$n.pca,
                   model$n.da), 
             file = "./DAPC_State_LOOVC.txt", 
             sep = "\t",quote = F ,row.names = F, col.names = F, append = T)