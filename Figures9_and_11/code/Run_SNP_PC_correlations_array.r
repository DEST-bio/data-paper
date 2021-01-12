#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(FactoMineR)
library(factoextra)
library(tidyverse)

print(args[1])

load("./LD500_naimp_PCA_40PCs_object.Rdata")

temp_dimdesc = dimdesc(LD500_naimp_PCA_40PCs_object, axes = as.numeric(args[1]), proba = 0.05)

temp_dimdesc[[c( paste("Dim",args[1], sep =".") )]] %>% .$quanti %>% as.data.frame() %>% mutate(snpID = rownames(.), PC = as.numeric(args[1])) -> corr_temp

write.table( corr_temp ,
             file = paste("PC",args[1],"correlations.txt", sep = "."),
             sep = "\t",quote = F ,row.names = F, col.names = F, append = T)