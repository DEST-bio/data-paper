# This script correspond to the filter 

#ijob -A jcbnunez -c 10 --mem=60G  --partition=standard
#module load gcc/7.1.0  
#module load openmpi/3.1.4
#module load gdal
#module load proj
#module load htslib/1.9
#module load bcftools/1.9
#module load intel/18.0  
#module load intelmpi/18.0
#module load R/4.0.0
#R

#Clean the environment
rm(list = ls())

##### libraries
#install.packages(c("data.table",
#                   "sp",
#                   "foreach",
#                   "tidyverse",
#                   "magrittr",
#                   "patchwork"))

library(data.table)
library(sp)
library(foreach)
library(tidyverse)
library(magrittr)
library(patchwork)

#User defined imput
priv.dt <- fread("SNAPE.goodSamps.0.001.delim") #file
ind.filter <- "good.samps" #good.samps, all.samps
caller <- "SNAPE" #PoolSNPs, SNAPE
maf_thresh <- 0.001 #0.001, 0.05

# "PoolSNPs good.samps 0.001" [Done]
# "PoolSNPs all.samps 0.001"  [Done]
# "PoolSNPs good.samps 0.05"  [Done]
# "PoolSNPs all.samps 0.05"   [Done]

# "SNAPE good.samps 0.001"    
# "SNAPE all.samps 0.001"    
# "SNAPE good.samps 0.05"     
# "SNAPE all.samps 0.05"     

##########################
#data <- list(
#  callers= c( "PoolSNPs", "SNAPE" ),
#  filters= c( "good.samps", "all.samps"),
#  maf_thresholds= c( 0.001, 0.05 )
#)
#
#data %>%
#  cross() %>%
#  map(lift(paste)) %>% 
#  unlist() %>% 
#  unique() 
##########################

if(caller=="PoolSNP") {
  setnames(priv.dt, names(priv.dt), paste("V", 2:9, sep=""))
}
priv.dt[,V8:=paste(V9, paste(V6, V7, V8, sep=""), sep=";")]

## Add names to the prov object
names(priv.dt) = c(
  "AF_bin_aob",
  "chr",
  "pos",
  "nPop",
  "missingPop",
  "ref",
  "alt",
  "afs_mut_pop",
  "afs_zero"
)

## Restrict any further analyses to the 4 main chromosomes				
priv.dt = priv.dt[which(priv.dt$chr %in% c("2L","2R","3L","3R")),]					

#Parse the allele frequency data 
unlist_vect = lapply( priv.dt[,9] , unlist )
split_vect = lapply( unlist_vect , strsplit, split="\\+" )
remove_1stobj = lapply( split_vect$afs_zero, function(x) x[-1])
num_vect = lapply( remove_1stobj , as.numeric )
mean_AF = lapply( num_vect , mean )
AF_vector = as.data.frame(do.call("rbind", mean_AF))

# Add the allele frequency bin vector to the priv object
priv.dt %<>% 
  mutate(AF = round(AF_vector$V1,2) ) %>% 
  .[which(.$AF >= maf_thresh & .$AF <= (1-maf_thresh) ),] %>% 
  mutate(Raster_data = paste(chr,nPop,AF, sep = "_"))

as.data.frame(table(priv.dt$Raster_data  , dnn = list("Raster_data")), responseName = "Freq") %>% 
  separate(Raster_data, into = c("chr","nPop","AF"), sep = "_") -> o

o[,c("nPop","AF","Freq")] = sapply(o[,c("nPop","AF","Freq")], as.numeric)

o %<>% mutate(caller = caller,
              MAF = maf_thresh,
              ind.filter = ind.filter)

save(o, file= paste(caller, 
                    maf_thresh, 
                    ind.filter, 
                    "Rdata", 
                    sep = "."))
