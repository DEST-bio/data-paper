#module load intel/18.0 intelmpi/18.0
#module load goolf/7.1.0_3.1.4
#module load gdal proj R/4.0.0
#R

### this code makes figures 10 -- Figure 11 and S13

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
### this first part of the code seeks to collate all runs of moments into a single object for analysis
### 

##
##  include the path to all ooutput files from the 4.0 code. In our case we tested various model, a final run with the best model, we also did a barrage of tests to compare poolSNP and SNAPE. hence why the paths are a bit messy....
  
  paths <- c(
	"/project/berglandlab/moments_jcbn_keric/Model_test/output",
	"/project/berglandlab/moments_jcbn_keric/final_runs/output"
	)

# This vector contains the labels of the analysis above. the first path are the runs made for testing models, the second path are the final runs, and all other are tests of PoolSNP vs SNAPE under various priors.
  
	path_meaning <- c(
	"Model_test",
	"final_runs"
	)

## Part 1. Use loop-list approach to load all runs into a single object		
o.l <- list()
for(j in 1:length(paths)){

	file_path = paths[j]
	
  ### run loop internal
    fs <- list.files(file_path, full.names=T)
    fsl <- gsub(file_path, "", fs)
    fsl <- gsub("bound.txt", "SM.bound.symm.txt", fsl)
    fsl <- gsub("bound.asymm.txt", "SM.bound.asymm.txt", fsl)
    fsl <- gsub("widebounds.txt", "widebounds.symm.txt", fsl)
    fsl <- gsub("widebounds", "bound", fsl)
    
        
    fsl %<>%
    	data.frame(files=.) %>% 
    	separate(files, into=c(
    	"caller",
    	"SFSmethod",
    	"pop1",
    	"pop2",
    	"Demo_cluster",
    	"Demo_model", 
    	"inference_method",
    	"migration_model",
    	"file"
    	), sep ="\\.")
    
    fsl$caller = gsub("/", "",  fsl$caller)
	  fsl$Demo_cluster = gsub("_output", "",  fsl$Demo_cluster)
  
    o = list()
    for(i in 1:length(fs)){ #Open i
    
	  ##Read the file
      tmp <- fread(fs[i])
      
      dim(tmp)[1] -> num_o_runs
      
	  ##Add metadata	  	  
      tmp[,caller:=fsl[i,"caller"]]
      tmp[,SFSmethod:=fsl[i,"SFSmethod"]]
      tmp[,pop1:=fsl[i,"pop1"]]
      tmp[,pop2:=fsl[i,"pop2"]]
      tmp[,Demo_cluster:=fsl[i,"Demo_cluster"]]
      tmp[,Demo_model:=fsl[i,"Demo_model"]] # <---- only activate for IM models
      tmp[,inference_method:=fsl[i,"inference_method"]]
      tmp[,migration_model:=fsl[i,"migration_model"]]
      tmp[,Chains_ran:=num_o_runs]
      tmp[,Run_Type:=path_meaning[j]]
      
      #if(fsl[i,"inference_method"] == "theta"){
      #names(tmp)[10] = "theta"
      #}
      
     aic_col = which(names(tmp) == "aic" ) 
     names(tmp)[aic_col] = "AIC" 
     

     ll_col = which(names(tmp) == "ll_model" ) 
     names(tmp)[ll_col] = "-2LL_model" 

     S_col = which(names(tmp) == "fs_sanitycheck" ) 
     names(tmp)[S_col] = "S" 
     
     File_col = which(names(tmp) == "fs_name" ) 
     names(tmp)[File_col] = "fs_file" 
     
     
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


#Add metadata
meatadata <- "/scratch/yey2sn/moments/data-paper/additionalAnalyses/moments/pairs_all.txt"
meatadata_df <- fread(meatadata)

names(meatadata_df)[3:5] = c("caller","pop1","pop2") 

left_join(o.l.all , meatadata_df) -> o.l.all_mapped

o.l.all_mapped %<>% 
  mutate(full_model = paste(Demo_model, 
                            migration_model, 
                            sep = "_") ) %>%
  mutate(Iteration_model = paste(full_model, 
                            Pair_name, 
                            sep = "_") )
  
o.l.all_mapped %<>% 
  group_by(Iteration_model) %>% 
  mutate(Within_Sim_id = row_number()) 

o.l.all_mapped %<>%
  mutate(formal_model_name = NA)


o.l.all_mapped$formal_model_name[o.l.all_mapped$full_model %in% c("IM_asymm")] = "1.S+BG+AsyM"
o.l.all_mapped$formal_model_name[o.l.all_mapped$full_model %in% c("IM_symm")] = "2.S+BG+SyM"
o.l.all_mapped$formal_model_name[o.l.all_mapped$full_model %in% c("SM_asymm")] = "3.S+AsyM"
o.l.all_mapped$formal_model_name[o.l.all_mapped$full_model %in% c("SM_symm")] = "4.S+SyM"


save(o.l.all_mapped, file = "AllData_IM_SM_Mig_models.Rdata")
