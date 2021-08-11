### scp aob2x@rivanna.hpc.virginia.edu:~/allResid.Rdata ~/.

### libraries
  library(data.table)
  library(ggplot2)
  library(viridis)
  library(patchwork)
  library(wesanderson)

### load data
  setwd("/Users/alanbergland/Documents/")
  load("GitHub/data-paper/Figure10_and_S13_S14/residuals/allResid.Rdata")

### aggregates
  allResid[,af1:=round((n1/max_n1)*40)/40]
  allResid[,af2:=round((n2/max_n2)*40)/40]

  allResid[model_type=="IMbg", model_type:="S+BG"]
  allResid[model_type=="split_mig", model_type:="S"]
  allResid[model_sym=="asymm", model_sym:="AsyM"]
  allResid[model_sym=="symm", model_sym:="SyM"]

  allResid.ag <- allResid[,list(mu=median(resid, na.rm=T)),
                        list(af1, af2, model_type, model_sym)]

  allResid.ag.ag <- allResid.ag[,list(delta=mu[model_type=="S"] - mu[model_type=="S+BG"]), list(af1, af2, model_sym)]

### plot
  cols <- wes_palette("Zissou1", 100, type = "continuous")

  sfs_plot <-
    ggplot(data=allResid.ag,
            aes(x=af1, y=af2, fill=mu)) +
    geom_raster() +
    facet_grid(model_type~model_sym) +
    scale_fill_gradientn(colours = cols) +
    xlab("Pop. 1 allele freq") +
    ylab("Pop. 2 allele freq") +
    labs(fill="Median\nResidual")

  delta_plot <-
    ggplot(data=allResid.ag.ag,
            aes(x=af1, y=af2, fill=delta)) +
    geom_raster() +
    facet_grid(.~model_sym) +
    scale_fill_gradientn(colours = cols) +
    xlab("Pop. 1 allele freq") +
    ylab("Pop. 2 allele freq") +
    labs(fill="Delta\nResidual\n(S)-(S+BG)")


### mega plot
  layout <- "
  AA
  AA
  BB"

  mega_plot <-
  sfs_plot + delta_plot +
  plot_layout(design = layout) +
  plot_annotation(tag_levels = 'A')

  ggsave(mega_plot, file="~/mega_sfs.png")
