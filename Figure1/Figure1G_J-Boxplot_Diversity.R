library(phyloseq)
library(microViz)
library(microbiome)
library(RColorBrewer)
library(vegan)
library(ggplot2)


# --------------------------------------------------------------------------------------------------------------
## Parameters and directories 
# --------------------------------------------------------------------------------------------------------------

project_name <- "MYCOBIOME_AsiaBRIDGE"
taxLevel <- "Genus"

## Directories 
wkdir <- file.path("/Users/thngkaixian/PhD_RESEARCH/RESEARCH_Projects", project_name, "SequencingData")
indir_ps <- file.path(wkdir, "FINAL_ANALYSIS")

outdir <- file.path(wkdir, "FINAL_ANALYSIS", "Analysis", "Figure_1")
if (!dir.exists(outdir)) try(dir.create(outdir, recursive= TRUE), silent= TRUE)


# --------------------------------------------------------------------------------------------------------------
## Read ps object
# --------------------------------------------------------------------------------------------------------------

ps <- readRDS(file.path(indir_ps, paste0("Final_psITS_", taxLevel, ".RData")))


# ==============================================================================================================
# SUBSET into subgroups
# ==============================================================================================================

  ## Asia | Europe
  # ------------------------------------------------------------------------------------------
  ps_Asia <- subset_samples(ps, Continent == "Asia")
  ps_Europe <- subset_samples(ps, Continent == "Europe")

  ## Regions - NWE, SE, UK 
  # ------------------------------------------------------------------------------------------
  ps.prop_NWE <- subset_samples(ps.prop, Region=="NWE")
  ps.prop_SE <- subset_samples(ps.prop, Region=="SE")
  ps.prop_UK <- subset_samples(ps.prop, Region=="UK")

  
# --------------------------------------------------------------------------------------------------------------
## ASSIGN ANALYSIS OF INTEREST
# --------------------------------------------------------------------------------------------------------------

  PS <- ps ; COHORT <- "AllPatients"
  
  # PS.PROP <- ps.prop_Asia ; COHORT <- "Asia"
  PS <- ps_Europe ;  COHORT <- "Europe"
  
  # PS.PROP <- ps.prop_NWE ; COHORT <- "RegionNWE"
  # PS.PROP <- ps.prop_SE ;  COHORT <- "RegionSE"
  # PS.PROP <- ps.prop_UK ;  COHORT <- "RegionUK"
  
  
# --------------------------------------------------------------------------------------------------------------
## DIVERSITY
# --------------------------------------------------------------------------------------------------------------
  
  diversity_var <- "shannonDiversity"
  
  groupvar <- "Continent"
  groupvar <- "Region"
  # groupvar <- "Country"
  
  
  df_div <- data.frame(sample_data(PS)[, c(groupvar, diversity_var)]) 
  colnames(df_div) <- c("Group", "aDiversity") # "SDI"
  df_div <- subset(df_div, !is.na(Group))

  
  ## Factor data and specify labels
  # --------------------------------------------------------------------------------------------------------------
  
  ## --- Continent
  df_div$Group <- factor(df_div$Group, levels = c("Asia", "Europe"))

  ## --- Region
  df_div$Group <- factor(df_div$Group, levels = c("NWE", "SE", "UK"))
  
  
  #### Colours ####
  # --------------------------------------------------------------------------------------------------------------
  
  ####===== Continent =====####
  Group_Colours <- c("Asia" = "darkorange", "Europe"="purple4")
  Group_Labels <- c("Asia"="South-East Asia", "Europe"="Europe")
  
  ####===== Region =====####
  Group_Colours <- c("NWE" = "darkolivegreen", "SE" = "royalblue", "UK" = "darkred" ) ; Group_Labels <- c("NWE"="Northern and\nWestern Europe", "SE"="Southern\nEurope", "UK"="United\nKingdom")
  

  
  ## Boxplot
  # --------------------------------------------------------------------------------------------------------------
  p_Aes <- ggplot(df_div, aes(x=df_div[ , "Group"], y=df_div[ , "aDiversity"], color=Group)) + #, color=Group
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(shape=16, position=position_jitter(0.1)) +
    labs(x=NULL,
      y = "Shannon Diversity Index") +
    ylim(0,2.5) +
    scale_colour_manual(values=Group_Colours) +
    scale_x_discrete(labels = Group_Labels) +
    theme(legend.position = "none",
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  p_Aes
  
  
  plot_outfile <- file.path(outdir, paste0("aDIVERSITY_", taxLevel, "_", COHORT, "_", groupvar, "-", diversity_var, ".png") )
  ggsave(plot_outfile, plot = p_Aes, height = 8, width = 8, units = "cm") 
  
  
  
  
  
  
  