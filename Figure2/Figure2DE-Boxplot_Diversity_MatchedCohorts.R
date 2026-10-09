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

outdir <- file.path(wkdir, "FINAL_ANALYSIS", "Analysis", "Figure_2")
if (!dir.exists(outdir)) try(dir.create(outdir, recursive= TRUE), silent= TRUE)


# --------------------------------------------------------------------------------------------------------------
## Read ps object
# --------------------------------------------------------------------------------------------------------------

ps <- readRDS(file.path(indir_ps, paste0("Final_psITS_", taxLevel, ".RData")))
ps <- prune_taxa(taxa_names(ps)[!grepl('unidentified',taxa_names(ps))], ps)

  ## Filter 
  RelAbundThreshold <- 0.01
  sampleThreshold <- 0.05
  
  ps.prop <- transform_sample_counts(ps, function(otu) {if (sum(otu)==0) otu else otu/sum(otu)})
  tax <- taxa_names(ps.prop)[apply(otu_table(ps.prop), 2, function(x) {sum(x >= RelAbundThreshold) >= sampleThreshold*nsamples(ps)})]
  ps <- prune_taxa(tax, ps)
  
  
  ## Add diversity information to sample data information
  for (method in c('shannon', 'simpson', 'invsimpson')) {
    sample_data(ps)[, paste0(method, 'Diversity')] <- diversity(otu_table(ps), index= method)
  }
  

# ==============================================================================================================
# SUBSET into subgroups
# ==============================================================================================================

  ## By matched cohorts
  # ------------------------------------------------------------------------------------------
  ps_AGB <- subset_samples(ps, !is.na(Matching_AgeGenderBSI))
  # unique(sample_data(ps_AGB)$Matching_AgeGenderBSI)
  
  ps_AGEL <- subset_samples(ps, !is.na(Matching_AgeGenderExacerbatorStatusFEV1))
  # unique(sample_data(ps_AGEL)$Matching_AgeGenderExacerbatorStatusFEV1)

  
# --------------------------------------------------------------------------------------------------------------
## ASSIGN ANALYSIS OF INTEREST
# --------------------------------------------------------------------------------------------------------------

  PS <- ps_AGB ; COHORT <- "MatchedAGB" #Matched Cohort
  # PS <- ps_AGEL ; COHORT <- "MatchedAGEL"  #Matched Cohort
  
  
# --------------------------------------------------------------------------------------------------------------
## DIVERSITY
# --------------------------------------------------------------------------------------------------------------
  
  diversity_var <- "shannonDiversity"
  groupvar <- "Continent"

  
  df_div <- data.frame(sample_data(PS)[, c(groupvar, diversity_var)]) 
  colnames(df_div) <- c("Group", "aDiversity") # "SDI"
  df_div <- subset(df_div, !is.na(Group))

  
  ## Factor data and specify labels
  # --------------------------------------------------------------------------------------------------------------
  
  ## --- Continent
  df_div$Group <- factor(df_div$Group, levels = c("Asia", "Europe"))

  
  #### Colours ####
  # --------------------------------------------------------------------------------------------------------------
  
  ####===== Continent =====####
  Group_Colours <- c("Asia" = "darkorange", "Europe"="purple4")
  Group_Labels <- c("Asia"="South-East Asia", "Europe"="Europe")
  
  
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
  
  
  
  
  
  
  