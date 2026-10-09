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
  
  ## Relative abundance
  ps.prop <- transform_sample_counts(ps, function(otu) {if (sum(otu)==0) otu else 100*otu/sum(otu)})
  

# ==============================================================================================================
# SUBSET into subgroups
# ==============================================================================================================

  ## By matched cohorts
  # ------------------------------------------------------------------------------------------
  ps.prop_AGB <- subset_samples(ps.prop, !is.na(Matching_AgeGenderBSI))
  # unique(sample_data(ps.prop_AGB)$Matching_AgeGenderBSI)
  
  ps.prop_AGEL <- subset_samples(ps.prop, !is.na(Matching_AgeGenderExacerbatorStatusFEV1))
  # unique(sample_data(ps.prop_AGEL)$Matching_AgeGenderExacerbatorStatusFEV1)


# --------------------------------------------------------------------------------------------------------------
## ASSIGN ANALYSIS OF INTEREST
# --------------------------------------------------------------------------------------------------------------
  
  PS.PROP <- ps.prop_AGB ; COHORT <- "MatchedAGB" #Matched Cohort
  PS.PROP <- ps.prop_AGEL ; COHORT <- "MatchedAGEL"  #Matched Cohort
  

# --------------------------------------------------------------------------------------------------------------
## PERMANOVA
# --------------------------------------------------------------------------------------------------------------

  groupvar <- "Continent"
  
  group <- as.matrix(sample_data(PS.PROP))[, groupvar]
  braydist <- as.matrix(vegdist(otu_table(PS.PROP), "bray"))
  PERMANOVA <- adonis2(braydist ~ group, permutations = 10000)
  pval <- PERMANOVA$`Pr(>F)`[1]
  pval
  
  ## Export p-value
  write.csv(pval, file = file.path(outdir, paste0("PERMANOVA_", taxLevel, "_", COHORT, "_", groupvar, ".csv")))
  
  
  ## PCOA
  # --------------------------------------------------------------------------------------------------------------
  ord.nmds.bray <- ordinate(PS.PROP, method= 'PCoA', distance= 'bray')
  gg <- cbind(ord.nmds.bray$vectors, sample_data(PS.PROP)[,groupvar])
  colnames(gg)[ncol(gg)] <- "GroupVar"
  centroids <- aggregate(cbind(Axis.1,Axis.2)~GroupVar, data=gg, mean)
  gg <- merge(gg,centroids,by="GroupVar",suffixes=c("",".centroid"))
  
  
  #### Colours ####
  # --------------------------------------------------------------------------------------------------------------
  
  ####===== Continent =====####
  Group_Colours <- c("Asia" = "darkorange", "Europe"="purple4")
  Group_Labels <- c("Asia"="South-East Asia", "Europe"="Europe")
  
  
  #### Plot PCOA
  # --------------------------------------------------------------------------------------------------------------
  p <- ggplot(gg) +
    scale_linetype_identity() +
    geom_segment(aes(x=Axis.1.centroid, y=Axis.2.centroid, xend=Axis.1, yend=Axis.2, colour = factor(GroupVar)), 
                 alpha = 0.3) +
    geom_point(aes(x=Axis.1,y=Axis.2, color = factor(GroupVar)), size=2, alpha = 0.5) + 
    geom_point(data=centroids, aes(x=Axis.1, y=Axis.2, color=GroupVar), size=5, shape = 16) +
    geom_point(data=centroids, aes(x=Axis.1, y=Axis.2, color=GroupVar), size=5, shape = 13, colour = "black") +
    scale_color_manual(values = Group_Colours, labels = Group_Labels) + #
    theme_bw() +
    theme(legend.position = "bottom",
          plot.subtitle = element_text(hjust = 1, size=9),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black")) +
    labs(colour="",
         subtitle = ifelse(pval<0.001, "PERMANOVA<0.001", ifelse(pval<0.01, "PERMANOVA<0.01", ifelse(pval<0.05, "PERMANOVA<0.05", "PERMANOVA>0.05"))),
         y = paste0("PC2 (",round(100*ord.nmds.bray$values[2,'Relative_eig'],2),"%)"),
         x = paste0("PC1 (",round(100*ord.nmds.bray$values[1,'Relative_eig'],2),"%)")
    )
  p
  
  
  #### Plot PCOA --- V2 remove legends
  p_NoLegend <- p + theme(legend.position = "none") 
  p_NoLegend
  
  p_filename <- file.path(outdir, paste0("PCoA_", taxLevel, "_", COHORT, "_", groupvar, "-V2-NoLegend", '.png'))
  ggsave(p_filename, plot = p_NoLegend, height = 8, width = 8, units = "cm") # Sqaure
  






