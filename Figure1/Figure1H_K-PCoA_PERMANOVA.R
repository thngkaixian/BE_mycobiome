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

  ## Asia | Europe
  # ------------------------------------------------------------------------------------------
  ps.prop_Asia <- subset_samples(ps.prop, Continent == "Asia")
  ps.prop_Europe <- subset_samples(ps.prop, Continent == "Europe")

  ## Regions - NWE, SE, UK 
  # ------------------------------------------------------------------------------------------
  ps.prop_NWE <- subset_samples(ps.prop, Region=="NWE")
  ps.prop_SE <- subset_samples(ps.prop, Region=="SE")
  ps.prop_UK <- subset_samples(ps.prop, Region=="UK")

  
# --------------------------------------------------------------------------------------------------------------
## ASSIGN ANALYSIS OF INTEREST
# --------------------------------------------------------------------------------------------------------------

  PS.PROP <- ps.prop ; COHORT <- "AllPatients"
  
  PS.PROP <- ps.prop_Asia ; COHORT <- "Asia"
  # PS.PROP <- ps.prop_Europe ;  COHORT <- "Europe"
  
  # PS.PROP <- ps.prop_NWE ; COHORT <- "RegionNWE"
  # PS.PROP <- ps.prop_SE ;  COHORT <- "RegionSE"
  # PS.PROP <- ps.prop_UK ;  COHORT <- "RegionUK"
  
  
# --------------------------------------------------------------------------------------------------------------
## PERMANOVA
# --------------------------------------------------------------------------------------------------------------
  
  groupvar <- "Continent"
  # groupvar <- "Region"
  groupvar <- "Country"
  
  
  ## PERMANOVA
  # --------------------------------------------------------------------------------------------------------------
  group <- as.matrix(sample_data(PS.PROP))[, groupvar]
  table(data.frame(sample_data(PS.PROP))[, groupvar], exclude=NULL)
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
  
    
    ## Factor data 
    # --------------------------------------------------------------------------------------------------------------

    ## --- Country
    gg$GroupVar <- factor(gg$GroupVar,
                          levels = c("Singapore", "KualaLumpur", "Belgium", "Germany", "Netherlands", "Greece", "Spain","England", "Scotland")) #,
    
    gg$Region <- sapply(gg$GroupVar, function(x){
      if(x=="Singapore" | x=="KualaLumpur") "Asia"
      else if (x=="Belgium" | x=="Germany" | x=="Netherlands") "Northern and Western Europe"
      else if (x=="Greece" | x=="Spain") "Southern Europe"
      else if (x=="England" | x=="Scotland") "UK"
    })
 
  
    #### Colours ####
    # --------------------------------------------------------------------------------------------------------------
    
    ####===== Continent =====####
    Group_Colours <- c("Asia" = "darkorange", "Europe"="purple4")
    Group_Labels <- c("Asia"="South-East Asia", "Europe"="Europe")
    
    ####===== Region =====####
    Group_Colours <- c("SG_KL" = "darkorange", "NWE" = "darkolivegreen", "SE" = "royalblue", "UK" = "darkred" ) # "SG_KL" = "orange",
    Group_Labels <- c("SG_KL"="Asia", "NWE"="NWE", "SE"="SE", "UK"="UK")
    
    
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

    # Output into .png file
    # p_filename <- file.path(outdir, "00_PCoA_FINAL", paste0("PCoA_", taxLevel, "_", COHORT, "_", groupvar, '.png'))
    # ggsave(p_filename, plot = p, height = 8, width = 10, units = "cm")
    
  
    #### Plot PCOA --- V2 remove legends
    p_NoLegend <- p + theme(legend.position = "none") 
    p_NoLegend
    
    p_filename <- file.path(outdir, paste0("PCoA_", taxLevel, "_", COHORT, "_", groupvar, "-V2-NoLegend", '.png'))
    ggsave(p_filename, plot = p_NoLegend, height = 8, width = 8, units = "cm") # Sqaure
    
    
    
    
    
    
    
    
   