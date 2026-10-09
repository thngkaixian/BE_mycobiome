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
  
  # Top taxa
  # ------------------------------------------------------------------------------------------
  top25 <- names(sort(taxa_sums(ps.prop), decreasing=TRUE))[1:min(25,ntaxa(ps.prop))] ; top25
  ps.prop.topN <- prune_taxa(top25, ps.prop)
  ps.prop.topN <- tax_reorder(ps.prop.topN, top25)
  ps.final <- ps.prop.topN 
  
  
# --------------------------------------------------------------------------------------------------------------
# Color palette
# --------------------------------------------------------------------------------------------------------------
  
  ## Genus
  colTaxa <- c("Candida" = "#468A4B" , # "468A4B" green4
    "Saccharomyces" = "#BFDBE5" , 
    "Penicillium" = "#CFDE7F" , 
    "Microidium" = "#7F4F99" ,
    "Aspergillus" = "#E16D53" , 
    "Malassezia" = "#3D898E" , 
    "Cladosporium" = "#9DCAEC" ,
    "Mycosphaerella" = "#EFBFCC" ,
    "Schizophyllum" = "#F3EB52" ,
    "Filobasidium" = "#2F6036" ,
    "Fusarium" = "#BB3630" ,
    "Heterobasidion" = "#F6E2C6" ,
    "Clavispora" = "#EEE592" ,
    "Cutaneotrichosporon"  = "#A8312D" ,
    "Wallemia" = "#665EA4", 
    "Trametes" = "#E5823A", 
    "Cryptococcus" = "#84C9DB", 
    "Lentinus" = "#5873B4" ,
    "Chaenotheca" = "#DE9883", 
    "Nakaseomyces" = "#2B2D76" , 
    "Rhodotorula" = "#7ABA57" )
  colTaxa<- data.frame(colTaxa)
  colnames(colTaxa) <- c("Hex_code")


# ==============================================================================================================
# SUBSET into subgroups
# ==============================================================================================================


  ## By matched cohorts
  # ------------------------------------------------------------------------------------------
  ps.final_AGB <- subset_samples(ps.final, !is.na(Matching_AgeGenderBSI))
  # unique(sample_data(ps.final_AGB)$Matching_AgeGenderBSI)
  
  ps.final_AGEL <- subset_samples(ps.final, !is.na(Matching_AgeGenderExacerbatorStatusFEV1))
  # unique(sample_data(ps.final_AGEL)$Matching_AgeGenderExacerbatorStatusFEV1)
  
  
# --------------------------------------------------------------------------------------------------------------
## ASSIGN ANALYSIS OF INTEREST
# --------------------------------------------------------------------------------------------------------------

  PS_plot <- ps.final_AGB ; COHORT <- "MatchedAGB"
  PS_plot <- ps.final_AGEL ; COHORT <- "MatchedAGEL"

    
# --------------------------------------------------------------------------------------------------------------
# Aggregated barplots - By Sample Source
# --------------------------------------------------------------------------------------------------------------

GROUP_VAR <- "Continent"

  
  #### Aggregate by clinical groups
  # --------------------------------------------------------------------------------------------------------------
  ps.merged <- transform_sample_counts(merge_samples(PS_plot, GROUP_VAR), function(otu) {if (sum(otu)==0) otu else 100*otu/sum(otu)})
  sample_data(ps.merged)[,GROUP_VAR] <- sample_names(ps.merged)
  
  
    ## Factor data and specify labels
    # --------------------------------------------------------------------------------------------------------------
    
    ## --- Continent
    sample_data(ps.merged)[,GROUP_VAR] <- factor(sample_data(ps.merged)$Continent,
                                                 levels = c("Asia", "Europe"),
                                                 labels = c("South-East Asia\n(N=106)", "Europe\n(N=106)"))
    
    
  ## Plot
  # --------------------------------------------------------------------------------------------------------------
  p_Agg <- plot_bar(ps.merged, x=GROUP_VAR, fill="Genus") +
    geom_bar(stat="identity", aes(color=Genus)) +
    xlab(NULL) +
    ylab('Relative abundance (%)') +
    guides(fill= guide_legend(ncol= 1)) +
    labs(fill = "Genus") +  
    # facet_grid(.~Region, scales = "free", space = "free") +
    theme(legend.text = element_text(face = "italic", size = 7),
          legend.key.height = unit(1, "mm"),
          axis.text.x = element_text(color="black", angle = 0, hjust = 0.5),
          axis.text.y = element_text(color="black"),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    scale_fill_manual(values=as.character(colTaxa[taxa_names(ps.merged),])) +
    scale_colour_manual(values=as.character(colTaxa[taxa_names(ps.merged),]))
  p_Agg$data$Genus <- factor(p_Agg$data$Genus, levels= taxa_names(ps.merged))
  p_Agg
  
  #### NO LEGENDS ####
  # --------------------------------------------------------------------------------------------------------------
  p_Agg_NoLegends <- p_Agg + theme(legend.position = "none")
  p_Agg_NoLegends
  
  # Output barplot into .png file
  output_barplot <- file.path(outdir, paste0("Barplot_", taxLevel, '-Aggregate', "-", COHORT, "-", GROUP_VAR, "-NoLegends", '.png'))
  
  ggsave(output_barplot, plot = p_Agg_NoLegends, height = 8, width = 10, units = "cm") # Continent
  
  
    

  
  
  
  
  
  

  
  