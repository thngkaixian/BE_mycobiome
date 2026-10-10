library(phyloseq)
library(vegan)
library(ggplot2)
library(reshape2)


# ------------------------------------------------------------------------------------------------------------------------
# Parameters and working directory
# ------------------------------------------------------------------------------------------------------------------------

taxLevel_Microbiome <- "Species"

wkdir <- "/Users/thngkaixian/PhD_RESEARCH/RESEARCH_Projects"
indir_ps <- file.path(wkdir, "MYCOBIOME_AsiaBRIDGE", "SequencingData", "FINAL_ANALYSIS")

## Out directories
outdir <- file.path(indir_ps, "Analysis", "Figure_5", "PLOTTING-Aspergillus-Species-Abundance-across-SampleTypes", taxLevel_Microbiome)
if (!dir.exists(outdir)) try(dir.create(outdir, recursive= TRUE), silent= TRUE)


# --------------------------------------------------------------------------------------------------------------
## Import data: HHP data
# --------------------------------------------------------------------------------------------------------------  

ps_microbiome_HHP <- readRDS(file.path(indir_ps, paste0("Final_ps_", taxLevel_Microbiome, "_HHP_35patients.RData")))

  ps_microbiome_HHP <- prune_taxa(taxa_names(ps_microbiome_HHP)[apply(otu_table(ps_microbiome_HHP), 2, sum) != 0], ps_microbiome_HHP) # remove empty taxa
  ps_microbiome_HHP <- prune_samples(names(which(apply(otu_table(ps_microbiome_HHP), 1, sum) != 0)), ps_microbiome_HHP) # remove patients without taxa
  # ps_microbiome_HHP.prop <- transform_sample_counts(ps_microbiome_HHP, function(otu) {if (sum(otu)==0) otu else 100*otu/sum(otu)})
  
  
# ------------------------------------------------------------------------------------------------------------------------
## Extract metadata and OTU from ps object
# ------------------------------------------------------------------------------------------------------------------------
  
## Subset Aspergillus species 
  
  #### Using Raw counts
  ps_AspergillusSpecies <- subset_taxa(ps_microbiome_HHP, Genus=="Aspergillus") 
  metadata <- data.frame(sample_data(ps_AspergillusSpecies))
  OTU <- data.frame(otu_table(ps_AspergillusSpecies))
  
  
  ## Prepare dataframe for bubble chart (long form)
  # ==============================================================================================================
  OTU$SequencingID <- row.names(OTU)
  OTU_long <- melt(OTU, id.vars = "SequencingID", variable.name = "Aspergillus_Species", 
                   value.name = "Relative_Abundance")
  OTU_long[which(OTU_long$Relative_Abundance == 0), "Relative_Abundance"] <- NA #Make 0 NA for the purpose of bubble plots. 
    
  df_metadataOTU <- merge(metadata, OTU_long, by.x = "row.names", by.y = "SequencingID", all.x = TRUE, all.y = TRUE)
  # row.names(df_metadataOTU) <- df_metadataOTU$Row.names ; df_metadataOTU$Row.names <- NULL
  
    ## Subset Air and Surfaces only  
    # --------------------------------------------------------------------------------------------------------------
    df_metadataOTU_Air_Surfaces <- subset(df_metadataOTU, SampleType != "Sputum") #Remove Sputum
    df_metadataOTU_Air_Surfaces$Aspergillus_Species <- factor(df_metadataOTU_Air_Surfaces$Aspergillus_Species,
                                                              levels = rev(sort(levels(df_metadataOTU_Air_Surfaces$Aspergillus_Species))))
    df_metadataOTU_Air_Surfaces$SampleType <- factor(df_metadataOTU_Air_Surfaces$SampleType, levels = c("Outdoor", "Indoor", "SwabS"))
    
    
# --------------------------------------------------------------------------------------------------------------
## PLOT 
# --------------------------------------------------------------------------------------------------------------

  ## Specify Facet labels
  Sample_Type_labels <- c("Indoor air", "Outdoor air", "Surfaces")
  names(Sample_Type_labels) <- c("Indoor", "Outdoor", "SwabS")
    
  ######## COLOURS #########
  Colour_palette_Aspergillus <- c("Aspergillus.aculeatinus" = "#E41A1C", 
                                  "Aspergillus.awamori" = "#377EB8",
                                  "Aspergillus.clavatus" = "#4DAF4A",
                                  "Aspergillus.fumigatus" = "#A65628",
                                  "Aspergillus.glaucus" = "#FF7F00",
                                  "Aspergillus.heteromorphus" = "#FFFF33",
                                  "Aspergillus.nidulans" = "#984EA3",
                                  "Aspergillus.niger" = "#F781BF",
                                  "Aspergillus.pseudoglaucus" = "#1B9E77",
                                  "Aspergillus.steynii" = "#D95F02",
                                  "Aspergillus.welwitschiae"="#7570B3")
  
    
  ## bubble chart --- Using Raw Read Counts
  # --------------------------------------------------------------------------------------------------------------
  p_bubblechart_Reads <- ggplot(df_metadataOTU_Air_Surfaces, 
                                aes(x = Sequencing_ID, y = Aspergillus_Species)) +
    geom_point(alpha = 0.7, shape=21, aes(fill=Aspergillus_Species, size=Relative_Abundance)) +
    scale_size_continuous(name="Metagenome reads",
                          limits = c(0, 2000),
                          range = c(1, 30)) +
    scale_fill_manual(values = Colour_palette_Aspergillus) +
    labs(y = NULL, x = NULL) + 
    facet_grid(.~SampleType, scales = "free", labeller = labeller(SampleType=Sample_Type_labels) ) + 
    guides(color = "none",
           fill = "none") + 
    theme_classic() +
    theme(axis.text.y = element_text(face = "italic"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black")) 
  p_bubblechart_Reads
  ## Export
  p_bubblechart_Reads_filename <- file.path(outdir, paste0("BubbleChart-", "AspergillusSpecies", "-", "Reads-",
                                                                  paste0(unique(df_metadataOTU_Air_Surfaces$SampleType), collapse = "_"), "-35patients", ".png"))
  ggsave(plot = p_bubblechart_Reads, filename = p_bubblechart_Reads_filename, 
         height = 12, width = 30, units = "cm")
  
  
  ## Remove legends
  # --------------------------------------------------------------------------------------------------------------
  p_bubblechart_Reads_NoLegend <- p_bubblechart_Reads + theme(legend.position = "none")
  ## Export
  p_bubblechart_Reads_NoLegend_filename <- file.path(outdir, paste0("BubbleChart-", "AspergillusSpecies", "-", "Reads-",
                                                                           paste0(unique(df_metadataOTU_Air_Surfaces$SampleType), collapse = "_"), "-35patients", "-NoLegend", ".png"))
  ggsave(plot = p_bubblechart_Reads_NoLegend, filename = p_bubblechart_Reads_NoLegend_filename, 
         height = 12, width = 28, units = "cm")
  
  
  

  
  