library(phyloseq)
library(vegan)
library(ggplot2)
library(reshape2)


# ------------------------------------------------------------------------------------------------------------------------
# Parameters and working directory
# ------------------------------------------------------------------------------------------------------------------------

taxLevel_Microbiome <- "Genus"

wkdir <- "/Users/thngkaixian/PhD_RESEARCH/RESEARCH_Projects"
indir_ps <- file.path(wkdir, "MYCOBIOME_AsiaBRIDGE", "SequencingData", "FINAL_ANALYSIS")

## Out directories
outdir <- file.path(indir_ps, "Analysis", "Figure_5", "PLOTTING-Fungal-Genus-RelAbund-across-SampleTypes", taxLevel_Microbiome)
if (!dir.exists(outdir)) try(dir.create(outdir, recursive= TRUE), silent= TRUE)


# --------------------------------------------------------------------------------------------------------------
## Import data: HHP data
# --------------------------------------------------------------------------------------------------------------  

ps_microbiome_HHP <- readRDS(file.path(indir_ps, "Final_ps_Genus_HHP_35patients.RData"))

  ps_microbiome_HHP <- prune_taxa(taxa_names(ps_microbiome_HHP)[apply(otu_table(ps_microbiome_HHP), 2, sum) != 0], ps_microbiome_HHP) # remove empty taxa
  ps_microbiome_HHP <- prune_samples(names(which(apply(otu_table(ps_microbiome_HHP), 1, sum) != 0)), ps_microbiome_HHP) # remove patients without taxa
  ps_microbiome_HHP.prop <- transform_sample_counts(ps_microbiome_HHP, function(otu) {if (sum(otu)==0) otu else 100*otu/sum(otu)})
  
  
# ------------------------------------------------------------------------------------------------------------------------
## Extract metadata and OTU from ps object
# ------------------------------------------------------------------------------------------------------------------------
  
## Subset Fungi
  
  #### Using Raw counts
  ps_Eukaryota <- subset_taxa(ps_microbiome_HHP, Superkingdom=="Eukaryota")
  metadata <- data.frame(sample_data(ps_Eukaryota))
  metadata <- metadata[, -c(53:304)] #Delete allergen data here
  OTU <- data.frame(otu_table(ps_Eukaryota))

  
  ## Prepare dataframe for bubble chart (long form)
  # ==============================================================================================================
  OTU$SequencingID <- row.names(OTU)
  OTU_long <- melt(OTU, id.vars = "SequencingID", variable.name = "Genus", 
                   value.name = "Relative_Abundance")
  OTU_long[which(OTU_long$Relative_Abundance == 0), "Relative_Abundance"] <- NA #Make 0 NA for the purpose of bubble plots. 
    
  df_metadataOTU <- merge(metadata, OTU_long, by.x = "row.names", by.y = "SequencingID", all.x = TRUE, all.y = TRUE)
  row.names(df_metadataOTU) <- df_metadataOTU$Row.names ; df_metadataOTU$Row.names <- NULL
  
    ## Subset Air and Surfaces only (Remove Sputum and Inhalers)
    # --------------------------------------------------------------------------------------------------------------
    df_metadataOTU_Air_Surfaces <- subset(df_metadataOTU, SampleType != "Sputum") #Remove Sputum
    df_metadataOTU_Air_Surfaces$Genus <- droplevels(df_metadataOTU_Air_Surfaces$Genus)
    df_metadataOTU_Air_Surfaces$Genus <- factor(df_metadataOTU_Air_Surfaces$Genus,
                                                levels = rev(sort(levels(df_metadataOTU_Air_Surfaces$Genus))))
    df_metadataOTU_Air_Surfaces$SampleType <- factor(df_metadataOTU_Air_Surfaces$SampleType, levels = c("Outdoor", "Indoor", "SwabS"))
    dim(df_metadataOTU_Air_Surfaces)
    
  
# --------------------------------------------------------------------------------------------------------------
## PLOT 
# --------------------------------------------------------------------------------------------------------------

  ## Specify Facet labels
  Sample_Type_labels <- c("Indoor air", "Outdoor air", "Surfaces")
  names(Sample_Type_labels) <- c("Indoor", "Outdoor", "SwabS")
    
  ## Colours
  Colour_palette <- c("indianred1","dodgerblue1","thistle","lightskyblue", "darkolivegreen1","cyan3","bisque", "royalblue", "khaki2", "pink",
             "yellow", "limegreen", "coral2", "violetred", "indianred2", "darkblue",
             "darkorange1", "cyan1", "royalblue4", "maroon2", "darkblue", "aquamarine",
             "deepskyblue2", "slateblue", "dodgerblue3", "darkolivegreen1", "darkgoldenrod1",
             "violetred", "grey","slategrey","black", "pink3",
             "bisque", "lightblue", "darkblue", "cadetblue", "indianred1", "turquoise2",
             "cyan1", "cyan3", "thistle", "salmon2", "blue2" )
  
    
  ## bubble chart --- Using Raw Read Counts
  # --------------------------------------------------------------------------------------------------------------
  min(df_metadataOTU_Air_Surfaces$Relative_Abundance, na.rm = TRUE) #bottom end of scale
  max(df_metadataOTU_Air_Surfaces$Relative_Abundance, na.rm = TRUE) #upper end of scale
  length(sort(unique(df_metadataOTU_Air_Surfaces$Genus))) #30 Genus
  
  
  p_bubblechart_Reads <- ggplot(df_metadataOTU_Air_Surfaces, 
                          aes(x = Sequencing_ID, y = Genus)) + #, size = Relative_Abundance, color = Genus
    geom_point(alpha = 0.7, shape=21, aes(fill=Genus, size=Relative_Abundance)) +
    scale_size_continuous(name="Metagenome reads", # Rename legend title here
                          limits = c(0, 20000), ##----------Ivan's comments about the circles.... 
                          range = c(1, 30)) +
      scale_fill_manual(values = Colour_palette) +
    labs(y = NULL, x = NULL) + 
    facet_grid(.~SampleType, scales = "free", labeller = labeller(SampleType=Sample_Type_labels) ) + #space = "free", | "Home.Sampling.Type" --- Specific Sample Type for surfaces  
    guides(color = "none",
           fill = "none") + 
    theme_classic() +
    theme(#legend.key = element_rect(fill = "white"), #-----Makes the legend background white
          axis.text.y = element_text(face = "italic"),
          # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black")) 
  p_bubblechart_Reads
  ## Export
  p_bubblechart_Reads_filename <- file.path(outdir, paste0("BubbleChart-", "Fungi", "-", "Reads-",
                                                            paste0(unique(df_metadataOTU_Air_Surfaces$SampleType), collapse = "_"), "-35patients", ".png"))
  ggsave(plot = p_bubblechart_Reads, filename = p_bubblechart_Reads_filename, 
         height = 18, width = 30, units = "cm")

  
    ## Remove legends
    # --------------------------------------------------------------------------------------------------------------
    p_bubblechart_Reads_NoLegend <- p_bubblechart_Reads + theme(legend.position = "none")
    ## Export
    p_bubblechart_Reads_NoLegend_filename <- file.path(outdir, paste0("BubbleChart-", "Fungi", "-", "Reads-",
                                                                    paste0(unique(df_metadataOTU_Air_Surfaces$SampleType), collapse = "_"), "-35patients", "-NoLegend", ".png"))
    ggsave(plot = p_bubblechart_Reads_NoLegend, filename = p_bubblechart_Reads_NoLegend_filename, 
           height = 18, width = 28, units = "cm")
    
  

  
  