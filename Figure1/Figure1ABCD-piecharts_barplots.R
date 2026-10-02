library(phyloseq)
library(microViz)
library(microbiome)
library(RColorBrewer)


# --------------------------------------------------------------------------------------------------------------
## Parameters and directories 
# --------------------------------------------------------------------------------------------------------------


project_name <- "MYCOBIOME_AsiaBRIDGE"

taxLevel <- "Phylum"
taxLevel <- "Genus"

## Directories 
wkdir <- file.path("/Users/thngkaixian/PhD_RESEARCH/RESEARCH_Projects", project_name, "SequencingData")
indir_ps <- file.path(wkdir, "FINAL_ANALYSIS")

outdir <- file.path(wkdir, "FINAL_ANALYSIS", "Analysis", "Figure_1")
if (!dir.exists(outdir)) try(dir.create(outdir, recursive= TRUE), silent= TRUE)


# --------------------------------------------------------------------------------------------------------------
## Read ps object
# --------------------------------------------------------------------------------------------------------------

ps.prop <- readRDS(file.path(indir_ps, paste0("Final_psITS_", taxLevel, "_RelativeAbundance.RData")))

  # Top taxa
  # ------------------------------------------------------------------------------------------
  top25 <- names(sort(taxa_sums(ps.prop), decreasing=TRUE))[1:min(25,ntaxa(ps.prop))] ; top25
  ps.prop.topN <- prune_taxa(top25, ps.prop)
  ps.prop.topN <- tax_reorder(ps.prop.topN, top25)
  
  ps.final <- ps.prop.topN 
  
  
# --------------------------------------------------------------------------------------------------------------
# Color palette
# --------------------------------------------------------------------------------------------------------------
  
  ## Phylum
  if(taxLevel=="Phylum"){
    colTaxa <- data.frame(Microbe= c("Ascomycota", "Basidiomycota"), #"Basidiomycota", "Ascomycota"
                          Hex_code=c("#F8D458", "#204E86")) # "#204E86", "#F8D458"
    row.names(colTaxa) <- colTaxa$Microbe ; colTaxa$Microbe <- NULL
  }
  
  ## Genus
  if(taxLevel=="Genus"){
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
  }
  

# --------------------------------------------------------------------------------------------------------------
## ASSIGN ANALYSIS OF INTEREST
# --------------------------------------------------------------------------------------------------------------

  PS_plot <- ps.final ; COHORT<-"AllPatients" ; COHORT


# --------------------------------------------------------------------------------------------------------------
# Piechart
# --------------------------------------------------------------------------------------------------------------
  
# GROUP_VAR <- "Disease"
  
  ## Dataframe for piechart
  agg <- data.frame(
    Geography= rep(COHORT, ntaxa(PS_plot)),
    Taxonomy= taxa_names(PS_plot),
    Relative_Abundance= as.numeric(apply(otu_table(PS_plot),2,mean))
  )
  agg$Taxonomy <- factor(agg$Taxonomy, levels= taxa_names(PS_plot))
  
  ## Plot pie chart 
  pie_chart <- ggplot(agg, aes(x="", y=Relative_Abundance, fill=Taxonomy)) +
    geom_col() +
    coord_polar(theta = "y") +
    guides(fill= guide_legend(title=taxLevel,ncol= 1)) +
    scale_fill_manual(values=as.character(colTaxa[taxa_names(PS_plot),])) +
    theme_bw() +
    theme(legend.text = element_text(face="italic"),
      legend.key.height = unit(1, "mm"),
      legend.background = element_rect(fill='transparent'),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.title.y=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank())
  pie_chart 
  
  # # Output barplot into .png file
  # output_pie_chart <- file.path(outdir_filter, paste0("PieChart_", taxLevel, '-Aggregate', "-", COHORT, "-", GROUP_VAR, '.png'))
  # ggsave(output_pie_chart, plot = pie_chart,  width = 12, height = 8, units = "cm",  dpi = 300)
  
  
  ## Remove legend
  # ----------------------------------------------------------------------
  pie_chart_noLegend <- pie_chart +
    theme(legend.position = "none",
          panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA) #transparent plot bg
          ) +
    ggtitle(NULL)
  pie_chart_noLegend
  
  # Output barplot into .png file
  output_pie_chart <- file.path(outdir, paste0("PieChart_", taxLevel, '-Aggregate', "-", COHORT, "-", GROUP_VAR,  "-removedLegendAndTitle", '.png'))
  ggsave(output_pie_chart, plot = pie_chart_noLegend,  bg="transparent", width = 4, height = 4, units = "cm",  dpi = 300)
   
  
# --------------------------------------------------------------------------------------------------------------
# Individual barplots - Facet Sample Source
# --------------------------------------------------------------------------------------------------------------
  
  
  ## Plot --- Phylum
  # --------------------------------------------------------------------------------------------------------------
  if(taxLevel=="Phylum"){
    sample_names_sorted <- sample_names(PS_plot)[order(otu_table(PS_plot)[, "Ascomycota"], decreasing = TRUE)]  # Get sorted sample names
    otu_table(PS_plot) <- otu_table(PS_plot)[, c("Ascomycota", "Basidiomycota")] #"Basidiomycota", "Ascomycota"
  }

  p_Indiv <- plot_bar(ps_reorder(PS_plot, sample_names_sorted), x="Sample_Name", fill="Phylum") +
    geom_bar(stat="identity") +
    guides(fill= guide_legend(ncol= 1)) +
    xlab("Bronchiectasis patients") +
    scale_x_discrete(limits = sample_names_sorted) +  # Manually set order
    ylab('Relative abundance (%)') +
    scale_y_continuous(expand = c(0, 0)) +  # Remove space between bars and X-axis
    scale_fill_manual(values= as.character(colTaxa[taxa_names(PS_plot), ]) ) +
    theme(legend.position = "right",
          legend.text = element_text(face = "italic"), #, size = 7
          legend.key.height = unit(1, "mm"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(color="black"),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
  p_Indiv$data$Phylum <- factor(p_Indiv$data$Phylum, levels= taxa_names(PS_plot))
  p_Indiv
  
  ## Output barplot into .png file
  output_barplot <- file.path(outdir, paste0("Barplot_", taxLevel, '-Individual', "-", COHORT, "-", GROUP_VAR, '.png'))
  ggsave(output_barplot, plot = p_Indiv, height = 10, width = 20, units = "cm") 
  
  

  
  ## Plot --- Genus
  # --------------------------------------------------------------------------------------------------------------
  
  ## Sort by Candida abundance 
  if(taxLevel=="Genus") sample_names_sorted <- sample_names(PS_plot)[order(otu_table(PS_plot)[, "Candida"], decreasing = TRUE)]  # Get sorted sample names
  
  p_Indiv <- plot_bar(ps_reorder(PS_plot, sample_names_sorted), x="Sample_Name", fill="Genus") + 
    geom_bar(stat="identity") +
    guides(fill= guide_legend(ncol= 1)) +
    xlab("Bronchiectasis patients") +
    scale_x_discrete(limits = sample_names_sorted) +  # Manually set order
    ylab('Relative abundance (%)') +
    scale_y_continuous(expand = c(0, 0)) +  # Remove space between bars and X-axis
    scale_fill_manual(values= as.character(colTaxa[taxa_names(PS_plot), ]) ) +
    theme(legend.position = "right",
          legend.text = element_text(face = "italic"), #, size = 7
          legend.key.height = unit(1, "mm"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(color="black"),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
  p_Indiv$data$Genus <- factor(p_Indiv$data$Genus, levels= taxa_names(PS_plot))
  p_Indiv
      
  ## Output barplot into .png file
  output_barplot <- file.path(outdir, paste0("Barplot_", taxLevel, '-Individual', "-", COHORT, "-", GROUP_VAR, '.png'))
  ggsave(output_barplot, plot = p_Indiv, height = 10, width = 20, units = "cm")  
      

  