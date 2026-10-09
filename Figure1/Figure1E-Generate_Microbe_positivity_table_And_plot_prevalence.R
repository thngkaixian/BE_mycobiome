## Classify microbe positivity (>1% RA threshold)

library(phyloseq)
library(dplyr)
library(ggplot2)


# --------------------------------------------------------------------------------------------------------------
## Directories 
# --------------------------------------------------------------------------------------------------------------

project_name <- "MYCOBIOME_AsiaBRIDGE"
taxLevel <- "Genus"

## Directories 
wkdir <- file.path("/Users/thngkaixian/PhD_RESEARCH/RESEARCH_Projects", project_name, "SequencingData")
indir_ps <- file.path(wkdir, "FINAL_ANALYSIS")

outdir <- file.path(wkdir, "FINAL_ANALYSIS", "Analysis", "Figure_1", "02 Microbe_Positivity")
if (!dir.exists(outdir)) try(dir.create(outdir, recursive= TRUE), silent= TRUE)


# --------------------------------------------------------------------------------------------------------------
## Read ps objects 
# --------------------------------------------------------------------------------------------------------------

ps <- readRDS(file.path(indir_ps, paste0("Final_psITS_", taxLevel, ".RData")))
ps <- prune_taxa(taxa_names(ps)[!grepl('unidentified',taxa_names(ps))], ps)

  ## Filter 
  RelAbundThreshold <- 0.01
  sampleThreshold <- 0.05
  
  ps.prop <- transform_sample_counts(ps, function(otu) {if (sum(otu)==0) otu else otu/sum(otu)})
  tax <- taxa_names(ps.prop)[apply(otu_table(ps.prop), 2, function(x) {sum(x >= RelAbundThreshold) >= sampleThreshold*nsamples(ps)})]
  ps <- prune_taxa(tax, ps)
  
  ## Relative abundance
  ps.prop <- transform_sample_counts(ps, function(otu) {if (sum(otu)==0) otu else 100*otu/sum(otu)})
  

# --------------------------------------------------------------------------------------------------------------
## Determine Positive or Negative (>1% Relative abudnance) 
# --------------------------------------------------------------------------------------------------------------
  
## Extract Relative abundance
df_rel_abun <- as.data.frame(otu_table(ps.prop))
  
## Positive or Negative 
RA_threshold <- 1
df_PositiveNegative <- ifelse(df_rel_abun>RA_threshold, "Positive", "Negative")
write.csv(df_PositiveNegative, file = file.path(outdir, paste0("Microbe_positive", RA_threshold, "RA", "_table.csv")))
  
  
  ## Create a dataframe with Phylum classification
  # --------------------------------------------------------------------------------------------------------------
  Ascomycota <- c("Candida", "Saccharomyces", "Penicillium", "Clavispora", "Aspergillus", 
                  "Microidium", 
                  "Cladosporium", "Fusarium", "Mycosphaerella", "Nakaseomyces", "Chaenotheca")
  Basidiomycota <- c("Schizophyllum", "Malassezia", "Wallemia", 
                     "Cutaneotrichosporon", "Trametes", "Lentinus", "Cryptococcus", 
                     "Heterobasidion", "Rhodotorula", "Filobasidium")
  ## Create function for Phylum classification
  Phylum_classifier <- function(Genus){
    if (Genus %in% Ascomycota)"Ascomycota"
    else if (Genus %in% Basidiomycota) "Basidiomycota"
    else "Not_Ascomycota_or_Basidiomycota"
  }
  
  TaxTable <- data.frame(Genus=colnames(df_PositiveNegative))
  TaxTable$Phylum <- sapply(TaxTable$Genus, Phylum_classifier)
  

# --------------------------------------------------------------------------------------------------------------
## Determine prevalence
# --------------------------------------------------------------------------------------------------------------
  
  ## All patients 
  df_Prevalence <- data.frame(colSums(df_PositiveNegative == "Positive"))
  colnames(df_Prevalence) <- paste0("MicrobePositive")
  df_Prevalence$Microbe <- row.names(df_Prevalence) ; row.names(df_Prevalence)<-NULL
  df_Prevalence <- merge(df_Prevalence, TaxTable, by.x = "Microbe", by.y = "Genus", all.x = TRUE, all.y = FALSE)
  df_Prevalence <- relocate(df_Prevalence, Microbe)
  df_Prevalence$Prevalence <- round((df_Prevalence$MicrobePositive/nrow(df_PositiveNegative))*100, 1)
  head(df_Prevalence)
  
  ## Export 
  df_Prevalence_outfile <- file.path(outdir, paste0("Genus_prevalance", RA_threshold, "_OVERALL", ".csv"))
  write.csv(df_Prevalence, df_Prevalence_outfile)
  

# --------------------------------------------------------------------------------------------------------------
## Mycobiome Colour palette  
# --------------------------------------------------------------------------------------------------------------

microbe_colors <- c("Candida" = "#E41A1C", "Saccharomyces" = "#377EB8", "Microidium" = "#4DAF4A",
                      "Penicillium" = "#984EA3","Aspergillus" = "#FF7F00", "Malassezia" = "#FFFF33",
                      "Cladosporium" = "#A65628", "Fusarium" = "#F781BF", "Cutaneotrichosporon" = "#999999",
                      "Mycosphaerella" = "#66C2A5", "Cryptococcus" = "#FC8D62", "Schizophyllum" = "#8DA0CB",
                      "Filobasidium" = "#E78AC3", "Clavispora" = "#A6D854", "Heterobasidion" = "#FFD92F",
                      "Alternaria" = "#E5C494", "Pichia" = "#B3B3B3", "Wallemia" = "#8DD3C7",
                      "Nakaseomyces" = "#FFFFB3", "Chaenotheca" = "#BEBADA", "Lentinus" = "#FB8072",
                      "Rhodotorula" = "#80B1D3", "Trichophyton" = "#FDB462","unidentified" = "Gray10")
  
microbe_colors <- c("Candida" = "#468A4B" ,
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

  
# --------------------------------------------------------------------------------------------------------------
## Prevalance barplot 
# --------------------------------------------------------------------------------------------------------------
  
df_plot <- df_Prevalence 
  
  ## Factor microbes by prevalence 
  # --------------------------------------------------------------------------------------------------------------
  df_plot <- df_plot[order(df_plot$Prevalence, decreasing = FALSE),]
  df_plot$Microbe <- factor(df_plot$Microbe, levels = df_plot$Microbe)
  
  df_plot <- subset(df_plot, Prevalence !=0) 
  
  
  ## Horizontal barplot
  # --------------------------------------------------------------------------------------------------------------
  plot_Prevalence <- ggplot(df_plot, aes(x = Prevalence, y = Microbe, fill = Microbe)) + # 
    geom_bar(stat="identity", position=position_dodge()) +
    xlim(0,100) +
    xlab("Prevalence (%)") +
    ylab(NULL) +
    # ggtitle(paste0("RA threshold: >", RA_threshold, "%")) +
    facet_grid(Phylum~. , scales = "free", space = "free") +
    scale_fill_manual(values = microbe_colors) +
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.y = element_text(face="italic"),
          axis.line = element_line(colour = "black")) +
    geom_text(aes(label=Prevalence), hjust = 0, 
              nudge_x = 2, size = 3) ##----Increased font size (% labels)
  plot_Prevalence
  
  ## Export --- By Continent
  plot_Prevalence_filename <- file.path(outdir, paste0(taxLevel, "_prevalance", RA_threshold, "_", "OVERALL", "_Facet_Phylum", ".png"))
  ggsave(filename = plot_Prevalence_filename, plot = plot_Prevalence, width = 9, height = 15, units = "cm")
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  