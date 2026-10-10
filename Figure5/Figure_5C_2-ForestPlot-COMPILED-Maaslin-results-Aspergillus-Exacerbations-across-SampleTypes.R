## Barplots of Maaslin results 

# Load libraries 
library(ggplot2)


# ------------------------------------------------------------------------------------------------------------------------
# Parameters and working directory
# ------------------------------------------------------------------------------------------------------------------------

project_name <- "HomeHospitals"
taxLevel <- "Species"

wkdir <- "/Users/thngkaixian/PhD_RESEARCH/RESEARCH_Projects/MYCOBIOME_AsiaBRIDGE/SequencingData/FINAL_ANALYSIS"
indir <- file.path(wkdir, "Analysis/Figure_5C/Maaslin-Aspergillus-Species-vs-ClinicalOutcome", taxLevel)

  ## Maaslin results directories
  # --------------------------------------------------------------------------------------------------------------

  #### Microbiome subset
  MicrobiomeSubset <- "AspergillusSpecies"
  Maaslin_Model <- "ModelCPLM_NormalisationNONE_TransformationNONE"
  Results_Folder <- paste0("COMPILED_across_SampleTypes", "_", MicrobiomeSubset)
  
  ## Create out directory
  outdir <- file.path(indir, Results_Folder, paste0(Maaslin_Model, "_RESULTS"))
  if (!dir.exists(outdir)) try(dir.create(outdir, recursive= TRUE), silent= TRUE)
  
  
# --------------------------------------------------------------------------------------------------------------
## Import data
# --------------------------------------------------------------------------------------------------------------
  
Clinical_Var <- "Exacerbations"

qval_cutoff <- 0.1

  
#### Outdoor ####
df_Results.Outdoor <- read.delim(file.path(indir, "Outdoor", paste0("Outdoor", "_", MicrobiomeSubset, "_", Clinical_Var), Maaslin_Model, "all_results.tsv"),  header = TRUE, sep = "\t")
df_Results.Outdoor <- subset(df_Results.Outdoor, df_Results.Outdoor$qval<qval_cutoff)
  
#### Indoor ####
df_Results.Indoor <- read.delim(file.path(indir, "Indoor", 
                                          paste0("Indoor", "_", MicrobiomeSubset, "_", Clinical_Var), 
                                          Maaslin_Model, "all_results.tsv"),  header = TRUE, sep = "\t")
df_Results.Indoor <- subset(df_Results.Indoor, df_Results.Indoor$qval<qval_cutoff)

#### SwabS ####
df_Results.SwabS <- read.delim(file.path(indir, "SwabS",
                                           paste0("SwabS", "_", MicrobiomeSubset, "_", Clinical_Var ),
                                           Maaslin_Model, "all_results.tsv"),  header = TRUE, sep = "\t")
  df_Results.SwabS <- subset(df_Results.SwabS, df_Results.SwabS$qval<qval_cutoff)
  

# --------------------------------------------------------------------------------------------------------------
## Combine data into single dataframe for plotting 
# --------------------------------------------------------------------------------------------------------------

  ## Create column - Indicate Sample Type 
  df_Results.Outdoor$SampleType <- "Outdoor"
  df_Results.Indoor$SampleType <- "Indoor"
  df_Results.SwabS$SampleType <- "SwabS"
  
  
  ## Combine data 
  # --------------------------------------------------------------------------------------------------------------
  df_Results_COMPILED <- rbind(df_Results.Outdoor, df_Results.Indoor, df_Results.SwabS)
  df_Results_COMPILED$SampleType <- factor(df_Results_COMPILED$SampleType, levels = c("Outdoor", "Indoor", "SwabS"), labels = c("Outdoor", "Indoor", "Surfaces"))  
  df_Results_COMPILED$feature <- factor(df_Results_COMPILED$feature, levels = rev(levels(factor(df_Results_COMPILED$feature))))

  
  ## Indicate Positive/Negative correlation 
  # --------------------------------------------------------------------------------------------------------------
  df_Results_COMPILED$coef_PositiveNegative <- ifelse(df_Results_COMPILED$coef > 0, "Positive", "Negative")
  
  
  ## Calculate 95% Confidence Intervals
  # --------------------------------------------------------------------------------------------------------------
  df_Results_COMPILED$ci_lower <- df_Results_COMPILED$coef - 1.96 * df_Results_COMPILED$stderr
  df_Results_COMPILED$ci_upper <- df_Results_COMPILED$coef + 1.96 * df_Results_COMPILED$stderr
  
  ## Export 
  write.csv(df_Results_COMPILED, file = file.path(outdir, paste0("COMPILED_Results", "_qval", qval_cutoff, "_",
                                                                 "Indoor-Outdoor-Surfaces", "-", Clinical_Var, ".csv")))
  

# --------------------------------------------------------------------------------------------------------------
## plot 
# --------------------------------------------------------------------------------------------------------------
  
  
  ## Forest Plot
  # --------------------------------------------------------------------------------------------------------------
  
  p_Forest <- ggplot(df_Results_COMPILED, aes(x = coef, y = feature)) +
    geom_point(shape = 16) +  
    geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), width = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    facet_grid(.~SampleType) + 
    xlab("Exacerbation frequency\n(Coefficient (95% CI))") +
    ylab(NULL) +
    xlim(-2,2) +
    theme(panel.border = element_blank(),
          panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.text.y = element_text(face = "italic")
    )
  p_Forest
  plot_filename <- file.path(outdir, paste0("ForestPlot", "_", Clinical_Var, "_qval", qval_cutoff, ".png"))
  ggsave(plot_filename, plot = p_Forest, width = 18, height = 8, units = "cm") 

  
  
  
  