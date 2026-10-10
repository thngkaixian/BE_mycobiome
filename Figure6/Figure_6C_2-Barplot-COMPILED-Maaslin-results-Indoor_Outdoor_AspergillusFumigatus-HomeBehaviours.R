## Barplots of Maaslin results 

# Load libraries 
library(ggplot2)


# ------------------------------------------------------------------------------------------------------------------------
# Parameters and working directory
# ------------------------------------------------------------------------------------------------------------------------

## Parameters
project_name <- "HomeHospitals"
taxLevel <- "Species"


## Directories 
wkdir <- "/Users/thngkaixian/PhD_RESEARCH/RESEARCH_Projects/MYCOBIOME_AsiaBRIDGE/SequencingData/FINAL_ANALYSIS"
indir <- file.path(wkdir, "Analysis", "Figure_6CD", "Maaslin-Aspergillus-Species-vs-HomeBehaviour", taxLevel)

  ## Maaslin results directories
  # --------------------------------------------------------------------------------------------------------------

  #### Microbiome subset
  MicrobiomeSubset <- "AspergillusSpecies"
  Maaslin_Model <- "ModelCPLM_NormalisationNONE_TransformationNONE"
  Results_OutFolder <- paste0("COMPILED_across_SampleTypes", "_", MicrobiomeSubset)
  
  ## Create out directory
  outdir <- file.path(indir, Results_OutFolder, paste0(Maaslin_Model, "_RESULTS"))
  if (!dir.exists(outdir)) try(dir.create(outdir, recursive= TRUE), silent= TRUE)
  
 
# --------------------------------------------------------------------------------------------------------------
## Import data
# --------------------------------------------------------------------------------------------------------------
  
qval_cutoff <- 0.1

  #### Outdoor ####
  # --------------------------------------------------------------------------------------------------------------
  df_Aircon_1234.Outdoor <- read.delim(file.path(indir, "Outdoor", paste0("Outdoor", "_", MicrobiomeSubset, "_", "Aircon_use_room"), Maaslin_Model, "all_results.tsv"),  header = TRUE, sep = "\t")
  df_Aircon_1234.Outdoor$SampleType <- "Outdoor"
  
  df_Fan.Outdoor <- read.delim(file.path(indir, "Outdoor", paste0("Outdoor", "_", MicrobiomeSubset, "_", "Use_of_fan_in_room"), Maaslin_Model, "all_results.tsv"),  header = TRUE, sep = "\t")
  df_Fan.Outdoor$SampleType <- "Outdoor"
  
  df_Window_Never_Yes.Outdoor <- read.delim(file.path(indir, "Outdoor", paste0("Outdoor", "_", MicrobiomeSubset, "_", "Window_open_room_NeverYes"), Maaslin_Model, "all_results.tsv"),  header = TRUE, sep = "\t")
  df_Window_Never_Yes.Outdoor$SampleType <- "Outdoor"
  
  df_Window_12h.Outdoor <- read.delim(file.path(indir, "Outdoor", paste0("Outdoor", "_", MicrobiomeSubset, "_", "Window_open_room_hrs_12h"), Maaslin_Model, "all_results.tsv"),  header = TRUE, sep = "\t")
  df_Window_12h.Outdoor$SampleType <- "Outdoor"
  

  #### Indoor ####
  # --------------------------------------------------------------------------------------------------------------
  df_Aircon_1234.Indoor <- read.delim(file.path(indir, "Indoor", paste0("Indoor", "_", MicrobiomeSubset, "_", "Aircon_use_room"), Maaslin_Model, "all_results.tsv"),  header = TRUE, sep = "\t")
  df_Aircon_1234.Indoor$SampleType <- "Indoor"
  
  df_Fan.Indoor <- read.delim(file.path(indir, "Indoor", paste0("Indoor", "_", MicrobiomeSubset, "_", "Use_of_fan_in_room"), Maaslin_Model, "all_results.tsv"),  header = TRUE, sep = "\t")
  df_Fan.Indoor$SampleType <- "Indoor"

  df_Window_Never_Yes.Indoor <- read.delim(file.path(indir, "Indoor", paste0("Indoor", "_", MicrobiomeSubset, "_", "Window_open_room_NeverYes"), Maaslin_Model, "all_results.tsv"),  header = TRUE, sep = "\t")
  df_Window_Never_Yes.Indoor$SampleType <- "Indoor"
  
  df_Window_12h.Indoor <- read.delim(file.path(indir, "Indoor", paste0("Indoor", "_", MicrobiomeSubset, "_", "Window_open_room_hrs_12h"), Maaslin_Model, "all_results.tsv"),  header = TRUE, sep = "\t")
  df_Window_12h.Indoor$SampleType <- "Indoor"
  
  
# --------------------------------------------------------------------------------------------------------------
## Combine data into single dataframe for plotting 
# --------------------------------------------------------------------------------------------------------------

  ## Combine data 
  # --------------------------------------------------------------------------------------------------------------
  df_Results_COMPILED <- rbind(df_Aircon_1234.Outdoor, df_Fan.Outdoor, df_Window_Never_Yes.Outdoor, df_Window_12h.Outdoor,
                               df_Aircon_1234.Indoor, df_Fan.Indoor, df_Window_Never_Yes.Indoor, df_Window_12h.Indoor)
  
  
  ## Factor 
  df_Results_COMPILED$metadata <- factor(df_Results_COMPILED$metadata, 
                                         levels = c("Aircon_use_room_1234", "Use_of_fan_in_room", "Window_open_room_Daily_012days", "Window_open_room_NeverYes", "Window_open_room_hrs_12h"),
                                         labels = c("Aircon use (hours)", "Use of fan", "Window use\n(Daily)", "Window use\n(At least 1 day per week)", "Window use\n(>12 hours per day)"))
  
  df_Results_COMPILED$SampleType <- factor(df_Results_COMPILED$SampleType, levels = c("Indoor", "Outdoor"))
  
  
  ## Indicate Positive/Negative correlation 
  # --------------------------------------------------------------------------------------------------------------
  df_Results_COMPILED$coef_PositiveNegative <- ifelse(df_Results_COMPILED$coef > 0, "Positive", "Negative")
  
  
  ## Calculate 95% Confidence Intervals
  # --------------------------------------------------------------------------------------------------------------
  df_Results_COMPILED$ci_lower <- df_Results_COMPILED$coef - 1.96 * df_Results_COMPILED$stderr
  df_Results_COMPILED$ci_upper <- df_Results_COMPILED$coef + 1.96 * df_Results_COMPILED$stderr
  
  
  ## Export 
  write.csv(df_Results_COMPILED, file = file.path(outdir, paste0("COMPILED_Results", "_qval", qval_cutoff, "_",
                                                                 "Indoor-Outdoor-AspFumigatus", "-V2-Window_NeverYes", ".csv")))
  
  
# --------------------------------------------------------------------------------------------------------------
## plot 
# --------------------------------------------------------------------------------------------------------------

df_Results_AspFumigatus <- subset(df_Results_COMPILED, feature=="Aspergillus.fumigatus")

  
  ## Barplot 
  p_bar<-ggplot(df_Results_AspFumigatus, aes(x = coef, y = feature, fill = SampleType)) +
    geom_bar(stat="identity", position=position_dodge()) + #position_dodge(width = 1)
    facet_grid(.~metadata) +
    scale_x_continuous(name = "Coefficient", breaks = seq(0, 2, by=0.5), limits = c(0, 2)) +
    labs(#x = "Coefficient (95% CI)",
         y = NULL,
         fill = NULL)  + 
    scale_fill_manual(values = c("Indoor" = "magenta4", "Outdoor" = "orange"), #Orange | purple
                      guide = guide_legend(reverse = TRUE)) +
    theme(legend.key.height = unit(4, "mm"), # Adjust key height
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text.y = element_text(face = "italic"),
          axis.line = element_line(colour = "black"))
  p_bar
  
  # plot_filename <- file.path(outdir, paste0("barplot", "_qval", qval_cutoff, "_",
  #                                           "Indoor-Outdoor-AspFumigatus", "-V1-Window_DailyUse", ".png")) #"-", paste0(unique(df_Results_COMPILED$metadata), collapse = "-"),
  plot_filename <- file.path(outdir, paste0("barplot", "_qval", qval_cutoff, "_",
                                            "Indoor-Outdoor-AspFumigatus", "-V2-Window_NeverYes", ".png")) #"-", paste0(unique(df_Results_COMPILED$metadata), collapse = "-"),
  
  ggsave(plot_filename, plot = p_bar, width = 25, height = 4, units = "cm") 
  

  
  
  
  
  
  

  
  
  
  