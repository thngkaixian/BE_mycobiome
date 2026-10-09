library(ggplot2)
library(ggtext)
library(dplyr)


# --------------------------------------------------------------------------------------------------------------
## Directories 
# --------------------------------------------------------------------------------------------------------------

Project_name <- "MYCOBIOME_AsiaBRIDGE"
taxLevel <- "Genus"

## Directories 
wkdir <- file.path("/Users/thngkaixian/PhD_RESEARCH/RESEARCH_Projects", Project_name, "SequencingData")
indir <- file.path(wkdir, "FINAL_ANALYSIS", "Analysis", "Figure_3A")

Maaslin_Model <- "ModelCPLM_NormalisationNONE_TransformationNONE"
outdir <- file.path(indir, "00_COMPILED_MAASLIN", "PLOTS_BubbleChart_Heatmap", paste0(Maaslin_Model, "_RESULTS"))
## Create out directory
if (!dir.exists(outdir)) try(dir.create(outdir, recursive= TRUE), silent= TRUE)


# --------------------------------------------------------------------------------------------------------------
## Import data
# --------------------------------------------------------------------------------------------------------------

FrequentExacerbator <- read.delim(file.path(indir, paste0("FrequentExacerbator", "_All_Patients"), Maaslin_Model, "all_results.tsv"), header = TRUE, sep = "\t")
BSI <- read.delim(file.path(indir, paste0("BSI", "_All_Patients"), Maaslin_Model, "all_results.tsv"), header = TRUE, sep = "\t")
MRC <- read.delim(file.path(indir, paste0("MRC", "_All_Patients"), Maaslin_Model, "all_results.tsv"), header = TRUE, sep = "\t")
FEV1 <- read.delim(file.path(indir, paste0("FEV1", "_All_Patients"), Maaslin_Model, "all_results.tsv"), header = TRUE, sep = "\t")
HospitalisationYesNo <- read.delim(file.path(indir, paste0("HospitalisationYesNo", "_All_Patients"), Maaslin_Model, "all_results.tsv"), header = TRUE, sep = "\t")


# --------------------------------------------------------------------------------------------------------------
## IMPORTANT !!!  Ensure there are 21 fungi --- for purpose of heatmap
# --------------------------------------------------------------------------------------------------------------
  
Ascomycota <- c("Candida", "Saccharomyces", "Penicillium", "Clavispora", "Aspergillus", 
                "Microidium", 
                "Cladosporium", "Fusarium", "Mycosphaerella", "Nakaseomyces", "Chaenotheca")
Basidiomycota <- c("Schizophyllum", "Malassezia", "Wallemia", 
                   "Cutaneotrichosporon", "Trametes", "Lentinus", "Cryptococcus", 
                   "Heterobasidion", "Rhodotorula", "Filobasidium")
AllFungi <- c(Ascomycota, Basidiomycota)


  ## Append with "excluded" / "missing" fungi
  # --------------------------------------------------------------------------------------------------------------

  ##-----BSI
  df_BSI_ExcludedFungi <- data.frame(feature = setdiff(AllFungi, BSI$feature), metadata = "BSI", value = "BSI", coef = NA , stderr = NA, N = NA, N.not.0 = NA, pval = NA, qval = NA)
  BSI <- rbind(BSI, df_BSI_ExcludedFungi)

  ##-----MRC
  df_MRC_ExcludedFungi <- data.frame(feature = setdiff(AllFungi, MRC$feature), metadata = "MRC_score", value = "MRC_score", coef = NA , stderr = NA, N = NA, N.not.0 = NA, pval = NA, qval = NA)
  MRC <- rbind(MRC, df_MRC_ExcludedFungi)

  ##-----FEV1
  df_FEV1_ExcludedFungi <- data.frame(feature = setdiff(AllFungi, FEV1$feature), metadata = "FEV1_percent_predicted", value = "FEV1_percent_predicted", coef = NA , stderr = NA, N = NA, N.not.0 = NA, pval = NA, qval = NA)
  FEV1 <- rbind(FEV1, df_FEV1_ExcludedFungi)
  
  
# --------------------------------------------------------------------------------------------------------------
## Combine data for heatmap 
# --------------------------------------------------------------------------------------------------------------

## Merge into single dataframe   
df_Combined <- rbind(FrequentExacerbator, BSI, MRC, FEV1, HospitalisationYesNo)


  ## Replace non-significant qval with NA
  # --------------------------------------------------------------------------------------------------------------
  qval_cutoff <- 0.1
  df_Combined$qval_Sigf <- sapply(df_Combined$qval, function(x){ifelse(x>=qval_cutoff | is.na(x), NA, x)})
  
  ## Replace coeff with 0 if qval non-significant
  # --------------------------------------------------------------------------------------------------------------
  df_Combined[,"coef_Sigf"] <- df_Combined$coef
  df_Combined[which(df_Combined$qval >= qval_cutoff) , "coef_Sigf"] <- 0
  sort(df_Combined$coef)
  # sort(df_Combined$coef_Sigf)
  
  
  ## Add Ascomycota/Basidiomycota classifier 
  # --------------------------------------------------------------------------------------------------------------
  df_Combined$Phylum <- sapply(df_Combined$feature, function(x){ if(x %in% Ascomycota ) "Ascomycota" else if (x %in% Basidiomycota) "Basidiomycota" })
  df_Combined$Phylum <- sapply(df_Combined$feature, function(x){ if(x %in% Ascomycota ) "Ascomycota" else if (x %in% Basidiomycota) "Basidiomycota" })
  

  ##----Export dataframe for viewing
  # df_Combined$qval_asterisks <- sapply(df_Combined$qval_Sigf, function(x){ if(x<0.001 & !is.na(x)) "***" else if(x<0.01 & !is.na(x)) "**" else if (x<0.1 & !is.na(x)) "*" else if(is.na(x)) "NS" })
  # df_Combined_filename <- file.path(outdir, paste0("COMPILED_MaaslinResults_", "ClinicalOutcome", "_in_AllPatients", "-", paste0(unique(df_Combined$metadata), collapse = "_"), ".csv"))
  # write.csv(df_Combined, df_Combined_filename)
  df_Combined$metadata
  

#--------------------------------------------------------------------------------------------------------------
# PLOT 
# --------------------------------------------------------------------------------------------------------------
  
  ## Factor Clinical Variables
  # --------------------------------------------------------------------------------------------------------------
  df_Combined$metadata <- factor(df_Combined$metadata,
                                 levels = rev(c("FrequentExacerbator", "Hospitalised_for_bronchiectasis_YesNo", "BSI" , "MRC_score", "FEV1_percent_predicted")),
                                 labels = rev(c("Frequent Exacerbator", "Severe Exacerbation", "BSI" , "MRC score", "FEV1 (% predicted)")))
  
  
  ## Heatmap
  # --------------------------------------------------------------------------------------------------------------
  p_heatmap <- ggplot(df_Combined, aes(x = feature, y = metadata, fill = coef)) + #fill = coef | log_qval_coef | coef_Sigf  
    geom_tile(color = "black") +
    coord_fixed() +
    xlab(NULL) + ylab(NULL) +
    labs(fill="Coefficient") +
    scale_fill_gradient2(low = "red", mid = "white", high = "green4", 
                         midpoint = 0,
                         limits = c(-1, 1),
                         na.value = "white",
                         oob = scales::squish
                         ) + 
    theme(legend.position = "top", legend.justification = "right",
          axis.text.x = element_markdown(angle = 45, hjust = 1, vjust = 1),
          panel.background = element_blank()) #+
  p_heatmap
  
  plot_filename <- file.path(outdir, paste0("Heatmap_", "ClinicalOutcome", "_in_AllPatients", "-", paste0(unique(df_Combined$metadata), collapse = "_"), ".png"))
  ggsave(plot_filename, plot = p_heatmap, width = 24, height = 10, units = "cm") # 5 rows
  
      #### Without Coefficient
      p_heatmap_NoLegend <- p_heatmap + theme(legend.position = "none")
      plot_filename <- file.path(outdir, paste0("Heatmap_", "ClinicalOutcome", "_in_AllPatients", "-", paste0(unique(df_Combined$metadata), collapse = "_"), "-NoLegend", ".png"))
      ggsave(plot_filename, plot = p_heatmap_NoLegend, width = 24, height = 8, units = "cm") # 5 rows

  
      
      
  ## Heatmp - with BOLD text for specific fungi 
  # --------------------------------------------------------------------------------------------------------------
      
    ## Specify fungi order
    df_Combined$feature <- factor(df_Combined$feature, levels = sort(unique(df_Combined$feature)))
    feature_order <- sort(unique(df_Combined$feature))
      
    ## Specify which fungi to bold       
    df_Combined <- df_Combined %>%
      mutate(feature_label = case_when(feature %in% c("Aspergillus", "Chaenotheca", "Heterobasidion", "Nakaseomyces", "Saccharomyces") ~ paste0("<b><i>", feature, "</i></b>"), TRUE ~ paste0("<i>", feature, "</i>")))
      
    ## create label map
    label_map <- df_Combined %>%
      distinct(feature) %>%
      mutate(feature_label = ifelse(
        feature %in% c("Aspergillus", "Chaenotheca", "Heterobasidion", "Nakaseomyces", "Saccharomyces"),
        paste0("<b><i>", feature, "</i></b>"),
        paste0("<i>", feature, "</i>")
      )) %>%
      tibble::deframe()
      
    
    p_heatmap_BOLD <- ggplot(df_Combined, aes(x = feature, y = metadata, fill = coef)) +
      geom_tile(color = "black") +
      coord_fixed() +
      scale_x_discrete(
        limits = feature_order,
        labels = label_map
      ) +
      xlab(NULL) + ylab(NULL) +
      labs(fill="Coefficient") +
      scale_fill_gradient2(
        low = "red", mid = "white", high = "green4",
        midpoint = 0,
        limits = c(-1, 1),
        na.value = "white",
        oob = scales::squish
      ) + 
      theme(
        legend.position = "top",
        legend.justification = "right",
        axis.text.x = ggtext::element_markdown(angle = 45, hjust = 1, vjust = 1),
        panel.background = element_blank()
      )      
      
    plot_filename <- file.path(outdir, paste0("Heatmap_", "ClinicalOutcome", "_in_AllPatients", "-", paste0(unique(df_Combined$metadata), collapse = "_"), "-BOLD", ".png"))
    ggsave(plot_filename, plot = p_heatmap_BOLD, width = 24, height = 10, units = "cm") # 5 rows
    
    #### Without Coefficient
    p_heatmap_BOLD_NoLegend <- p_heatmap_BOLD + theme(legend.position = "none")
    plot_filename <- file.path(outdir, paste0("Heatmap_", "ClinicalOutcome", "_in_AllPatients", "-", paste0(unique(df_Combined$metadata), collapse = "_"), "-BOLD", "-NoLegend", ".png"))
    ggsave(plot_filename, plot = p_heatmap_BOLD_NoLegend, width = 24, height = 8, units = "cm") # 5 rows
    
    
      
    
  
  
  
  
  
  
  