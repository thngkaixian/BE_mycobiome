library(ggplot2)


# --------------------------------------------------------------------------------------------------------------
## Directories 
# --------------------------------------------------------------------------------------------------------------

Project_name <- "MYCOBIOME_AsiaBRIDGE"
Library_Number <- "Library1234"
taxLevel <- "Genus"

Filter <- "Filtered"

## Directories 
wkdir <- file.path("/Users/thngkaixian/PhD_RESEARCH/RESEARCH_Projects", Project_name)
indir <- file.path(wkdir, "SequencingData", Library_Number, "Analysis")


  ## Maaslin results directories
  # --------------------------------------------------------------------------------------------------------------
  Co_variate <- "Country"

  Maaslin_Model <- "ModelCPLM_NormalisationNONE_TransformationNONE"

  ## In-directories 
  indir_Maaslin_Results <- file.path(indir, '03 Maaslin_ByRelAbundance_Fixed_Random_Effects', Filter)
  ## Create out directory
  outdir <- file.path(indir_Maaslin_Results, "00_COMPILED_PLOTS", 
                      paste0(Maaslin_Model, "_RESULTS"))
  if (!dir.exists(outdir)) try(dir.create(outdir, recursive= TRUE), silent= TRUE)
    
    
# --------------------------------------------------------------------------------------------------------------
## Import data
# --------------------------------------------------------------------------------------------------------------

qval_cutoff <- 0.1


  ## ASIA (SG/KL)
  # --------------------------------------------------------------------------------------------------------------
  Asia_FrequentExacerbator <- read.delim(file.path(indir_Maaslin_Results, paste0('FrequentExacerbator', '_RandomEffectAs_Country', '_in_Continent_Asia'), Maaslin_Model, 'all_results.tsv'),  header = TRUE, sep = "\t")
  Asia_MRC <- read.delim(file.path(indir_Maaslin_Results, paste0('MRC_Score', '_RandomEffectAs_Country', '_in_Continent_Asia'), Maaslin_Model, 'all_results.tsv'),  header = TRUE, sep = "\t" )
  Asia_BSI <- read.delim(file.path(indir_Maaslin_Results, paste0('BSI', '_RandomEffectAs_Country', '_in_Continent_Asia'), Maaslin_Model, 'all_results.tsv'), header = TRUE, sep = "\t")
  Asia_FEV1 <- read.delim(file.path(indir_Maaslin_Results, paste0('FEV1', '_RandomEffectAs_Country', '_in_Continent_Asia'), Maaslin_Model, 'all_results.tsv'), header = TRUE, sep = "\t")
  Asia_Hospitalisation <- read.delim(file.path(indir_Maaslin_Results, paste0('Hospitalisation', '_RandomEffectAs_Country', '_in_Continent_Asia'), Maaslin_Model, 'all_results.tsv'), header = TRUE, sep = "\t")
  
  Asia <- rbind(Asia_FrequentExacerbator,
                Asia_MRC,
                Asia_BSI, Asia_FEV1,
                Asia_Hospitalisation)
  Asia.qval <- subset(Asia, qval < qval_cutoff)
    
    
  ## European Regions-----NWE
  # --------------------------------------------------------------------------------------------------------------
  NWE_FrequentExacerbator <- read.delim(file.path(indir_Maaslin_Results, paste0('FrequentExacerbator', '_RandomEffectAs_Country', '_in_Region_NWE'), Maaslin_Model, 'all_results.tsv'),  header = TRUE, sep = "\t")
  NWE_MRC <- read.delim(file.path(indir_Maaslin_Results, paste0('MRC_Score', '_RandomEffectAs_Country', '_in_Region_NWE'), Maaslin_Model, 'all_results.tsv'),  header = TRUE, sep = "\t" )
  NWE_BSI <- read.delim(file.path(indir_Maaslin_Results, paste0('BSI', '_RandomEffectAs_Country', '_in_Region_NWE'), Maaslin_Model, 'all_results.tsv'), header = TRUE, sep = "\t")
  NWE_FEV1 <- read.delim(file.path(indir_Maaslin_Results, paste0('FEV1', '_RandomEffectAs_Country', '_in_Region_NWE'), Maaslin_Model, 'all_results.tsv'), header = TRUE, sep = "\t")
  NWE_Hospitalisation <- read.delim(file.path(indir_Maaslin_Results, paste0('Hospitalisation', '_RandomEffectAs_Country', '_in_Region_NWE'), Maaslin_Model, 'all_results.tsv'), header = TRUE, sep = "\t")
  
  NWE <- rbind(NWE_FrequentExacerbator,
              NWE_MRC,
               NWE_BSI, NWE_FEV1,
                NWE_Hospitalisation)
  NWE.qval <- subset(NWE, qval < qval_cutoff)
  
  
  ## European Regions-----SE
  # --------------------------------------------------------------------------------------------------------------
  SE_FrequentExacerbator <- read.delim(file.path(indir_Maaslin_Results, paste0('FrequentExacerbator', '_RandomEffectAs_Country', '_in_Region_SE'), Maaslin_Model, 'all_results.tsv'),  header = TRUE, sep = "\t")
  SE_MRC <- read.delim(file.path(indir_Maaslin_Results, paste0('MRC_Score', '_RandomEffectAs_Country', '_in_Region_SE'), Maaslin_Model, 'all_results.tsv'),  header = TRUE, sep = "\t" )
  SE_BSI <- read.delim(file.path(indir_Maaslin_Results, paste0('BSI', '_RandomEffectAs_Country', '_in_Region_SE'), Maaslin_Model, 'all_results.tsv'), header = TRUE, sep = "\t")
  SE_FEV1 <- read.delim(file.path(indir_Maaslin_Results, paste0('FEV1', '_RandomEffectAs_Country', '_in_Region_SE'), Maaslin_Model, 'all_results.tsv'), header = TRUE, sep = "\t")
  SE_Hospitalisation <- read.delim(file.path(indir_Maaslin_Results, paste0('Hospitalisation', '_RandomEffectAs_Country', '_in_Region_SE'), Maaslin_Model, 'all_results.tsv'), header = TRUE, sep = "\t")
  
  SE <- rbind(SE_FrequentExacerbator,
              SE_MRC,
              SE_BSI, SE_FEV1,
              SE_Hospitalisation)
  SE.qval <- subset(SE, qval < qval_cutoff)
  
  
  ## European Regions-----UK
  # --------------------------------------------------------------------------------------------------------------
  UK_FrequentExacerbator <- read.delim(file.path(indir_Maaslin_Results, paste0('FrequentExacerbator', '_RandomEffectAs_Country', '_in_Region_UK'), Maaslin_Model, 'all_results.tsv'),  header = TRUE, sep = "\t")
  UK_MRC <- read.delim(file.path(indir_Maaslin_Results, paste0('MRC_Score', '_RandomEffectAs_Country', '_in_Region_UK'), Maaslin_Model, 'all_results.tsv'),  header = TRUE, sep = "\t" )
  UK_BSI <- read.delim(file.path(indir_Maaslin_Results, paste0('BSI', '_RandomEffectAs_Country', '_in_Region_UK'), Maaslin_Model, 'all_results.tsv'), header = TRUE, sep = "\t")
  UK_FEV1 <- read.delim(file.path(indir_Maaslin_Results, paste0('FEV1', '_RandomEffectAs_Country', '_in_Region_UK'), Maaslin_Model, 'all_results.tsv'), header = TRUE, sep = "\t")
  UK_Hospitalisation <- read.delim(file.path(indir_Maaslin_Results, paste0('Hospitalisation', '_RandomEffectAs_Country', '_in_Region_UK'), Maaslin_Model, 'all_results.tsv'), header = TRUE, sep = "\t")
  
  UK <- rbind(UK_FrequentExacerbator,
              UK_MRC,
              UK_BSI, UK_FEV1,
              UK_Hospitalisation)
  UK.qval <- subset(UK, qval < qval_cutoff)

  


# --------------------------------------------------------------------------------------------------------------
## IF Plotting ALL FUNGI --- Ensure there are 21 fungi
## - Append with "excluded" / "missing" fungi
# --------------------------------------------------------------------------------------------------------------
  
Ascomycota <- c("Candida", "Saccharomyces", "Penicillium", "Clavispora", "Aspergillus", 
                "Microidium", 
                "Cladosporium", "Fusarium", "Mycosphaerella", "Nakaseomyces", "Chaenotheca")
Basidiomycota <- c("Schizophyllum", "Malassezia", "Wallemia", 
                   "Cutaneotrichosporon", "Trametes", "Lentinus", "Cryptococcus", 
                   "Heterobasidion", "Rhodotorula", "Filobasidium")
AllFungi <- c(Ascomycota, Basidiomycota)
 

  ##-----NWE (All clinical Variables)
  # --------------------------------------------------------------------------------------------------------------
  df_NWE_FrequentExacerbator_ExcludedFungi <- data.frame(feature = setdiff(AllFungi, NWE_FrequentExacerbator$feature), metadata = "FrequentExacerbator", value = "FreqExac", coef = NA , stderr = NA, N = NA, N.not.0 = NA, pval = NA, qval = NA)
  df_NWE_MRC_ExcludedFungi <- data.frame(feature = setdiff(AllFungi, NWE_MRC$feature), metadata = "MRC_score", value = "MRC_score", coef = NA , stderr = NA, N = NA, N.not.0 = NA, pval = NA, qval = NA)
  df_NWE_BSI_ExcludedFungi <- data.frame(feature = setdiff(AllFungi, NWE_BSI$feature), metadata = "BSI", value = "BSI", coef = NA , stderr = NA, N = NA, N.not.0 = NA, pval = NA, qval = NA)
  df_NWE_FEV1_ExcludedFungi <- data.frame(feature = setdiff(AllFungi, NWE_FEV1$feature), metadata = "FEV1_percent_predicted", value = "FEV1_percent_predicted", coef = NA , stderr = NA, N = NA, N.not.0 = NA, pval = NA, qval = NA)
  df_NWE_Hospitalisation_ExcludedFungi <- data.frame(feature = setdiff(AllFungi, NWE_Hospitalisation$feature), metadata = "Hospitalised_for_bronchiectasis_YesNo", value = "Yes", coef = NA , stderr = NA, N = NA, N.not.0 = NA, pval = NA, qval = NA)
  
  NWE <- rbind(NWE, 
               df_NWE_FrequentExacerbator_ExcludedFungi, 
               df_NWE_MRC_ExcludedFungi,
               df_NWE_BSI_ExcludedFungi,
               df_NWE_FEV1_ExcludedFungi,
               df_NWE_Hospitalisation_ExcludedFungi)
  NWE <- NWE[order(NWE$metadata),] #Order dataframe
  
  
  ##-----SE (All clinical Variables)
  # --------------------------------------------------------------------------------------------------------------
  df_SE_FrequentExacerbator_ExcludedFungi <- data.frame(feature = setdiff(AllFungi, SE_FrequentExacerbator$feature), metadata = "FrequentExacerbator", value = "FreqExac", coef = NA , stderr = NA, N = NA, N.not.0 = NA, pval = NA, qval = NA)
  df_SE_MRC_ExcludedFungi <- data.frame(feature = setdiff(AllFungi, SE_MRC$feature), metadata = "MRC_score", value = "MRC_score", coef = NA , stderr = NA, N = NA, N.not.0 = NA, pval = NA, qval = NA)
  df_SE_BSI_ExcludedFungi <- data.frame(feature = setdiff(AllFungi, SE_BSI$feature), metadata = "BSI", value = "BSI", coef = NA , stderr = NA, N = NA, N.not.0 = NA, pval = NA, qval = NA)
  df_SE_FEV1_ExcludedFungi <- data.frame(feature = setdiff(AllFungi, SE_FEV1$feature), metadata = "FEV1_percent_predicted", value = "FEV1_percent_predicted", coef = NA , stderr = NA, N = NA, N.not.0 = NA, pval = NA, qval = NA)
  df_SE_Hospitalisation_ExcludedFungi <- data.frame(feature = setdiff(AllFungi, SE_Hospitalisation$feature), metadata = "Hospitalised_for_bronchiectasis_YesNo", value = "Yes", coef = NA , stderr = NA, N = NA, N.not.0 = NA, pval = NA, qval = NA)
  
  SE <- rbind(SE,
               df_SE_FrequentExacerbator_ExcludedFungi, 
               df_SE_MRC_ExcludedFungi,
               df_SE_BSI_ExcludedFungi,
               df_SE_FEV1_ExcludedFungi,
               df_SE_Hospitalisation_ExcludedFungi)
  SE <- SE[order(SE$metadata),] #Order dataframe


# --------------------------------------------------------------------------------------------------------------
## Combine data into dataframe for plotting 
# --------------------------------------------------------------------------------------------------------------

  ## Create column for geographic region  
  # --------------------------------------------------------------------------------------------------------------
  Asia$Geography <- "SEA" #"SG_KL"
  NWE$Geography <- "NWE"
  SE$Geography <- "SE"
  UK$Geography <- "UK"
  

  ## Combine data 
  # --------------------------------------------------------------------------------------------------------------
  df_Results_COMPILED <- rbind(Asia, NWE, SE, UK)
  
  
  ## Factor
  # --------------------------------------------------------------------------------------------------------------
  df_Results_COMPILED$Geography <- factor(df_Results_COMPILED$Geography, 
                                          levels = c("SEA", "NWE", "SE", "UK"))

  df_Results_COMPILED$feature <- factor(df_Results_COMPILED$feature, levels = rev(levels(factor(df_Results_COMPILED$feature)))) #Reverse the order
  
  df_Results_COMPILED$metadata <- factor(df_Results_COMPILED$metadata, 
                                         levels = c("FrequentExacerbator", "Hospitalised_for_bronchiectasis_YesNo", "BSI", "MRC_score", "FEV1_percent_predicted") ,
                                          labels = c("Frequent Exacerbator", "Severe Exacerbation", "BSI", "MRC score", "FEV1 (% predicted)"))

  
  ## Replace non-significant qval with NA
  # --------------------------------------------------------------------------------------------------------------
  df_Results_COMPILED$qval_Sigf <- sapply(df_Results_COMPILED$qval, function(x){ifelse(x>=qval_cutoff | is.na(x), NA, x)})
  
  
  ## Replace coeff with 0 if qval non-significant
  # --------------------------------------------------------------------------------------------------------------
  df_Results_COMPILED[,"coef_Sigf"] <- df_Results_COMPILED$coef
  df_Results_COMPILED[which(df_Results_COMPILED$qval >= qval_cutoff) , "coef_Sigf"] <- 0
  sort(df_Results_COMPILED$coef)
  sort(df_Results_COMPILED$coef_Sigf)
  
  
  ## Add Ascomycota/Basidiomycota classifier 
  # --------------------------------------------------------------------------------------------------------------
  df_Results_COMPILED$Phylum <- sapply(df_Results_COMPILED$feature, function(x){ if(x %in% Ascomycota ) "Ascomycota" else if (x %in% Basidiomycota) "Basidiomycota" })

    
  ## Indicate Positive/Negative correlation 
  # --------------------------------------------------------------------------------------------------------------
  df_Results_COMPILED$coef_PositiveNegative <- ifelse(df_Results_COMPILED$coef > 0, "Positive", "Negative")


  ## Calculate 95% Confidence Intervals
  # --------------------------------------------------------------------------------------------------------------
  df_Results_COMPILED$ci_lower <- df_Results_COMPILED$coef - 1.96 * df_Results_COMPILED$stderr
  df_Results_COMPILED$ci_upper <- df_Results_COMPILED$coef + 1.96 * df_Results_COMPILED$stderr


  ##----Export dataframe for viewing
  # df_Results_COMPILED$qval_asterisks <- sapply(df_Results_COMPILED$qval_Sigf, function(x){ if(x<0.001 & !is.na(x)) "***" else if(x<0.01 & !is.na(x)) "**" else if (x<0.1 & !is.na(x)) "*" else if(is.na(x)) "NS" })
  # df_Results_COMPILED_filename <- file.path(outdir, paste0("COMPILED_MaaslinResults_", "qval", qval_cutoff, "_", "ClinicalOutcome", "_", "FacetClinicalVariables", "_Across_", "Regions", "-",
  #                                                          paste0(unique(df_Results_COMPILED$metadata), collapse = "_"), ".csv"))
  # write.csv(df_Results_COMPILED, df_Results_COMPILED_filename)
  
  
  
  
#--------------------------------------------------------------------------------------------------------------
## plot 
# --------------------------------------------------------------------------------------------------------------

## Select dataset for plotting ##
DF_PLOT <- df_Results_COMPILED
sort(DF_PLOT$coef)


  ## Heatmap
  # --------------------------------------------------------------------------------------------------------------
  p_heatmap <-
    ggplot(DF_PLOT, aes(y = feature, x = Geography, fill = coef_Sigf)) + 
    geom_tile(color = "black") +
    facet_grid(.~metadata, scales = "free_y" ) +
    xlab(NULL) +
    ylab(NULL) +
    labs(fill="Coefficient") +
    scale_fill_gradient2(low = "red", mid = "white", high = "green4", 
                         midpoint = 0,          
                         limits = c(-2, 2),    
                         na.value = "white",
                         oob = scales::squish 
                         ) +
    theme(legend.position = "top",
          legend.justification = "right",
          panel.background = element_blank(), 
          axis.text.y = element_text(face = "italic", size = 11), #Fungi
          axis.text.x = element_text(size = 11),
          strip.text = element_text(size = 12)) 
  p_heatmap
  
  ## Export
  p_heatmap_filename <- file.path(outdir, paste0("Heatmap_", "qval", qval_cutoff, "_", "FacetClinicalVariables", "_Across_", "Regions", "-",
                                            paste0(unique(df_Results_COMPILED$metadata), collapse = "-"), ".png"))
  ggsave(p_heatmap_filename, plot = p_heatmap, width = 34, height = 18, units = "cm") # 5 facets
  
  

  
