library(ggplot2)
# library(vioplot)


# --------------------------------------------------------------------------------------------------------------
## Directories 
# --------------------------------------------------------------------------------------------------------------

Project_name <- "MYCOBIOME_AsiaBRIDGE"

taxLevel <- "Genus"

## Directories 
wkdir <- file.path("/Users/thngkaixian/PhD_RESEARCH/RESEARCH_Projects", Project_name)
indir <- file.path(wkdir, "SequencingData", "FINAL_ANALYSIS", "Analysis")


  ## Maaslin results directories
  # --------------------------------------------------------------------------------------------------------------
  Maaslin_Model <- "ModelCPLM_NormalisationNONE_TransformationNONE"

  ## In-directories
  indir_Maaslin_Results <- file.path(indir, 'Figure_3C_MatchedCohort')
  ## Create out directory
  outdir <- file.path(indir_Maaslin_Results, "00_COMPILED_PLOTS", paste0(Maaslin_Model, "_RESULTS"))
  if (!dir.exists(outdir)) try(dir.create(outdir, recursive= TRUE), silent= TRUE)
  
  
# --------------------------------------------------------------------------------------------------------------
## Import data
# --------------------------------------------------------------------------------------------------------------

qval_cutoff <- 0.1

  ## With Covariate (Random Effect)
  # --------------------------------------------------------------------------------------------------------------
  Co_variate <- "Country"
  

    #### Asia ####
    # --------------------------------------------------------------------------------------------------------------
    df_Asia_FrequentExacerbator <-  read.delim(file.path(indir_Maaslin_Results, paste0('FrequentExacerbator', '_RandomEffectAs_', Co_variate, '_in_Continent_', "Asia"), Maaslin_Model, 'all_results.tsv'), header = TRUE, sep = "\t")
    df_Asia_BSI <-                  read.delim(file.path(indir_Maaslin_Results, paste0('BSI', '_RandomEffectAs_', Co_variate, '_in_Continent_', "Asia"), Maaslin_Model, 'all_results.tsv'), header = TRUE, sep = "\t")
    df_Asia_FEV1 <-                 read.delim(file.path(indir_Maaslin_Results, paste0('FEV1', '_RandomEffectAs_', Co_variate, '_in_Continent_', "Asia"), Maaslin_Model, 'all_results.tsv'), header = TRUE, sep = "\t")
    df_Asia_MRC <-                  read.delim(file.path(indir_Maaslin_Results, paste0('MRC_Score', '_RandomEffectAs_', Co_variate, '_in_Continent_', "Asia"), Maaslin_Model, 'all_results.tsv'),  header = TRUE, sep = "\t" )
    df_Asia_HospitalisationYesNo <- read.delim(file.path(indir_Maaslin_Results, paste0('HospitalisationYesNo', '_RandomEffectAs_', Co_variate, '_in_Continent_', "Asia"), Maaslin_Model, 'all_results.tsv'), header = TRUE, sep = "\t")
    
    df_Asia <- rbind(df_Asia_FrequentExacerbator,
                  df_Asia_BSI,
                  df_Asia_FEV1, df_Asia_MRC,
                  df_Asia_HospitalisationYesNo)
    df_Asia.qval <- subset(df_Asia, qval < qval_cutoff)  
    
    
    #### Europe-MatchedAGB ####
    # --------------------------------------------------------------------------------------------------------------
    df_EuropeAGB_FrequentExacerbator <-  read.delim(file.path(indir_Maaslin_Results, paste0('FrequentExacerbator', '_RandomEffectAs_', Co_variate, '_in_Continent_', "Europe-MatchedAGB"), Maaslin_Model, 'all_results.tsv'), header = TRUE, sep = "\t")
    df_EuropeAGB_BSI <-                  read.delim(file.path(indir_Maaslin_Results, paste0('BSI', '_RandomEffectAs_', Co_variate, '_in_Continent_', "Europe-MatchedAGB"), Maaslin_Model, 'all_results.tsv'), header = TRUE, sep = "\t")
    df_EuropeAGB_FEV1 <-                 read.delim(file.path(indir_Maaslin_Results, paste0('FEV1', '_RandomEffectAs_', Co_variate, '_in_Continent_', "Europe-MatchedAGB"), Maaslin_Model, 'all_results.tsv'), header = TRUE, sep = "\t")
    df_EuropeAGB_MRC <-                  read.delim(file.path(indir_Maaslin_Results, paste0('MRC_Score', '_RandomEffectAs_', Co_variate, '_in_Continent_', "Europe-MatchedAGB"), Maaslin_Model, 'all_results.tsv'),  header = TRUE, sep = "\t" )
    df_EuropeAGB_HospitalisationYesNo <- read.delim(file.path(indir_Maaslin_Results, paste0('HospitalisationYesNo', '_RandomEffectAs_', Co_variate, '_in_Continent_', "Europe-MatchedAGB"), Maaslin_Model, 'all_results.tsv'), header = TRUE, sep = "\t")
  
    df_EuropeAGB <- rbind(,df_EuropeAGB_FrequentExacerbator,
                       df_EuropeAGB_BSI,
                       df_EuropeAGB_FEV1, df_EuropeAGB_MRC,
                       df_EuropeAGB_HospitalisationYesNo)
    df_EuropeAGB.qval <- subset(df_EuropeAGB, qval < qval_cutoff)  
    
    
    #### Europe-MatchedAGEL ####
    # --------------------------------------------------------------------------------------------------------------
    df_EuropeAGEL_FrequentExacerbator <-  read.delim(file.path(indir_Maaslin_Results, paste0('FrequentExacerbator', '_RandomEffectAs_', Co_variate, '_in_Continent_', "Europe-MatchedAGEL"), Maaslin_Model, 'all_results.tsv'), header = TRUE, sep = "\t")
    df_EuropeAGEL_BSI <-                  read.delim(file.path(indir_Maaslin_Results, paste0('BSI', '_RandomEffectAs_', Co_variate, '_in_Continent_', "Europe-MatchedAGEL"), Maaslin_Model, 'all_results.tsv'), header = TRUE, sep = "\t")
    df_EuropeAGEL_FEV1 <-                 read.delim(file.path(indir_Maaslin_Results, paste0('FEV1', '_RandomEffectAs_', Co_variate, '_in_Continent_', "Europe-MatchedAGEL"), Maaslin_Model, 'all_results.tsv'), header = TRUE, sep = "\t")
    df_EuropeAGEL_MRC <-                  read.delim(file.path(indir_Maaslin_Results, paste0('MRC_Score', '_RandomEffectAs_', Co_variate, '_in_Continent_', "Europe-MatchedAGEL"), Maaslin_Model, 'all_results.tsv'),  header = TRUE, sep = "\t" )
    df_EuropeAGEL_HospitalisationYesNo <- read.delim(file.path(indir_Maaslin_Results, paste0('HospitalisationYesNo', '_RandomEffectAs_', Co_variate, '_in_Continent_', "Europe-MatchedAGEL"), Maaslin_Model, 'all_results.tsv'), header = TRUE, sep = "\t")
    
    df_EuropeAGEL <- rbind(df_EuropeAGEL_FrequentExacerbator,
                           df_EuropeAGEL_BSI,
                           df_EuropeAGEL_FEV1, df_EuropeAGEL_MRC,
                           df_EuropeAGEL_HospitalisationYesNo)
    df_EuropeAGEL.qval <- subset(df_EuropeAGEL, qval < qval_cutoff)  
    
    
    
    
# --------------------------------------------------------------------------------------------------------------
## Combine data into dataframe for plotting 
# --------------------------------------------------------------------------------------------------------------

  ## Create column - Indicate Geography 
  # --------------------------------------------------------------------------------------------------------------
  
  #### using qval < 0.1 ####
  df_Asia.qval$Geography <- "Asia"
  df_EuropeAGB.qval$Geography <- "Europe_AGB"
  df_EuropeAGEL.qval$Geography <- "Europe_AGEL"
  
  
  ## Combine data 
  # --------------------------------------------------------------------------------------------------------------
  df_Results_COMPILED <- rbind(df_Asia.qval, df_EuropeAGB.qval, df_EuropeAGEL.qval)
  # df_Results_COMPILED <- df_AllPatients_AGB
  # df_Results_COMPILED <- df_AllPatients_AGEL
  
  
    ## Replace non-significant qval with NA
    # --------------------------------------------------------------------------------------------------------------
    qval_cutoff <- 0.1
    df_Results_COMPILED$qval_Sigf <- sapply(df_Results_COMPILED$qval, function(x){ifelse(x>=qval_cutoff | is.na(x), NA, x)})
    
    ## Replace coeff with 0 if qval non-significant
    # --------------------------------------------------------------------------------------------------------------
    df_Results_COMPILED[,"coef_Sigf"] <- df_Results_COMPILED$coef
    df_Results_COMPILED[which(df_Results_COMPILED$qval >= qval_cutoff) , "coef_Sigf"] <- 0
    sort(df_Results_COMPILED$coef)
    sort(df_Results_COMPILED$coef_Sigf)
    
  
  ## Factor 
  # --------------------------------------------------------------------------------------------------------------
  unique(df_Results_COMPILED$metadata)
  df_Results_COMPILED$metadata <- factor(df_Results_COMPILED$metadata, levels = c("Exacerbations",
                                                                                  "FrequentExacerbator",
                                                                                  "Hospitalised_for_bronchiectasis_YesNo",
                                                                                  "Exacerbator_YesNo",
                                                                                  "BSI", 
                                                                                  "FEV1_percent_predicted",
                                                                                  "MRC_score" # "No..of.lobes.involved", "Severe_NonSevere"
                                                                                  )#,
                                                                       # labels = c("Exacerbation\nfrequency", 
                                                                       #            "Frequent Exacerbator",
                                                                       #            "Hospitalisation",
                                                                       #            "Exacerbator",
                                                                       #            "BSI", 
                                                                       #            "FEV1\n(%\npredicted)",
                                                                       #            "MRC\nscore"# "No..of.lobes.involved",
                                                                       #            )
                                         )
  # df_Results_COMPILED$metadata <- factor(df_Results_COMPILED$metadata, levels = rev(levels(factor(df_Results_COMPILED$metadata)))) #so that order starts from the top
  
  
  df_Results_COMPILED$feature <- factor(df_Results_COMPILED$feature, levels = rev(levels(factor(df_Results_COMPILED$feature)))) #so that order starts from the top
  df_Results_COMPILED$Geography <- factor(df_Results_COMPILED$Geography,
                                          levels=c("Asia", "Europe_AGB", "Europe_AGEL")#,
                                          # labels = c("Asia", "Europe\n(Severity-Matched)", "Europe\n(Exacerbation and\nFEV1-Matched)")
                                          )

  
  ## Indicate Positive/Negative correlation 
  # --------------------------------------------------------------------------------------------------------------
  df_Results_COMPILED$coef_PositiveNegative <- ifelse(df_Results_COMPILED$coef > 0, "Positive", "Negative")


  ## Calculate 95% Confidence Intervals
  # --------------------------------------------------------------------------------------------------------------
  df_Results_COMPILED$ci_lower <- df_Results_COMPILED$coef - 1.96 * df_Results_COMPILED$stderr
  df_Results_COMPILED$ci_upper <- df_Results_COMPILED$coef + 1.96 * df_Results_COMPILED$stderr
  
  
  ## Calculate association_score
  # --------------------------------------------------------------------------------------------------------------
  # df_Results_COMPILED$log_qval_coef <- -log10(df_Results_COMPILED$qval_Sigf)*sign(df_Results_COMPILED$coef)  
  # df_Results_COMPILED$log_qval_coef <- -log10(df_Results_COMPILED$qval)*sign(df_Results_COMPILED$coef)  
  # sort(df_Results_COMPILED$log_qval_coef)
  # vioplot(df_Results_COMPILED$log_qval_coef)
  # var(df_Results_COMPILED$log_qval_coef)
  
  
      ##----Export dataframe for viewing
      df_Results_COMPILED$qval_asterisks <- sapply(df_Results_COMPILED$qval_Sigf, function(x){ if(x<0.001 & !is.na(x)) "***" else if(x<0.01 & !is.na(x)) "**" else if (x<0.1 & !is.na(x)) "*" else if(is.na(x)) "NS" })
      df_Results_COMPILED_filename <- file.path(outdir, paste0("COMPILED_MaaslinResults_", "ClinicalOutcome", "_in_Asia_EuropeAGB_EuropeAGEL", "-", paste0(unique(df_Results_COMPILED$metadata), collapse = "_"), ".csv"))
      # df_Results_COMPILED_filename <- file.path(outdir, paste0("COMPILED_MaaslinResults_", "ClinicalOutcome", "_in_AllPatients_AGB", "-", paste0(unique(df_Results_COMPILED$metadata), collapse = "_"), ".csv"))
      # df_Results_COMPILED_filename <- file.path(outdir, paste0("COMPILED_MaaslinResults_", "ClinicalOutcome", "_in_AllPatients_AGEL", "-", paste0(unique(df_Results_COMPILED$metadata), collapse = "_"), ".csv"))
      write.csv(df_Results_COMPILED, df_Results_COMPILED_filename)
      
  
  
#--------------------------------------------------------------------------------------------------------------
## plot 
# --------------------------------------------------------------------------------------------------------------

## Select dataset for plotting ##
DF_PLOT <- df_Results_COMPILED

  ## Forest plot  
  # --------------------------------------------------------------------------------------------------------------

  # New facet label names
  Clinical.labs <- c("Frequent Exacerbator", "Severe Exacerbation", "BSI", "FEV1\n(%\npredicted)", "MRC\nscore") 
  names(Clinical.labs) <- c("FrequentExacerbator", "Hospitalised_for_bronchiectasis_YesNo", "BSI", "FEV1_percent_predicted", "MRC_score") 

  Geography.labs <- c("South-East Asian", "European:\nAge, Sex, Severity-Matched", "European:\nAge, Sex, Exacerbation,\nFEV1-Matched")
  names(Geography.labs) <- c("Asia", "Europe_AGB", "Europe_AGEL")
  
  
  ## Plot
  p_Forest <- ggplot(DF_PLOT, aes(x = coef, y = feature)) +
    geom_point(aes(color = coef_PositiveNegative), shape = 16) +
    geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper, color = coef_PositiveNegative), width = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed") + 
    ggtitle(paste0("Random Effect: ", Co_variate)) +
    facet_grid(metadata~Geography, scales ="free_y", space = "free",
               labeller = labeller(metadata = Clinical.labs, Geography = Geography.labs)) + #-----Rename Facet Labels
    xlab("Coefficient (95% CI)") +
    ylab(NULL) +
    scale_color_manual(values = c("Positive" = "green4", "Negative" = "red3")) + 
    theme(legend.position = "none",
          panel.border = element_blank(),
          plot.background = element_rect(fill='transparent', color=NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text.y = element_text(face = "italic"),
          strip.text.x = element_text(size = 11), ##---Text size of X-axis facet labels
          # strip.text.y = element_text(size = 11),  ##---Text size of Y-axis facet labels
          strip.background = element_rect(fill = "white", colour = "black", size = 1)
    )
  p_Forest
  
  
    #### Plot --- BLANK OUT SPECIFIC PANEL(S) by Plot a white box over it.
    # --------------------------------------------------------------------------------------------------------------
  
    ## Panels to blank out --- To prevent re-ordering
    blank_rect <- data.frame(metadata  = factor(c("FrequentExacerbator", "FEV1_percent_predicted", "BSI"),
                                                levels = levels(DF_PLOT$metadata)),
                             Geography = factor(c("Europe_AGEL", "Europe_AGEL", "Europe_AGB"),
                                                levels = levels(DF_PLOT$Geography)))
  
    # blank_rect$label <- c("N.A. as matched\nfor exacerbation", "N.A. as matched\nfor lung function", "N.A. as matched\nfor severity")  # match number of boxes

    p_Forest_BLANKED <- p_Forest + 
      geom_rect(data = blank_rect,
                aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
                inherit.aes = FALSE,
                fill = "white",   # makes it look blank
                color = NA)
    p_Forest_BLANKED
  

    #### Plot --- Remove title
    # --------------------------------------------------------------------------------------------------------------
    p_Forest_without_Title <- p_Forest + ggtitle(NULL)
    p_Forest_without_Title
    
    p_Forest_BLANKED_without_Title <- p_Forest_BLANKED + ggtitle(NULL)
    p_Forest_BLANKED_without_Title
    
  
    ##---Export 
    plot_filename <- file.path(outdir, paste0("COMPILED_", "ForestPlot", "_", "qval", qval_cutoff, "_RandomEffectAs_", Co_variate, 
                                              "_FacetBy_Asia_EuropeanMatchedCohorts", ".png"))
    ggsave(plot_filename, plot = p_Forest, bg='transparent', width = 25, height = 22, units = "cm")
  
    
    p_Forest_without_Title_filename <- file.path(outdir, paste0("COMPILED_", "ForestPlot", "_", "qval", qval_cutoff, "_RandomEffectAs_", Co_variate,
                                                                "_FacetBy_Asia_EuropeanMatchedCohorts", "-V2-NoTitle", ".png"))
    ggsave(p_Forest_without_Title_filename, plot = p_Forest_without_Title, bg='transparent', width = 20, height = 25, units = "cm")

    
    
    p_Forest_BLANKED_without_Title_filename <- file.path(outdir, paste0("COMPILED_", "ForestPlot", "_", "qval", qval_cutoff, "_RandomEffectAs_", Co_variate,
                                                                "_FacetBy_Asia_EuropeanMatchedCohorts", "-V2-NoTitle-BLANKEDoutFacets", ".png"))
    ggsave(p_Forest_BLANKED_without_Title_filename, plot = p_Forest_BLANKED_without_Title, bg='transparent', width = 20, height = 25, units = "cm")
    
    
        
  
  
  
  
  
  
  
  
  
  
