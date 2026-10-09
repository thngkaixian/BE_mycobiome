## Maaslin: Mycobiome ~ Clinical outcome

library(phyloseq)
library(vegan)
library(Maaslin2)


# --------------------------------------------------------------------------------------------------------------
## Directories 
# --------------------------------------------------------------------------------------------------------------

Project_name <- "MYCOBIOME_AsiaBRIDGE"

taxLevel <- "Genus"

wkdir <- file.path("/Users/thngkaixian/PhD_RESEARCH/RESEARCH_Projects", Project_name, "SequencingData")
indir_ps <- file.path(wkdir, "FINAL_ANALYSIS")


  ## Non-matched cohort
  outdir <- file.path(indir_ps, "Analysis", "Figure_3B_ByRegion")
  
  ## Matched cohort
  # outdir <- file.path(wkdir, "SequencingData", Library_Number, "Analysis", "03 Maaslin_ByRelAbundance_Fixed_Random_Effects-MatchedCohort")
  
  #Create out directory
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
  
  
# ==============================================================================================================
# SUBSET into subgroups 
# ==============================================================================================================

  ## Continent --- Asia | Europe 
  # ------------------------------------------------------------------------------------------
  ps.prop_Asia <- subset_samples(ps.prop, Continent=="Asia")
  
  ## Regions --- NWE | SE | UK 
  # ------------------------------------------------------------------------------------------
  ps.prop_NWE <- subset_samples(ps.prop, Region=="NWE")
  ps.prop_SE <- subset_samples(ps.prop, Region=="SE")
  ps.prop_UK <- subset_samples(ps.prop, Region=="UK")

  
# ==============================================================================================================
## Assign data to analyse
# ==============================================================================================================
  

  ##-----Continent
  ps_subgroup <- ps.prop_Asia

  ##-----Region 
  ps_subgroup <- ps.prop_NWE
  ps_subgroup <- ps.prop_SE
  ps_subgroup <- ps.prop_UK
  

# ----------------------------------------------------------------------------------------------------------------------------------------
## Preparing metadata for input to maaslin
# ----------------------------------------------------------------------------------------------------------------------------------------  
  
## Extract metadata and OTU table
metadata <- data.frame(sample_data(ps_subgroup))
OTU <- data.frame(otu_table(ps_subgroup))

  #### INSANITY CHECK ####
  table(metadata$Continent) ; table(metadata$Country) ; table(metadata$Region)
  unique(metadata$Matching_AgeGenderBSI) #To Check Matched Patients
  unique(metadata$Matching_AgeGenderExacerbatorStatusFEV1)  #To Check Matched Patients


# ==============================================================================================================
## run masslin (PARAMETERS 2) - CPLM
# ==============================================================================================================

  ## update out-directory with parameters
  outdir_MaaslinParameters <- file.path(outdir,
                                        paste0("Model", "CPLM", "_",
                                               "Normalisation", "NONE", "_",
                                               "Transformation", "NONE"))
  if (!dir.exists(outdir_MaaslinParameters)) try(dir.create(outdir_MaaslinParameters, recursive= TRUE), silent= TRUE)


# ----------------------------------------------------------------------------------------------------------------------------------------
## ANALYSIS
# ----------------------------------------------------------------------------------------------------------------------------------------


  ### Continuous ####
  # ==============================================================================================================
  Maaslin2(input_data = OTU,
           input_metadata = metadata,

           fixed_effects = c("BSI"),
           # fixed_effects = c("FEV1_percent_predicted"),
           # fixed_effects = c("MRC_score"),

           ##----- COVARIATES CORRECTION
           random_effects = c("Country"),

           output = outdir_MaaslinParameters,
           min_abundance = 0,
           min_prevalence = 0, 
           analysis_method = "CPLM", 
           normalization = "NONE",
           transform = "NONE", 
           correction = "BH",
           standardize = TRUE,
           cores = 5)

  
  
  ### Categorical ####
  # ==============================================================================================================
  
  ##----- FrequentExacerbator
  Maaslin2(input_data = OTU,
           input_metadata = metadata,
           fixed_effects = "FrequentExacerbator",
           reference = c("FrequentExacerbator,Non_FE"),

           ##-----COVARIATES
           random_effects = c("Country"),

           output = outdir_MaaslinParameters,
           min_abundance = 0, 
           min_prevalence = 0, 
           analysis_method = "CPLM",
           normalization = "NONE", 
           transform = "NONE", 
           correction = "BH",
           standardize = TRUE,
           cores = 5)
  

  ##----- Hospitalised_for_bronchiectasis_YesNo
  Maaslin2(input_data = OTU,
           input_metadata = metadata,
           fixed_effects = "Hospitalised_for_bronchiectasis_YesNo",
           reference = c("Hospitalised_for_bronchiectasis_YesNo,No"),

           ##-----COVARIATES
           random_effects = c("Country"),

           output = outdir_MaaslinParameters,
           min_abundance = 0, 
           min_prevalence = 0, 
           analysis_method = "CPLM",
           normalization = "NONE", 
           transform = "NONE",
           correction = "BH",
           standardize = TRUE,
           cores = 5)

  

