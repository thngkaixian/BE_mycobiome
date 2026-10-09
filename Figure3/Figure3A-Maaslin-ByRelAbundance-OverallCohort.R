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

outdir <- file.path(wkdir, "FINAL_ANALYSIS", "Analysis", "Figure_3A")
if (!dir.exists(outdir)) try(dir.create(outdir, recursive= TRUE), silent= TRUE)


# --------------------------------------------------------------------------------------------------------------
## Read ps objects 
# --------------------------------------------------------------------------------------------------------------

ps <- readRDS(file.path(indir_ps, paste0("Final_psITS_", taxLevel, ".RData")))
ps <- prune_taxa(taxa_names(ps)[!grepl('unidentified',taxa_names(ps))], ps)

  ## Filtering
  # --------------------------------------------------------------------------------------------------------------
  # Filter<-"Filtered"
  
  RelAbundThreshold <- 0.01
  sampleThreshold <- 0.05
  
  ps.prop <- transform_sample_counts(ps, function(otu) {if (sum(otu)==0) otu else otu/sum(otu)})
  tax <- taxa_names(ps.prop)[apply(otu_table(ps.prop), 2, function(x) {sum(x >= RelAbundThreshold) >= sampleThreshold*nsamples(ps)})]
  ps <- prune_taxa(tax, ps)
  
  ## Relative abundance 
  ps.prop <- transform_sample_counts(ps, function(otu) {if (sum(otu)==0) otu else 100*otu/sum(otu)})
  
  
# ----------------------------------------------------------------------------------------------------------------------------------------
## Preparing metadata for input to maaslin
# ----------------------------------------------------------------------------------------------------------------------------------------  
  
metadata <- data.frame(sample_data(ps.prop))
OTU <- data.frame(otu_table(ps.prop))


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
           # fixed_effects = "BSI",
           fixed_effects = "MRC_score",
           # fixed_effects = "FEV1_percent_predicted",
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
           output = outdir_MaaslinParameters,
           min_abundance = 0, 
           min_prevalence = 0, 
           analysis_method = "CPLM", 
           normalization = "NONE",
           transform = "NONE", 
           correction = "BH", 
           standardize = TRUE,
           cores = 5)
  
  
  Maaslin2(input_data = OTU,
           input_metadata = metadata,
           fixed_effects = "Hospitalised_for_bronchiectasis_YesNo",
           reference = c("Hospitalised_for_bronchiectasis_YesNo,No"),
           output = outdir_MaaslinParameters,
           min_abundance = 0, 
           min_prevalence = 0, 
           analysis_method = "CPLM", 
           normalization = "NONE", 
           transform = "NONE",
           correction = "BH",
           standardize = TRUE,
           cores = 5)









