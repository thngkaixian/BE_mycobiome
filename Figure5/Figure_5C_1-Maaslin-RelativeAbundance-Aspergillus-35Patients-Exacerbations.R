## Maaslin analysis for ITS with homes sequenced 
library(phyloseq)
library(Maaslin2)


# --------------------------------------------------------------------------------------------------------------
## Directories 
# --------------------------------------------------------------------------------------------------------------

taxLevel_Microbiome <- "Species"

wkdir <- "/Users/thngkaixian/PhD_RESEARCH/RESEARCH_Projects"
indir_ps <- file.path(wkdir, "MYCOBIOME_AsiaBRIDGE", "SequencingData", "FINAL_ANALYSIS")

## Out directories
outdir <- file.path(indir_ps, "Analysis", "Figure_5C", "Maaslin-Aspergillus-Species-vs-ClinicalOutcome", taxLevel_Microbiome)
if (!dir.exists(outdir)) try(dir.create(outdir, recursive= TRUE), silent= TRUE)

  
# --------------------------------------------------------------------------------------------------------------
## Import data: HHP data
# --------------------------------------------------------------------------------------------------------------

SampleType <- "Outdoor"
# SampleType <- "Indoor"
# SampleType <- "SwabS"
  
  
  ## Import microbiome ps object 
  ps_microbiome_HHP <- readRDS(file.path(indir_ps, paste0("Final_ps_", taxLevel_Microbiome, "_HHP_35patients.RData")))
  ps_microbiome_HHP <- prune_samples(sample_names(ps_microbiome_HHP)[sample_data(ps_microbiome_HHP)$SampleType==SampleType], ps_microbiome_HHP)
    nsamples(ps_microbiome_HHP)
    sample_names(ps_microbiome_HHP)
    table(sample_data(ps_microbiome_HHP)$LKC_ID)
    
    ## Summary of samples analysed
    # ----------------------------------------------------------------------------------------------------
    dfSummary_SamplesAnalysed <- table(sample_data(ps_microbiome_HHP)[, c("Cohort", "Home.Sampling.Type")]) ; dfSummary_SamplesAnalysed

    
    ## SET LKC IDs as sample name 
    # ----------------------------------------------------------------------------------------------------
    sample_names(ps_microbiome_HHP) <- sample_data(ps_microbiome_HHP)$LKC_ID
    # write.csv(sample_names(ps_microbiome_HHP), file = file.path(outdir, paste0("Microbiome_IDs_from_HHP", ".csv")))
    nsamples(ps_microbiome_HHP)
    sample_names(ps_microbiome_HHP)
    
    
# --------------------------------------------------------------------------------------------------------------=
## metadata wrangling
# --------------------------------------------------------------------------------------------------------------=
    
## Import additional clinical data
Additional_metadata  <- read.csv(file.path(wkdir, "HomeHospitals", "metadata", "HomeHospital_Patient_Metadata-AdditionalClinicalMetadata-Bronchiectasis.csv"), row.names = 1)
    
  ## Extract metadata from ps
  metadata <- data.frame(sample_data(ps_microbiome_HHP))

  ## Merge Additional metadata
  # all(metadata$LKC_ID %in% Additional_metadata$Patient.ID)
  # all(Additional_metadata$Patient.ID %in% metadata$LKC_ID)
  metadata <- merge(metadata, Additional_metadata, by.x = "LKC_ID", by.y="Patient.ID", all.x = TRUE, all.y = FALSE)

  
  ## Wrangle metadata 
  # --------------------------------------------------------------------------------------------------------------=
  metadata$Exacerbation_yes_no <- sapply(metadata$Number_of_exacerbations_in_the_previous_year, function(x){if(x==0) "Exac_No" else if (x>0) "Exac_Yes"})
  
    
  ## factor environemnt variables 
  # --------------------------------------------------------------------------------------------------------------=
  
  #### Windows ####
  metadata$Window_open_room_hrs <- factor(metadata$Window_open_room_hrs, levels = c("0-6h", "6-12h", "12-18h", "18-24h"))
  metadata$Window_open_room_hrs_1234 <- sapply(metadata$Window_open_room_hrs, function(x){ if(x=="0-6h") 1 else if (x=="6-12h") 2 else if (x=="12-18h") 3 else if (x=="18-24h") 4 })
  
  metadata$Window_open_room <- factor(metadata$Window_open_room, levels = c("Never", "< 3 times a week", "Daily"), labels = c("Never", "1_2_days", "Daily"))
  metadata$Window_open_room_123 <- sapply(metadata$Window_open_room, function(x){ if(x=="Never") 1 else if (x=="1_2_days") 2 else if (x=="Daily") 3 })
  
  # table(metadata$Window_open_room, metadata$Window_open_room_hrs)
  # table(metadata$Window_open_room)
  # table(metadata$Window_open_room_hrs)
  
  
      metadata$Window_open_room_NeverYes <-sapply(metadata$Window_open_room, function(x){ if(x=="Never") "Never" else if(x=="1_2_days" | x=="Daily") "Yes" })
      metadata$Window_open_room_NeverYes <- factor(metadata$Window_open_room_NeverYes, levels = c("Never", "Yes"))
      metadata$Window_open_room_NeverYes_12 <- sapply(metadata$Window_open_room_NeverYes, function(x){ if(x=="Never") 1 else if (x=="Yes") 2 })
      
      metadata$Window_open_room_Daily_012days <-sapply(metadata$Window_open_room, function(x){ if(x=="Never" | x=="1_2_days") "0_1_2_days" else if(x=="Daily") "Daily" })
      metadata$Window_open_room_Daily_012days <- factor(metadata$Window_open_room_Daily_012days, levels = c("0_1_2_days", "Daily"))
      table(metadata$Window_open_room_Daily_012days)
      

      metadata$Window_open_room_hrs_12h <- sapply(metadata$Window_open_room_hrs, function(x){ if(x=="0-6h" | x=="6-12h") "Less_than_12h" else if(x=="12-18h" | x=="18-24h") "More_than_12h"})
      metadata$Window_open_room_hrs_12h <- factor(metadata$Window_open_room_hrs_12h, levels = c("Less_than_12h", "More_than_12h"))
      table(metadata$Window_open_room_hrs_12h)
      metadata$Window_open_room_hrs_12h_12 <- sapply(metadata$Window_open_room_hrs_12h, function(x){ if(x=="Less_than_12h") 1 else if (x=="More_than_12h") 2 })
      
      metadata$Window_open_room_hrs_6_18_24 <- sapply(metadata$Window_open_room_hrs, function(x){ if(x=="0-6h") "0-6h" else if(x=="6-12h" | x=="12-18h") "6-18h" else if(x=="18-24h") "18-24h"})
      metadata$Window_open_room_hrs_6_18_24 <- factor(metadata$Window_open_room_hrs_6_18_24, levels = c("0-6h", "6-18h", "18-24h"))
      table(metadata$Window_open_room_hrs_6_18_24)
      
      metadata$Window_12h_Daily <- paste0(metadata$Window_open_room_hrs_12h, "_", metadata$Window_open_room_Daily_012days)
      metadata$Window_12h_Daily <- factor(metadata$Window_12h_Daily, levels = c("Less_than_12h_0_1_2_days", "Less_than_12h_Daily", "More_than_12h_0_1_2_days", "More_than_12h_Daily"))
      table(metadata$Window_12h_Daily)
      
      
  #### Aircon_use_room ####
  metadata$Aircon_use_room <- factor(metadata$Aircon_use_room, levels = c("Never", "0-3h", "3-12h", "12-24h"))
  table(metadata$Aircon_use_room)
  metadata$Aircon_use_room_1234 <- sapply(metadata$Aircon_use_room, function(x){ if(x=="Never") 1 else if (x=="0-3h") 2 else if (x=="3-12h") 3 else if(x=="12-24h") 4 })
      
      metadata$Aircon_use_room_YesNo <- sapply(metadata$Aircon_use_room, function(x){ if(x=="Never") "Never" else if(x=="0-3h" | x=="3-12h" | x=="12-24h") "Yes"})
      metadata$Aircon_use_room_YesNo <- factor(metadata$Aircon_use_room_YesNo, levels = c("Never", "Yes"))
      table(metadata$Aircon_use_room_YesNo)
      metadata$Aircon_use_room_YesNo_12 <- sapply(metadata$Aircon_use_room_YesNo, function(x){ if (x=="Never") 1 else if (x=="Yes") 2})
      
      metadata$Aircon_use_room_Never_3_24 <- sapply(metadata$Aircon_use_room, function(x){ if(x=="Never") "Never" else if (x=="0-3h") "0-3h" else if(x=="3-12h" | x=="12-24h") "more_than_3h"})
      metadata$Aircon_use_room_Never_3_24 <- factor(metadata$Aircon_use_room_Never_3_24, levels = c("Never", "0-3h", "more_than_3h"))
  
      
  ##### Return into ps object
  # --------------------------------------------------------------------------------------------------------------
  row.names(metadata) <- metadata$Sequencing_ID
      
  sample_data(ps_microbiome_HHP) <- metadata
      
      
# ----------------------------------------------------------------------------------------------------
## SELECT PS TO BE USED 
# ----------------------------------------------------------------------------------------------------
    
PS <- ps_microbiome_HHP

  PS <- prune_taxa(taxa_names(PS)[apply(otu_table(PS), 2, sum) != 0], PS) # remove empty taxa
  PS <- prune_samples(names(which(apply(otu_table(PS), 1, sum) != 0)), PS) # remove patients without taxa
  PS.prop <- transform_sample_counts(PS, function(otu) {if (sum(otu)==0) otu else 100*otu/sum(otu)})

  
  ## Subset specific taxa for analysis 
  PS.prop_Aspergillus <- subset_taxa(PS.prop, Genus=="Aspergillus")
  # taxa_names(PS.prop_Aspergillus)

  ps_subgroup <- PS.prop_Aspergillus

  
# ----------------------------------------------------------------------------------------------------------------------------------------
## Preparing metadata for input to maaslin
# ----------------------------------------------------------------------------------------------------------------------------------------  
  
metadata <- data.frame(sample_data(ps_subgroup))
OTU <- data.frame(otu_table(ps_subgroup))


# ==============================================================================================================
## run masslin (PARAMETERS 2) - CPLM
# ==============================================================================================================
  
  ## update out-directory with parameters
  outdir_MaaslinParameters <- file.path(outdir, 
                                        # outdir_Filter,
                                        paste0("Model", "CPLM", "_",
                                               "Normalisation", "NONE", "_",
                                               "Transformation", "NONE"))
  if (!dir.exists(outdir_MaaslinParameters)) try(dir.create(outdir_MaaslinParameters, recursive= TRUE), silent= TRUE)
  
  
# ----------------------------------------------------------------------------------------------------------------------------------------
## Exacerbations
# ----------------------------------------------------------------------------------------------------------------------------------------
  
  unique(metadata$SampleType)
  metadata$Home.Sampling.Type
  
  ### Continuous ####
  # ==============================================================================================================
  Maaslin2(input_data = OTU,
           input_metadata = metadata,

           ## Outdoor --- Exacerbation
           # ----------------------------------------------------------------------------------------------------
           fixed_effects = "Number_of_exacerbations_in_the_previous_year",

           ## Indoor --- Exacerbation
           # ----------------------------------------------------------------------------------------------------
           # fixed_effects = "Number_of_exacerbations_in_the_previous_year",
           # random_effects = "Aircon_use_room,Window_open_room_Daily_012days,Window_open_room_hrs_12h,Use_of_fan_in_room",

           ## Surfaces --- Exacerbation
           # ----------------------------------------------------------------------------------------------------
           # fixed_effects = "Number_of_exacerbations_in_the_previous_year",
           # random_effects = "Home.Sampling.Type", #SwabS correction

           
           output = outdir_MaaslinParameters,
           min_abundance = 0,
           min_prevalence = 0,
           analysis_method = "CPLM",
           normalization = "NONE",
           transform = "NONE",
           correction = "BH",
           standardize = TRUE,
           cores = 5)
  





  
  
