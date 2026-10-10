library(phyloseq)


# --------------------------------------------------------------------------------------------------------------
## Directories 
# --------------------------------------------------------------------------------------------------------------

taxLevel_Microbiome <- "Species"

wkdir <- "/Users/thngkaixian/PhD_RESEARCH/RESEARCH_Projects"
indir_ps <- file.path(wkdir, "MYCOBIOME_AsiaBRIDGE", "SequencingData", "FINAL_ANALYSIS")

## Out directories
outdir <- file.path(indir_ps, "Analysis", "Figure_6CD", "Indoor-AspergillusFumigatus-vs-WindowUse", taxLevel_Microbiome)
if (!dir.exists(outdir)) try(dir.create(outdir, recursive= TRUE), silent= TRUE)

  
# --------------------------------------------------------------------------------------------------------------
## Import data: HHP data
# --------------------------------------------------------------------------------------------------------------

SampleType <- "Indoor"

  ## Import microbiome ps object 
  ps_microbiome_HHP <- readRDS(file.path(indir_ps, paste0("Final_ps_", taxLevel_Microbiome, "_HHP_35patients.RData")))
  ps_microbiome_HHP <- prune_samples(sample_names(ps_microbiome_HHP)[sample_data(ps_microbiome_HHP)$SampleType==SampleType], ps_microbiome_HHP)
    nsamples(ps_microbiome_HHP)  
    sample_names(ps_microbiome_HHP)
    table(sample_data(ps_microbiome_HHP)$LKC_ID)
    
    
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
  
  
  ## factor environemnt variables 
  # --------------------------------------------------------------------------------------------------------------=
  
  #### Windows ####
  metadata$Window_open_room_hrs <- factor(metadata$Window_open_room_hrs, levels = c("0-6h", "6-12h", "12-18h", "18-24h"))
  metadata$Window_open_room_hrs_1234 <- sapply(metadata$Window_open_room_hrs, function(x){ if(x=="0-6h") 1 else if (x=="6-12h") 2 else if (x=="12-18h") 3 else if (x=="18-24h") 4 })
  
  metadata$Window_open_room <- factor(metadata$Window_open_room, levels = c("Never", "< 3 times a week", "Daily"), labels = c("Never", "1_2_days", "Daily"))
  metadata$Window_open_room_123 <- sapply(metadata$Window_open_room, function(x){ if(x=="Never") 1 else if (x=="1_2_days") 2 else if (x=="Daily") 3 })

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
      
      
      #### Never, <12h, >12h windows ####
      metadata[, c("Window_open_room", "Window_open_room_hrs_12h", "Window_open_room_hrs")]
      table(metadata$Window_open_room_Daily_012days, metadata$Window_open_room_hrs_12h)
      table(metadata$Window_open_room, metadata$Window_open_room_hrs_12h)
      
      metadata$Window_open_room_FreqDuration <- paste0(metadata$Window_open_room, "_", metadata$Window_open_room_hrs_12h)
      table(metadata$Window_open_room_FreqDuration, exclude = NULL)
      
      metadata$Window_open_room_Never_MtLt12h <-sapply(metadata$Window_open_room_FreqDuration, function(x){ if(x=="Never_Less_than_12h") "Never" else if(x=="1_2_days_Less_than_12h" | x=="Daily_Less_than_12h") "less_12h" else if (x=="1_2_days_More_than_12h" | x=="Daily_More_than_12h") "more_12h" })
      table(metadata$Window_open_room_Never_MtLt12h, exclude = NULL)
      
      
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
      
  ## RE-update ps object 
  sample_data(ps_microbiome_HHP) <- metadata
      
      
# ----------------------------------------------------------------------------------------------------
## SELECT PS TO BE USED 
# ----------------------------------------------------------------------------------------------------
    
PS <- ps_microbiome_HHP

  PS <- prune_taxa(taxa_names(PS)[apply(otu_table(PS), 2, sum) != 0], PS) # remove empty taxa
  PS <- prune_samples(names(which(apply(otu_table(PS), 1, sum) != 0)), PS) # remove patients without taxa
  PS.prop <- transform_sample_counts(PS, function(otu) {if (sum(otu)==0) otu else 100*otu/sum(otu)})


# ----------------------------------------------------------------------------------------------------------------------------------------
## Extract metadata from ps object --- create dataframe for analysis
# ----------------------------------------------------------------------------------------------------------------------------------------  
  
  metadata <- data.frame(sample_data(PS.prop))
  metadata <- metadata[, -c(53:304)] 
  OTU <- data.frame(otu_table(PS.prop))

  ## merge 
  metadataOTU <- merge(metadata, OTU, by.x = "row.names", by.y = "row.names", all.x = TRUE, all.y = TRUE)
  row.names(metadataOTU) <- metadataOTU$Row.names ; metadataOTU$Row.names <- NULL
 
 
# ----------------------------------------------------------------------------------------------------------------------------------------  
## ANALYSIS 
# ----------------------------------------------------------------------------------------------------------------------------------------  
library(dunn.test)

metadataOTU$Window_open_room_Never_MtLt12h <- factor(metadataOTU$Window_open_room_Never_MtLt12h, levels = c("Never", "less_12h", "more_12h"))
metadataOTU$Window_open_room_hrs_12h

  #### Window_open_room_NeverYes ####
  # ----------------------------------------------------------------------------------------------------------------------------------------  
  Wil_test <- wilcox.test(Aspergillus.fumigatus ~ Window_open_room_NeverYes, data = metadataOTU)
  pval <- Wil_test$p.value
  
  p3 <- ggplot(metadataOTU, aes(x=Window_open_room_NeverYes, y=Aspergillus.fumigatus)) +
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(aes(colour=Window_open_room_Never_MtLt12h), width = 0.1, height = 0.1) + #Window_open_room_Never_MtLt12h | Window_open_room_hrs
    ggtitle("Indoor", subtitle = paste0("p=", round(pval, 3))) +
    scale_y_log10()
  p3
  p3_out_filename <- file.path(outdir, paste0("Boxplot-", "IndoorAF-", "Window_open_room_NeverYes", ".png"))
  ggsave(plot = p3, filename = p3_out_filename, width=14, height = 10, units="cm")
  
  p3_Final <- p3 +
    labs(title = NULL, subtitle = NULL,
         y="Indoor *A. fumigatus*<br>relative abundance (log<sub>10</sub>)"#,
         ) + 
    scale_x_discrete(name="Window use frequency", labels=c("Never"="Never\nopen", "Yes"="At least 1 day\nper week")) + 
    scale_color_manual(name="Window use duration",
                       labels=c("Never"="0 hour", "less_12h"="<12 hours per day", "more_12h"=">12 hours per day"),
                         values=c("Never"="mediumblue", "less_12h"="violetred", "more_12h"="darkgreen")) +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.title.y = ggtext::element_markdown(),
          axis.line = element_line(colour="black"))
  p3_Final
  p3_Final_out_filename <- file.path(outdir, paste0("Boxplot-", "IndoorAF-", "Window_open_room_NeverYes", "-FINAL", ".png"))
  ggsave(plot = p3_Final, filename = p3_Final_out_filename, width=12, height = 8, units="cm")
  
  
  
  
  
  
  
