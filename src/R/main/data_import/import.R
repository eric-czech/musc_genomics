#'-----------------------------------------------------------------------------
#' Genomics data importation script
#'
#' This script will collect data from the cBioPortal database and merge this
#' information with sensitivity data in the COSMIC and CTD2 databases.
#' 
#' @author eczech
#'-----------------------------------------------------------------------------
source('utils.R')
source('data_import/import_lib.R')
lib('dplyr')


GetRawData <- function(){
  d.biop <- GetBioPortalData() # Fetch cell line genomic predictors
  d.ctd2 <- GetCTD2V2Data()    # Fetch AUC response data
  d.cosmic <- GetCOSMICData()  # Fetch IC 50 response data
  
  # Join everything together and return merged result
  d <- d.biop %>% 
    # Add CTD and COSMIC data
    left_join(d.ctd2, by='tumor_id') %>% 
    left_join(d.cosmic, by='tumor_id') %>%
    # Replace NA tissue values which at TOW, were only present for COSMIC
    mutate(tissue=ifelse(is.na(tissue), 'UNKNOWN', tissue)) %>% 
    # Rearrange columns in result
    select(tumor_id, origin, tissue, ic_50, auc, everything()) 
}

# d.biop <- GetBioPortalData()
# d <- GetRawData()


