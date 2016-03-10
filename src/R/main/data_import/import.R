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
    left_join(d.ctd2, by='tumor_id') %>% 
    left_join(d.cosmic, by='tumor_id') %>%
    select(tumor_id, origin, ic_50, auc, everything())
}

# d.biop <- GetBioPortalData()
# d <- GetRawData()


