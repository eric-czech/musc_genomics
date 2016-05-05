#'-----------------------------------------------------------------------------
#' CGDS Exploration Script
#'
#' This module contains several function calls for verifying results and formats
#' for data from the CGDS API
#' 
#' @author eczech
#'-----------------------------------------------------------------------------

source('data_import/import_lib.R')

cgds <- GetCGDS()

gp <- getGeneticProfiles(cgds, 'cellline_ccle_broad')

cl <- getCaseLists(cgds, 'cellline_ccle_broad')

pd <- getProfileData(cgds, c('KRAS'), 'cellline_ccle_broad_mrna_median_Zscores', "cellline_ccle_broad_all") 