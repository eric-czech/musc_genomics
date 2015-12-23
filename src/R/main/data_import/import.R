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
lib('cgdsr')
lib('foreach')
lib('dplyr')

#---------------------------#
# Load cBioPortal CCLE Data #
#---------------------------#

# Function used to put CCLE tumor ids into columns and prefix all
# previous columns with a short string indicating the data type
custom.tumorids <- c('TT_OESOPHAGUS'='TT1', 'TT_THYROID'='TT2')
prep <- function(d, prefix) d %>% 
  # Add column prefix
  setNames(., paste(prefix, names(.), sep='.')) %>% 
  # Put row name in column named 'tumor_id'
  add_rownames(var='tumor_id') %>%
  # Match any known conflicts in tumor IDs to custom IDs
  mutate(tumor_id=ifelse(tumor_id %in% names(custom.tumorids), custom.tumorids[tumor_id], tumor_id)) %>%
  # Restrict tumor IDs to only substring before first underscore
  mutate(tumor_id=str_replace(tumor_id, '_.*', ''))

# Get list of unique gene symbols to collect data for
gene.symbols <- GetHugoGeneNames()

# Load copy number, gene expression, and mutation data
biop.cn <- GetBioPortalData(gene.symbols, GEN_PROF_COPY_NUMBER) %>% prep('cn')
biop.ge <- GetBioPortalData(gene.symbols, GEN_PROF_GENE_EXPRESSION) %>% prep('ge')
biop.mu <- GetBioPortalData(gene.symbols, GEN_PROF_MUTATION) %>% prep('mu')

# Merge all datasets above into one data frame
biop.data <- biop.cn %>% 
  inner_join(biop.ge, by='tumor_id') %>%
  inner_join(biop.mu, by='tumor_id')

#------------------#
# Load COSMIC Data #
#------------------#
# TBD

#----------------#
# Load CTD2 Data #
#----------------#
# TBD
