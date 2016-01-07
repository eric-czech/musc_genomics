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

# options(enable.cache=T)
# options(enable.cache=F)

#---------------------------#
# Load cBioPortal CCLE Data #
#---------------------------#

# Function used to put CCLE tumor ids into columns and prefix all
# previous columns with a short string indicating the data type
custom.tumorids <- c('TT_OESOPHAGUS'='TT1', 'TT_THYROID'='TT2')
prep <- function(d, prefix) d %>% 
  # Make tumor id first column
  select(tumor_id, everything()) %>%
  # Add column prefix indicating origin of dataset (e.g. 'cn' for Copy Number)
  setNames(., c('tumor_id', paste(prefix, names(.), sep='.')[-1])) %>% 
  # Remove "X" prefixed to all tumor IDs by CGDS API
  mutate(tumor_id=str_replace(tumor_id, '^X', '')) %>% 
  # Match any known conflicts in tumor IDs to custom IDs
  mutate(tumor_id=ifelse(tumor_id %in% names(custom.tumorids), custom.tumorids[tumor_id], tumor_id)) %>%
  # Restrict tumor IDs to only substring before first underscore
  mutate(tumor_id=str_replace(tumor_id, '_.*', ''))

# Get list of unique gene symbols to collect data for
gene.symbols <- GetHUGOGeneNames()

# Load copy number, gene expression, and mutation data
biop.cn <- GetBioPortalData(gene.symbols, GEN_PROF_COPY_NUMBER) %>% prep('cn')
biop.ge <- GetBioPortalData(gene.symbols, GEN_PROF_GENE_EXPRESSION) %>% prep('ge')
biop.mu <- GetBioPortalData(gene.symbols, GEN_PROF_MUTATION) %>% prep('mu')

biop.ct <- sapply(list(biop.cn, biop.ge, biop.mu), nrow)
if (length(unique(biop.ct)) != 1)
  stop('BioPortal datasets have differing numbers of rows')

# Merge all datasets above into one data frame
biop.data <- biop.cn %>% 
  inner_join(biop.ge, by='tumor_id') %>%
  inner_join(biop.mu, by='tumor_id')

# Verify that all tumor IDs were found in each dataset
if (nrow(biop.data) != nrow(biop.cn))
  stop('Some records lost after join on tumor ID (review and retry)')

#----------------#
# Load CTD2 Data #
#----------------#

ctd2.data <- GetCTD2Data()

data <- biop.data %>% left_join(ctd2.data, by='tumor_id')

#------------------#
# Load COSMIC Data #
#------------------#
# TBD



