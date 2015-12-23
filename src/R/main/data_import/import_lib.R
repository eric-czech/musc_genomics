#'-----------------------------------------------------------------------------
#' Utilities for data importation scripts
#'
#' This module contains several functions for fetching or processing raw data
#' from HUGO, cBioPortal, and COSMIC databases
#' 
#' @author eczech
#'-----------------------------------------------------------------------------
source('utils.R')
source('cache.R')

lib('dplyr')
lib('stringr')
lib('RCurl')
lib('cgdsr')


#------------------------------#
# HUGO Constants and Functions #
#------------------------------#

HUGO_GENE_NAMES_URL <- 'https://s3-us-west-2.amazonaws.com/eric.a.czech/Public/musc_genomics/hugo_genes.csv'

GetHUGOGeneNames <- function(){
  #' Returns a vector of standardized gene symbols from HUGO (see http://www.genenames.org/).
  #' 
  #' Note that all genes are returned except for this with Approved.Name of 
  #' 'entry withdrawn' or 'symbol withdrawn'
  
  # Load raw HUGO data frame, remove withdrawn genes, and return vector of unique symbols
  loader <- function() read.csv(textConnection(getURL(HUGO_GENE_NAMES_URL)), sep=',', stringsAsFactors=F)
  Fetch('hugo_gene_names', loader) %>%
    filter(!str_detect(tolower(Approved.Name), 'symbol withdrawn|entry withdrawn')) %>%
    .$Approved.Symbol %>% unique
}

#------------------------------------#
# cBioPortal Constants and Functions #
#------------------------------------#

GetCGDS <- function() CGDS("http://www.cbioportal.org/public-portal/")
GEN_PROF_COPY_NUMBER <- 'cellline_ccle_broad_log2CNA'
GEN_PROF_GENE_EXPRESSION <- 'cellline_ccle_broad_mrna_median_Zscores'
GEN_PROF_MUTATION <- 'cellline_ccle_broad_mutations'

GetBioPortalData <- function(gene.symbols, genetic.profile, chunk.size=50){
  #' Returns cBioPortal data for the given genetic profile
  #' 
  #' Note that the list of genes to get data for will be split into chunks
  #' and data for each chunk will be collected in separate calls
  
  cgds <- GetCGDS()  
  loader <- function(){
    # Split the gene symbol list into chunks of size no greater than #chunk.size
    gene.partitions <- split(gene.symbols, ceiling(seq_along(gene.symbols)/CHUNK_SIZE))
    
    # Fetch CCLE data for the given #genetic.profile, for each chunk
    temp.chunks <- gene.partitions[1:5] # TODO: Remove this limit later
    data <- foreach(genes=temp.chunks)%do%{
      getProfileData(cgds, genes, genetic.profile, "cellline_ccle_broad_all") 
    }
    
    # Verify that row names from each chunk are equivalent.  This is necessary to
    # make sure that non only does the data for each chunk have the same dimensions
    # but also that the row names (tumor IDs) are sorted in the same order
    rows.equal <- Reduce(function(d1, d2) all(row.names(d1) == row.names(d2)), data)
    if (!rows.equal) stop('Data for some chunks of genes returned differing row names')
    
    # Combine all data chunks into a single data frame
    foreach(d=data, .combine=cbind) %do% d
  }
  # Lazy-load these results (they're expensive to compute), saving them
  # on disk or loading from disk if previously created
  FetchFromDisk(genetic.profile, loader)
}

