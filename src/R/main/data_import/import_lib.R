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
lib('foreach')
lib('iterators')
lib('reshape2')
select <- dplyr::select


#--------------------------------------#
##### HUGO Constants and Functions #####
#--------------------------------------#

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

#---------------------------------------------#
##### cBioPortal Constants and Functions ######
#---------------------------------------------#

GetCGDS <- function() CGDS("http://www.cbioportal.org/public-portal/")
GEN_PROF_COPY_NUMBER <- 'cellline_ccle_broad_log2CNA'
GEN_PROF_GENE_EXPRESSION <- 'cellline_ccle_broad_mrna_median_Zscores'
GEN_PROF_MUTATION <- 'cellline_ccle_broad_mutations'


GetGeneticProfileData <- function(gene.symbols, genetic.profile, chunk.size=50){
  #' Returns cBioPortal data for the given genetic profile
  #' 
  #' Note that the list of genes to get data for will be split into chunks
  #' and data for each chunk will be collected in separate calls
  
  cgds <- GetCGDS()  
  loader <- function(){
    # Split the gene symbol list into chunks of size no greater than #chunk.size
    gene.partitions <- split(gene.symbols, ceiling(seq_along(gene.symbols)/chunk.size))
    
    # Utility function used to convert factor data into character data
    # with one minor provision for when a factor has a 'NaN' level (not NA); 
    # this seems to happen only for cell lines with at least one present mutation
    to.char <- function(x) {
      if (is.factor(x)) {
        r <- as.character(x) 
        ifelse(tolower(r) == 'nan', NA, r)
      } else x
    }
    
    # Fetch CCLE data for the given #genetic.profile, for each chunk
    #temp.chunks <- gene.partitions[1:Inf] # TODO: Remove this limit later
    # Note: The following data comes from the cancer study named 'cellline_ccle_broad'
    data <- foreach(genes=gene.partitions)%do%{
      getProfileData(cgds, genes, genetic.profile, "cellline_ccle_broad_all") %>%
        add_rownames(var='tumor_id') %>% mutate_each(funs(to.char)) 
    }
    
    
    # Verify that tumor IDs from each chunk are equivalent.  This is necessary to
    # make sure that non only does the data for each chunk have the same dimensions
    # but also that the tumor IDs are sorted in the same order
    ids.equal <- foreach(i=1:(length(data)-1), .combine=c) %do% {
      all.equal(data[[i]]$tumor_id, data[[i+1]]$tumor_id)
    } %>% all
    if (!ids.equal) stop('Data for some chunks of genes returned differing tumor ids')
    
    # Combine all data chunks into a single data frame and remove NA-only columns
    tumor.ids <- data[[1]]$tumor_id
    (foreach(d=data, .combine=cbind) %do% {d %>% select(-tumor_id)}) %>% 
      Filter(function(x)!all(is.na(x)), .) %>%
      mutate(tumor_id=tumor.ids)
  }
  # Lazy-load these results (they're expensive to compute), saving them
  # on disk or loading from disk if previously created
  RAW_CACHE$load(genetic.profile, loader)
  #FetchFromDisk(genetic.profile, loader) 
}

GetBioPortalData <- function(){
  #' Returns cBioPortal data for all desired genetic profiles.  
  #' 
  #' In other words, this function returns a unified data frame containing
  #' all gene expression, copy number, and mutation data where each of those
  #' have columns corresponding to different genes.  In the resulting data frame,
  #' those columns containing gene names are prefixed by a short string describing
  #' which kind of data they contain (e.g. 'cn.' for copy number)
  #' 
  #' Note that this function may take a long time to run.  Previous run times:
  #' Trial 1 - 1hr 19mins (for all genes -- 38,649 features in total)
  
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
    # Extract tumor origin from tumor_id (e.g. 1321N1_CENTRAL_NERVOUS_SYSTEM -> CENTRAL_NERVOUS_SYSTEM)
    mutate(origin=str_extract(tumor_id, '(?<=_).*$')) %>%
    mutate(origin=ifelse(is.na(origin), 'Unknown', origin)) %>%
    # Restrict tumor IDs to only substring before first underscore
    mutate(tumor_id=str_replace(tumor_id, '_.*', ''))
  
  # Get list of unique gene symbols to collect data for
  gene.symbols <- GetHUGOGeneNames()
  
  # Load copy number, gene expression, and mutation data
  biop.cn <- GetGeneticProfileData(gene.symbols, GEN_PROF_COPY_NUMBER) %>% prep('cn')
  biop.ge <- GetGeneticProfileData(gene.symbols, GEN_PROF_GENE_EXPRESSION) %>% prep('ge')
  biop.mu <- GetGeneticProfileData(gene.symbols, GEN_PROF_MUTATION) %>% prep('mu')
  
  biop.ct <- sapply(list(biop.cn, biop.ge, biop.mu), nrow)
  if (length(unique(biop.ct)) != 1)
    stop('BioPortal datasets have differing numbers of rows')
  
  # Merge all datasets above into one data frame
  biop.data <- biop.cn %>% 
    inner_join(biop.ge, by=c('tumor_id', 'origin')) %>%
    inner_join(biop.mu, by=c('tumor_id', 'origin'))
  
  # Verify that all tumor IDs were found in each dataset
  if (nrow(biop.data) != nrow(biop.cn))
    stop('Some records lost after join on tumor ID (review and retry)')
  
  biop.data
}


#--------------------------------------#
##### CTD2 Constants and Functions #####
#--------------------------------------#

# Links found by browsing here: https://ctd2.nci.nih.gov/dataPortal/
CTD2_V1_URL <- 'ftp://caftpd.nci.nih.gov/pub/dcc_ctd2/Broad/CTRPv1.0_2013_pub_Cell_154_1151/CTRPv1.0_2013_pub_Cell_154_1151.zip'  
CTD2_V2_URL <- 'ftp://caftpd.nci.nih.gov/pub/dcc_ctd2/Broad/CTRPv2.0_2015_ctd2_ExpandedDataset/CTRPv2.0_2015_ctd2_ExpandedDataset.zip'  

GetCTD2V2Data <- function(){
  loader <- function(){
    file.path <- RAW_CACHE$download('ctd2_auc_v2.zip', CTD2_V2_URL, mode="wb")
    
    # Load raw AUC data for experiments
    d.auc <- read.csv(unz(file.path, 'v20.data.curves_post_qc.txt'), sep='\t', stringsAsFactors=F) %>%
      select(area_under_curve, master_cpd_id, experiment_id)
    
    # Load compound/drug meta data
    d.cmpd <- read.csv(unz(file.path, 'v20.meta.per_compound.txt'), sep='\t', stringsAsFactors=F) %>%
      select(cpd_name, master_cpd_id)
    
    # Load experiment meta data (note that some of these records are duplicated in full -- no idea why)
    # * Example experiment with duplicates = 85
    d.exp <- read.csv(unz(file.path, 'v20.meta.per_experiment.txt'), sep='\t', stringsAsFactors=F) %>%
      select(experiment_id, master_ccl_id) %>%
      group_by(experiment_id) %>% do({.[1,]}) %>% ungroup
    
    # Load cell line meta data
    d.cline <- read.csv(unz(file.path, 'v20.meta.per_cell_line.txt'), sep='\t', stringsAsFactors=F) %>%
      select(master_ccl_id, ccl_name)
    
    # TODO: Find out if 'navitoclax' filter is appropriate or if it should be expanded to include these
    # Navitoclax variants
    #   [343] "navitoclax"                               "navitoclax:birinapant (1:1 mol/mol)"     
    #   [345] "navitoclax:gemcitabine (1:1 mol/mol)"     "navitoclax:MST-312 (1:1 mol/mol)"        
    #   [347] "navitoclax:piperlongumine (1:1 mol/mol)"  "navitoclax:pluripotin (1:1 mol/mol)"     
    #   [349] "navitoclax:PLX-4032 (1:1 mol/mol)"        "navitoclax:sorafenib (1:1 mol/mol)"  
    
    # Join all datasets and select only relevant fields
    d.ctd <- d.auc %>% 
      inner_join(d.cmpd, by='master_cpd_id') %>%
      inner_join(d.exp, by='experiment_id') %>%
      inner_join(d.cline, by='master_ccl_id') %>%
      filter(tolower(cpd_name) == 'navitoclax') %>%
      #filter(str_detect(tolower(cpd_name), '^navitoclax')) %>%
      select(tumor_id=ccl_name, cpd_name, area_under_curve, experiment_id)
    
    # Verify that records at this point are unique to the cell line and experiment
    # (this was previously a problem when experiment meta data records were repeated)
    if (d.ctd %>% group_by(tumor_id, experiment_id) %>% tally %>% .$n %>% unique != 1)
      stop('Found CTD2 record duplicates')
    
    # Plotting cell lines with multiple AUC values
    # ids <- d.ctd %>% group_by(tumor_id) %>% tally %>% filter(n >= 2) %>% .$tumor_id %>% unique
    # d.ctd %>% filter(tumor_id %in% ids) %>% ggplot(aes(x=tumor_id, y=area_under_curve)) + geom_point() + 
    #    theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    
    d.ctd %>% 
      mutate(tumor_id=toupper(str_trim(tumor_id))) %>%
      group_by(tumor_id) %>% 
      summarise(auc=mean(area_under_curve)) %>% ungroup
  }
  
  # Lazy-load these results (they're expensive to compute), saving them
  # on disk or loading from disk if previously created
  RAW_CACHE$load('ctd2_auc_v2', loader)
}

GetCTD2V1Data <- function(){
  loader <- function(){
    file.path <- RAW_CACHE$download('ctd2_v1_expanded_dataset.zip', CTD2_V1_URL, mode="wb")
    
    read.csv(unz(file.path, 'v10.D3.area_under_conc_curve.txt'), sep='\t', stringsAsFactors=F) %>%
      filter(str_detect(tolower(cpd_name), 'navitoclax')) %>% 
      rename(tumor_id=ccl_name) %>% select(-cpd_name) %>%
      mutate(tumor_id=toupper(str_trim(tumor_id))) %>%
      group_by(tumor_id) %>% summarise(area_under_curve=mean(area_under_curve)) %>% ungroup
  }
  # Lazy-load these results (they're expensive to compute), saving them
  # on disk or loading from disk if previously created
  RAW_CACHE$load('ctd2_auc_v1', loader)
}


#----------------------------------------#
##### COSMIC Constants and Functions #####
#----------------------------------------#

COSMIC_V1_URL <- 'ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-5.0/gdsc_manova_input_w5.csv'

GetCOSMICData <- function(){
  loader <- function(){
    file.path <- RAW_CACHE$download('cosmic_v1_manova_input.csv', COSMIC_V1_URL, mode="wb")
    
    # Read in downloaded data frame
    read.csv(file.path, sep=',', stringsAsFactors=F) %>% 
      # Rename necessary fields
      select(tumor_id=Cell.Line, ic_50=ABT.263_IC_50) %>% 
      # Convert IC 50 to numeric and replace hyphens in Cell Line ID (w/o hyphens, they match HUGO)
      mutate(
        ic_50 = suppressWarnings(as.numeric(ic_50)), 
        tumor_id=str_replace_all(tumor_id, '\\-', '')
      ) %>% 
      # Remove any rows with IC 50 NA
      filter(!is.na(ic_50))
  }
  # Lazy-load these results (they're expensive to compute), saving them
  # on disk or loading from disk if previously created
  RAW_CACHE$load('cosmic_v1', loader)
}
