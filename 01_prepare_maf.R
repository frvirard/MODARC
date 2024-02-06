### setup ----------------------------------------------------------------------
source("config.R")
remotes::install_github("frvirard/signatR", auth_token = token)

# library ######################################################################
lib_lst <- c("dplyr",
             "purrr",
             "stringr",
             "readr",
             "reticulate")

invisible(purrr::walk(lib_lst, library, character.only = TRUE))
rm(lib_lst)

# functions ####################################################################

# variable #####################################################################
maf_dr <- file.path(data_dr, "maf")

## Prepare maf files -----------------------------------------------------------
#' Samples (T N pairs) from 21 patients were processed by 2 different 
#' persons (Berenice Chavanel and Marie Pierre Cros) 
#' One sample seemed to have been removed due to discrepancy between clinical data table and aliquot number.
#' 2 batches were identified in Vincent caller folder, however all sample 
#' were recalled together into the second batch (GWZ_20221025_MODARC)

signatR::prepare_maf(file.path(file_dr, "vcf"),
                     vcf_extension = ".2.tsv", 
                     get_vaf = FALSE)

## Sample name cleaning --------------------------------------------------------
maf_lst <- list.files(maf_dr)

walk(maf_lst, function(maf){
  vroom::vroom(file.path(maf_dr, maf), show_col_types = FALSE)|>
    mutate(sample = str_remove(sample, "_calls.2.tsv"))|>
    vroom::vroom_write(file.path(maf_dr, maf))
})

## create a liftover version ---------------------------------------------------

chain <- rtracklayer::import.chain(file.path(file_dr,"hg19ToHg38.over.chain"))

dataset_lst <- c("RECA-EU", "KIRC-US")

walk(dataset_lst, function(set){
  message("oo Processing ", set)
  file <- list.files(file_dr, full.names = TRUE)|>
    str_subset(set)
  output <- basename(file)|>
    str_replace(".tsv.gz", "_GRCh38.maf.gz")
  
  if(!file.exists(file.path(maf_dr, output))){
    dt <-  vroom::vroom(file, show_col_types = FALSE)
    message("converting to GRCh38 genome ...")
    
    dt <- dt |>
      mutate(chromosome = if_else(chromosome == "MT", "M", chromosome))|>
      mutate(chromosome = paste0("chr", chromosome))|>
      GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE, 
                                              ignore.strand = TRUE,
                                              seqnames.field = "chromosome",
                                              start.field = "chromosome_start",
                                              end.field = "chromosome_end")|>
      rtracklayer::liftOver(chain)|>
      as_tibble()|>
      mutate(assembly_version = "GRCh38")|>
      select(-c(group, group_name, width, strand))|>
      rename(chromosome = seqnames, 
             chromosome_start = start, 
             chromosome_end = end)
    
    message("done")
    
    vroom::vroom_write(dt, file.path(maf_dr, output))
  }else{
    message("file already detected, skipping ...")
  }
})