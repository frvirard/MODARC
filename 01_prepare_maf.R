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
maf_lst <- list.files(maf_dr)|>
  str_subset("wo_intergenic", negate = TRUE)

walk(maf_lst, function(maf){
  vroom::vroom(file.path(maf_dr, maf), show_col_types = FALSE)|>
    mutate(sample = str_remove(sample, "_calls.2.tsv"))|>
    vroom::vroom_write(file.path(maf_dr, maf))
})
