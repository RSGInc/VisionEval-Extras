# this module contains helper funcs for settings/configuration management


#' Load settings
#'
#' Helper function to load ithim settings
#'
#' This function performs the following steps:
#'
#' \itemize{
#'  \item loads settings
#'    }
#'
#' @param module_dir local ithim module directory
#'
#' @return data.table containing ithim settings
#'
#' @export


# loads csv file containing ithim parameters
loadSettings <- function(module_dir){
  
  # read settings
  settings <- fread(file.path(module_dir, "settings.csv"))
  
  return(settings)
}
  

#' Parse year
#'
#' Helper function to extract list of years from VE output csv directory and 
#' scrub years before the reference year (ie, backcasting years)
#'
#' This function performs the following steps:
#'
#' \itemize{
#'  \item identifies all household files in a VE output csv directory
#'  
#'  \item parses year from each file in the list
#'  
#'  \item removes years that are before the reference year
#'    }
#'
#' @param base_model local VE model directory
#' @param ref_year reference base) year for ithim model
#'
#' @return list of years for which a VE model exists (excluding those before 
#' the reference year)
#'
#' @export


parseYears <- function(base_model, ref_year){
  
  # parse years from output dir
  files <- Sys.glob(file.path(base_model, "results", "output", "*", "Household_*"))
  getYear <- function(item){
    as.integer(tail(strsplit(basename(file_path_sans_ext(item)), "_")[[1]], 1))
  }
  years <- lapply(files, getYear)
  
  # if there is a backcast year (ie, a year before the reference year) scrub it
  years <- years[years >= ref_year]
  return(years)
}