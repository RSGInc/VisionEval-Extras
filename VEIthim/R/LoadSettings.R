# loads csv file containing ithim parameters
loadSettings <- function(module_dir){
  
  # read settings
  settings <- fread(file.path(module_dir, "ithim_settings.csv"))
  
  return(settings)
}
  

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