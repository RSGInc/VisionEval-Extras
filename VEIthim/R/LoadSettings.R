# loads csv file containing ithim parameters

loadSettings <- function(module_dir){
  
  # read settings
  settings <- fread(file.path(module_dir, "ithim_settings.csv"))
  
  return(settings)
}
  
