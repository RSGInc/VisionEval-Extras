## development code for ITHIM model implementation

# #Load packages and test functions
library(filesstrings)
library(fields)
library(tools)
library(data.table)
library(pracma)
source("R/HarmonizeData.R")
source("R/GenerateScenarios.R")
source("R/ITHIMCore.R")

# params
gender_ratio <- 0.50  # M:F gender ratio

# path to completed VE model run csvs
ve_dir <- file.path("C:", "VisionEval-Dev", "built", "visioneval", "4.2.3", "runtime")
model_dir <- file.path(ve_dir, "models", "VERSPM-metro-reference")
output_dir <- file.path(model_dir, "results", "output", "VERSPM-metro-reference_CSV_202401030956")
files <- Sys.glob(file.path(output_dir, "Household_*"))

# path to ITHIM data inputs (Bzone PM, PA estimates)
input_dir <- file.path(ve_dir, "models", "VERSPM-metro-reference", "inputs")
zonePM <- fread(file.path(input_dir, "bzone_pm.csv"))
setnames(zonePM, c("Geo"), c("Bzone"))
zonePA <- fread(file.path(input_dir, "bzone_pa.csv"))
setnames(zonePA, c("Geo"), c("Bzone"))

# parse years from output dir
getYear <- function(item){
  tail(strsplit(basename(file_path_sans_ext(item)), "_")[[1]], 1)
}
years <- lapply(files, getYear)

# for development, we'll work on first file
file <- files[1]
year <- years[1][[1]]

### step 1: harmonize data
# transform households table into persons table
read_cols = c(
  "HhId", "Bzone",  "HhSize",
  "Age0to14", "Age15to19", "Age20to29", "Age30to54", "Age55to64", "Age65Plus",
  "Dvmt", "TransitTrips", "WalkPMT", "BikePMT"
)
Hhs <- fread(file, select = read_cols)
persons <- hhsToPersons(Hhs, gender_ratio, input_dir)

# join PM2.5, leisure PA estimates to persons tables
persons <- appendZoneAttributesToPersons(persons, zonePM, zonePA, year)

# finally, append counterfactual scenario to persons
trips <- appendCounterfactualScenario(persons)
persons <- trips[[2]]
trips <- trips[[1]]

### step 2: run ITHIM core
# AP pathway
pm_conc <- scenario_pm_calculations(trips, persons)
 
# # TODO: are these needed?
# scenario_pm <- pm_conc$scenario_pm
# pm_conc_pp <- pm_conc$pm_conc_pp
# pm_conc <- NULL
# 
# # Assign relative risks to each person in the synthetic population for each disease
# # related to air pollution and each scenario based on the individual PM exposure levels
# RR_AP_calculations <- gen_ap_rr(pm_conc_pp)
# 
# # TODO: needed?
# if(!constant_mode) pm_conc_pp <- NULL
# 
# # PA pathway
# # calculate total mMETs for each person in the synthetic population
# mmets_pp <- total_mmet(trip_scen_sets)
# 
# # TODO: needed?
# trip_scen_sets <- NULL
# 
# # assign a relative risk to each person in the synthetic population for each disease
# # related to physical activity levels and each scenario based on the individual mMET values
# RR_PA_calculations <- gen_pa_rr(mmets_pp, 
#                                 conf_int = ifelse(constant_mode, TRUE, FALSE))
# 
# # TODO: needed?
# if(!constant_mode) mmets_pp <- NULL
# 
# # crash pathway
# 
# # summary

### step 3: roll-up health impact estimates back to Hhs