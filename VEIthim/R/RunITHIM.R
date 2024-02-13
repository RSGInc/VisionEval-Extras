## development code for ITHIM model implementation

# #Load packages and test functions
library(filesstrings)
library(fields)
library(tools)
library(data.table)
library(dplyr)
library(pracma)
library(tidyr)
source("R/LoadSettings.R")
source("R/HarmonizeData.R")
source("R/GenerateScenarios.R")
source("R/ITHIMCore.R")
source("R/RiskFuncs.R")
source("R/HealthBurden.R")
source("R/Crashes.R")


# path to completed VE model run csvs
ve_dir <- file.path("C:", "VisionEval-Dev", "built", "visioneval", "4.2.3", "runtime")
model_dir <- file.path(ve_dir, "models", "VERSPM-metro-reference")
output_dir <- file.path(model_dir, "results", "output", "VERSPM-metro-reference_CSV_202401030956")
module_dir <- file.path(getwd(), "inst", "extdata")

# load ithim settings
settings <- loadSettings(module_dir)

# path to completed VE model run csvs
files <- Sys.glob(file.path(output_dir, "Household_*"))

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
  "Dvmt", "TransitPMT", "WalkPMT", "BikePMT"
)
Hhs <- fread(file, select = read_cols)
persons <- hhsToPersons(Hhs, input_dir, settings)

# join PM2.5, leisure PA estimates to persons tables
# persons dt is now nested, second level in persons list
persons$persons <- appendZoneAttributes(persons$persons, year, input_dir)

# finally, append counterfactual scenario to persons
persons$persons <- appendCounterfactualScenario(persons$persons)

# and parse scenario names from persons
scenarios <- unique(persons$persons$scenario)

# and extract trips
trips <- extractTrips(persons$persons)

### step 2: run ITHIM core
# AP pathway
pm_exposure <- scenario_pm_calculations(trips, persons, input_dir, settings)

# Assign relative risks to each person in the synthetic population for each disease
# related to air pollution and each scenario based on the individual PM exposure levels
RR_AP_calculations <- gen_ap_rr(pm_exposure, scenarios, module_dir)

# PA pathway
# calculate total mMETs for each person in the synthetic population
mmets_pp <- total_mmet(persons, trips, scenarios, input_dir)

# assign a relative risk to each person in the synthetic population for each disease
# related to physical activity levels and each scenario based on the individual mMET values
RR_PA_calculations <- gen_pa_rr(mmets_pp, scenarios, module_dir)

# create one dataframe containing both the PA, the AP and the combined PA and AP relative risks
# (for those diseases affected by both PA and AP) for all people in the synthetic population and all scenarios
RR_PA_AP_calculations <- combined_rr_ap_pa(
  ind_pa = RR_PA_calculations, ind_ap = RR_AP_calculations, scenarios)

# calculate the health burden (Yll and deaths) for each disease and age and sex category
# by combining the AP and PA pathways for diseases affected by both AP and PA
hb_AP_PA <- health_burden(ind_ap_pa = RR_PA_AP_calculations, module_dir, input_dir, scenarios,
                          conf_int = ifelse(constant_mode, TRUE, FALSE))

# crash pathway
injuryRates <- prepareCrashData(persons, trips, module_dir)
injuries <- injuries_function(trips, injuryRates, scenarios)
browser()

# summary

### step 3: roll-up health impact estimates back to Hhs