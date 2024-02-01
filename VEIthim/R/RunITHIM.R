## development code for ITHIM model implementation

# #Load packages and test functions
library(filesstrings)
library(fields)
library(tools)
library(data.table)
library(dplyr)
library(pracma)
library(tidyr)
source("R/HarmonizeData.R")
source("R/GenerateScenarios.R")
source("R/ITHIMCore.R")
source("R/RiskFuncs.R")
source("R/HealthBurden.R")
source("R/Crashes.R")


# params
gender_ratio <- 0.50  # M:F gender ratio

# path to completed VE model run csvs
ve_dir <- file.path("C:", "VisionEval-Dev", "built", "visioneval", "4.2.3", "runtime")
model_dir <- file.path(ve_dir, "models", "VERSPM-metro-reference")
output_dir <- file.path(model_dir, "results", "output", "VERSPM-metro-reference_CSV_202401030956")
module_dir <- file.path(getwd(), "inst", "extdata")
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
pm_exposure <- scenario_pm_calculations(trips, persons, input_dir)
scenarios <- pm_exposure[[2]]
pm_exposure <- pm_exposure[[1]]

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
injuryTables <- prepareCrashData(persons, trips, module_dir)
get_all_distances(trips, persons, scenarios)

# summary

### step 3: roll-up health impact estimates back to Hhs