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
ref_year <- 2021
base_model <- "VE-skatsref-dl-base"
scenarios <- c(
  "VE-skatsref-dl-density", "VE-skatsref-dl-vmt_charge", "VE-skatsref-dl-sov", "VE-skatsref-dl-transit")
base_dir <- file.path("C:", "VisionEval-Dev", "built", "visioneval", "4.2.3", "runtime", "models")

# load ithim settings
settings <- loadSettings(module_dir)

# file paths to base model run csvs, etc.
base_model <- file.path(base_dir, base_model)
module_dir <- file.path(getwd(), "inst", "extdata")
files <- Sys.glob(file.path(base_model, "results", "output", "*", "Household_*"))

# parse years from output dir
getYear <- function(item){
  as.integer(tail(strsplit(basename(file_path_sans_ext(item)), "_")[[1]], 1))
}
years <- lapply(files, getYear)

# if there is a backcast year (ie, a year before the reference year) scrub it
years <- years[years >= ref_year]

# we'll need to loop through each year
for (year in years){
  
  # identify correct hh table within base scenario outputs directory
  print(year)
  hh_year <- paste0("Household_*", year, "*")
  file <- Sys.glob(file.path(base_model, "results", "output", "*", hh_year))
  input_dir <- file.path(base_model, "inputs", "ithim")  # use base model inputs
  
  # load base model households table
  read_cols = c(
    "HhId", "Bzone",  "HhSize",
    "Age0to14", "Age15to19", "Age20to29", "Age30to54", "Age55to64", "Age65Plus",
    "Dvmt", "TransitPMT", "WalkPMT", "BikePMT"
  )
  Hhs <- fread(file, select = read_cols)
  
  # transform households table into persons table
  persons <- hhsToPersons(Hhs, input_dir, settings)
  
  # join baseline  PM2.5, leisure PA estimates to persons tables
  # persons dt is now nested; second level in persons list
  persons$refPersons <- appendZoneAttributes(persons$persons, year, input_dir)
  
  # now we have established the baseline (reference) population
  # we will now work on scenarios one at a time
  # comparing each to the reference scenario population
  # but first, instantiate a list where we will store results
  scenarioEstimates <- list(
    hhDeaths = NULL,
    hhYlls = NULL,
    zoneDeaths = NULL,
    zoneYlls = NULL)
  
  # loop through scenarios
  for (scenario in scenarios){
    
    ### PREPARE SCENARIO OUTPUTS
    # instantiate list of (baseline, scenario)
    print(paste(year, scenario))
    ref_scenarios <- c("baseline", scenario)
    
    # identify correct hh table within base scenario outputs directory
    scenario_dir <- file.path(base_dir, scenario)
    file <- Sys.glob(file.path(scenario_dir, "results", "output", "*", hh_year))
    scenarioHhs <- fread(file, select = read_cols)
    
    # append this counterfactual scenario to persons
    persons$scenarioPersons <- appendCounterfactualScenario(
      persons$refPersons, scenarioHhs, scenario, settings)
    
    # extract pseudo trips from persons
    trips <- extractTrips(persons$scenarioPersons)
    
    # when looping through the first year, need to develop crash rates
    # we will use these rates for all future time periods
    if (year==years[1]){
      injuryRates <- prepareCrashData(trips, input_dir)
      }
    
    ### RUN CORE ITHIM STEPS
    ## AP pathway
    pm_exposure <- scenario_pm_calculations(trips, persons, input_dir, settings)
    
    # Assign relative risks to each person in the synthetic population for each disease
    # related to air pollution and each scenario based on the individual PM exposure levels
    RR_AP_calculations <- gen_ap_rr(pm_exposure, ref_scenarios, module_dir)
    
    ## PA pathway
    # calculate total mMETs for each person in the synthetic population
    mmets_pp <- total_mmet(persons, trips, ref_scenarios, input_dir)
    
    # assign a relative risk to each person in the synthetic population for each disease
    # related to physical activity levels and each scenario based on the individual mMET values
    RR_PA_calculations <- gen_pa_rr(mmets_pp, ref_scenarios, module_dir)
    
    # create one dataframe containing both the PA, the AP and the combined PA and AP relative risks
    # (for those diseases affected by both PA and AP) for all people in the synthetic population and all scenarios
    RR_PA_AP_calculations <- combined_rr_ap_pa(
      ind_pa = RR_PA_calculations, ind_ap = RR_AP_calculations, ref_scenarios)
    
    # calculate the health burden (Yll and deaths) for each disease and age and sex category
    # by combining the AP and PA pathways for diseases affected by both AP and PA
    healthBurden <- health_burden(
      ind_ap_pa = RR_PA_AP_calculations, module_dir, input_dir, ref_scenarios,
      conf_int = ifelse(constant_mode, TRUE, FALSE))
    
    ## crash pathway
    # this is a modified version of the California ITHIM code
    # i.e., rates based, non-probabilistic approach
    injuries <- injuries_function(trips, injuryRates, ref_scenarios)
    injuryBurden <- injury_death_to_yll(
      persons$scenarioPersons, injuries, input_dir, ref_scenarios)
    
    ### SUMMARIZE IMPACTS
    healthBurden <- join_hb_and_injury(healthBurden, injuryBurden)
    scenarioEstimates <- aggregateHealthBurden(
      healthBurden, persons$refPersons, scenario, scenarioEstimates)
  }  # end of scenarios loop
  
  # write single csv containing all health impact estimates for year
  # TODO: something more elegant?
  write.csv(scenarioEstimates$hhDeaths,
            file.path(base_dir, paste0("Household_mortality_", year, ".csv")))
  write.csv(scenarioEstimates$hhYlls,
            file.path(base_dir, paste0("Household_ylls_", year, ".csv")))
  write.csv(scenarioEstimates$zoneDeaths,
            file.path(base_dir, paste0("Bzone_mortality_", year, ".csv")))
  write.csv(scenarioEstimates$zoneYlls,
            file.path(base_dir, paste0("Bzone_ylls_", year, ".csv")))
}  # end of years loop
