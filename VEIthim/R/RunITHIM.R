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
source("R/ITHIMCoreAP.R")
source("R/ITHIMCorePA.R")
source("R/RiskFuncsAP.R")
source("R/RiskFuncsPA.R")
source("R/HealthBurden.R")
source("R/Crashes.R")


### local file paths for base VE model run
base_model <- "VE-skatsref-dl-base"  # base model name
ref_year <- 2021  # model base year
base_dir <- file.path(
  "C:", "VisionEval-Dev", "built", "visioneval", "4.2.3", "runtime", "models")

### names of counterfactual VE model scenarios
scenarios <- c(
  "VE-skatsref-dl-density", "VE-skatsref-dl-vmt_charge",
  "VE-skatsref-dl-sov", "VE-skatsref-dl-transit")

# file paths to base model run csvs, etc.
base_model <- file.path(base_dir, base_model)
module_dir <- file.path(getwd(), "inst", "extdata")
input_dir <- file.path(base_model, "inputs", "ithim")  # use base model inputs

# load ithim settings and parse base, forecast years
settings <- loadSettings(module_dir)
years <- parseYears(base_model, ref_year)

# we'll need to loop through each year
for (year in years){
  
  ### LOAD REQUIRED VE OUTPUTS
  print(year)

  # load VE Hhs table
  Hhs <- loadVEHhs(base_model, year)
  
  # transform households table into persons table
  persons <- hhsToPersons(Hhs, input_dir, settings)
  
  # join baseline PM2.5, leisure PA estimates to persons tables
  # persons dt is now nested; second level in persons list
  persons$refPersons <- appendZoneAttributes(persons$persons, year, input_dir)
  
  # now we have established the baseline (reference) population
  # we will now work on scenarios one at a time
  # comparing each to the reference scenario population
  for (scenario in scenarios){
    
    ### PREPARE SCENARIO OUTPUTS
    print(paste(year, scenario))
    
    # instantiate list of (baseline, scenario)
    ref_scenarios <- c("baseline", scenario)
    scenario_dir <- file.path(base_dir, scenario)
    
    # load VE Hhs table for scenario
    scenarioHhs <- loadVEHhs(scenario_dir, year)
    
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
    mmets <- total_mmet(persons, trips, ref_scenarios, input_dir)
    
    # assign a relative risk to each person in the synthetic population for each disease
    # related to physical activity levels and each scenario based on the individual mMET values
    RR_PA_calculations <- gen_pa_rr(mmets, ref_scenarios, module_dir)
    
    # create one dataframe containing both the PA, the AP and the combined PA and AP relative risks
    # (for those diseases affected by both PA and AP) for all people in the synthetic population and all scenarios
    RR_PA_AP_calculations <- combined_rr_ap_pa(
      ind_pa = RR_PA_calculations, ind_ap = RR_AP_calculations, ref_scenarios)
    
    # calculate the health burden (Yll and deaths) for each disease and age and sex category
    # by combining the AP and PA pathways for diseases affected by both AP and PA
    healthBurden <- health_burden(
      ind_ap_pa = RR_PA_AP_calculations, module_dir, input_dir, ref_scenarios)

    ## crash pathway
    # this is a modified version of the California ITHIM code
    # i.e., rates based, non-probabilistic approach
    injuries <- injuries_function(trips, injuryRates, ref_scenarios)
    injuryBurden <- injury_death_to_yll(
      persons$scenarioPersons, injuries, input_dir, ref_scenarios)
    
    ### SUMMARIZE IMPACTS
    healthBurden <- join_hb_and_injury(healthBurden, injuryBurden)
    aggregateHealthBurden(
      healthBurden, persons$refPersons, scenario, scenario_dir)
    
  }  # end of scenarios loop
}  # end of years loop
