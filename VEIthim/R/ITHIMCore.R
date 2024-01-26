#' Calculate total AP exposure per person
#' 
#' Calculate total AP exposure per person based on population and personal travel
#' 
#' This function performs the following steps:
#' 
#' \itemize{
#' \item the ventilation rates per mode are defined - these parameters are fixed
#' 
#' \item the exposure factor rates by activity are defined - these parameters are fixed
#' 
#' \item calculate pm concentration not related to transport
#' 
#' \item calculate PM emission factors for each mode by dividing total emissions by distances travelled
#' 
#' \item calculate PM emissions for each mode in each scenario by multiplying the scenario distance
#'   by the emission factors
#' 
#' \item for modes without any assigned distance, use the PM emissions from the VEHICLE_INVENTORY instead
#' 
#' \item calculate the total PM concentrations for each scenario
#' 
#' \item add ventilation and exposure factors to the trip set by stage mode
#' 
#' \item add total scenario PM concentrations to the trip set
#' 
#' \item add the inhaled air and total PM (in micro grams) to the trip set
#' 
#' \item define the amount of time per day spent as leisure sedentary screen time, 
#'   non-discretionary time and other time - fixed
#' 
#' \item add total time spent travelling by each participant to the trip set
#'
#' \item calculate the sleep and the rest ventilation rates
#' 
#' \item for each participant in the synthetic population (with travel component), 
#'   calculate the total air inhaled, the total PM inhaled and 
#'   the total PM concentration inhaled for each scenario
#'
#' \item assign all participants in the synthetic population without travel component,
#'   the baseline or scenario PM concentrations
#'
#' \item join all people with and without travel in the synthetic population
#' }
#' 
#' 
#' @param dist total distance travelled by mode by the population for all scenarios
#' @param trip_scen_sets trips data frame of all trips from all scenarios
#' 
#' @return background PM concentration for baseline and all scenarios
#' @return total AP exposure per person in the synthetic population (for baseline and scenarios)
#' 
#' @export


scenario_pm_calculations <- function(trips, Persons){
  
  # TODO: move into setting
  PM_TRANS_SHARE= 0.101
  
  # concentration contributed by non-transport share (remains constant across the scenarios)
  background_pm <- Persons[scenario=="baseline", PM25]
  non_transport_pm <- background_pm*(1 - PM_TRANS_SHARE)
  
  # Exposure factor rate by activity (the ratio between that mode’s PM2.5 and the background’s PM2.5)
  exp_facs <- list(
    auto = 2.5, 
    transit = 1.9, 
    bike = 2.0, 
    walk = 1.6
  )
  
  # add mode speeds
  mode_speeds <- list(
    auto = 14.4,
    transit = 10, 
    bike = 7.2,
    walk = 2.5
    )
  
  # and pm emissions inventory values
  pm_inventory <- list(
    auto = 0.081980353,
    transit = 0.122466774, 
    bike = 0,
    walk = 0
  )
  
  mode_inventory <- data.frame(
    exp_fac = unlist(exp_facs),
    speed = unlist(mode_speeds),
    PM_emission_inventory = unlist(pm_inventory),
    stringsAsFactors = F)

  # we assume that the travel predicted by VE represents all travel in the region
  # ie, we do not need to add travel not covered in synthetic trip set
  ## adding in travel not covered in the synthetic trip set, based on distances traveled relative to car, set in VEHICLE_INVENTORY

  # total distance traveled by each mode
  browser()
  distances_by_mode <- Persons[, lapply(.SD, sum, na.rm=TRUE),by=.(scenario), .SDcols=c("auto", "transit", "walk", "bike") ]
  scenarios <- distances_by_mode$scenario  # stash scenario names
  counterfactual_scenarios <- scenarios[2:length(scenarios)]  # scrub "baseline"
  distances_by_mode <- transpose(
    distances_by_mode[, c("auto", "transit", "walk", "bike")],
    keep.names = "mode")
  
  # transpose creates temporary names we want to replace
  temp_names <- names(distances_by_mode)
  temp_names <- temp_names[2:length(temp_names)]  # scrub "mode"
  setnames(distances_by_mode, temp_names, scenarios)

  # convert to data.frame for indexing
  distances_by_mode <- as.data.frame(distances_by_mode)
  rownames(distances_by_mode) <- distances_by_mode$mode
  distances_by_mode$mode <- NULL
  
  # and now we can append to mode_inventory
  mode_inventory$total_dist <- distances_by_mode$total_dist

  ## get emission factor by dividing inventory by baseline distance. (We don't need to scale to a whole year, as we are just scaling the background concentration.)
  # this is just emission inventory divided by total travel
  baseline_emission_factors <- mode_inventory$PM_emission_inventory / distances_by_mode$baseline
  trans_emissions <- distances_by_mode[scenarios] * t(repmat(baseline_emission_factors, length(scenarios), 1))
  
  # TODO: figure out something more elegant
  # transit emissions shouldn't scale directly with transit trips (i.e., based on service, not person-trips)
  # so replace all scenario value with baseline value
  for(counterfactual in counterfactual_scenarios){
    trans_emissions["transit" , counterfactual] <- trans_emissions["transit" , "baseline"]
  }
  
  ## scenario travel pm2.5 calculated as relative to the baseline
  baseline_pm <- sum(trans_emissions["baseline"], na.rm = T)
  
  ## in this sum, the non-transport pm is constant; the transport emissions scale the transport contribution (PM_TRANS_SHARE) to the base level (PM_CONC_BASE)
  for(counterfactual in counterfactual_scenarios){
    Persons[scenario==counterfactual, "PM25"] <- 
      non_transport_pm + PM_TRANS_SHARE * background_pm * sum(trans_emissions[[counterfactual]], na.rm = T)/baseline_pm
  }
  
  # Join trip_set and exponent factors df
  browser()
  trips <- dplyr::left_join(trips, exp_facs, 'stage_mode')
  
  #----
  # Dan: These new lines of code are for the ventilation rate
  # Dan: MET values for each mode/activity. These values come from the Compendium
  # Dan: Lambed told me that we need v-rates for everyone. I assigned to 
  # auto_rickshaw and other the same MET values as bus and cycle.
  # TODO: I don't know if we should add MET values for drivers in ghost trips
  # (bus, car, motorcycle, truck). Right now these are not included because they
  # have participant_id = 0.
  met_df <- data.frame(
    stage_mode = c(
      "car", "transit", "bike", "walk", "sleep",
      "moderate", "vigorous", "leisure", "light_activities"
    ),
    met = c(
      CAR_DRIVER_MET, PASSENGER_MET, CYCLING_MET, WALKING_MET, 0.95,
      MODERATE_PA_MET, VIGOROUS_PA_MET,
      SEDENTARY_ACTIVITY_MET, LIGHT_ACTIVITY_MET
    ),
    compendium_code = c(
      "16010", "16016", "01011", "16060", "07030",
      "16030", "GPAQ", "GPAQ", "05080", "05080"
    )
  )
  
  # Dan: Extract people from the synthetic population to calculate their 
  # ventilation rates
  people_for_vent_rates <- trip_set %>% filter(participant_id != 0) %>% 
    distinct(participant_id, .keep_all=T) %>% 
    dplyr::select(participant_id, age, sex)
  
  # Dan: Adding these new columns to the trip set
  trip_set <- trip_set %>% 
    left_join(people_for_vent_rates, by = 'participant_id')
  
  # Dan: Adding MET values [dimensionless] for each mode
  trip_set <- trip_set %>% left_join(met_df %>% dplyr::select(stage_mode, met), 
                                     by = 'stage_mode')
  
  # Dan: Calculate new variables for each stage
  #   - vo2 = oxygen uptake [lt/min]
  #   - pct_vo2max = upper limit of VO2max in percentage [%]
  #   - upper_vo2max = permitted upper limit of VO2 [lt/min]
  #   - adj_vo2 = adjusted oxygen uptake [lt/min]
  #   - e_ijk = error to account variance between activities
  #   - log_vent_rate = log of ventilation rate (empirical equation)
  #   - vent_rate = ventilation rate by removing the log in the empirical equation [lt/min]
  #   - v_rate = ventilation rate in different units of measurement [m3/h]
  trip_set <- trip_set %>% 
    rowwise() %>% 
    mutate(
      # Ventilation rate for travel modes
      vo2 = ecf * met * rmr, 
      pct_vo2max = ifelse(stage_duration < 5, 100,
                          ifelse(stage_duration > 540, 33,
                                 121.2 - (14 * log(stage_duration)))),
      upper_vo2max = vo2max * (pct_vo2max / 100),
      adj_vo2 = ifelse(vo2 > upper_vo2max, upper_vo2max, vo2),
      # Draw sample from normal distribution taking into account the variance 
      # between activities
      e_ijk = rnorm(1, 0, sd_test_level),
      log_vent_rate = intercept_a + (slope_b * log(adj_vo2/body_mass)) + d_k + e_ijk,
      vent_rate = exp(log_vent_rate) * body_mass,
      v_rate = vent_rate * 60 / 1000
    )
  #----
  
  # Create df with scenarios and concentration
  conc_pm_df <- data.frame(scenario = unique(trip_set$scenario),
                           conc_pm = conc_pm)
  
  # Join trip_set with PM concentration df
  trip_set <- left_join(trip_set, conc_pm_df, by = 'scenario')
  
  # cubic meters of air inhaled are the product of the ventilation rate and the 
  # time (hours/60) spent travelling by that mode
  trip_set$air_inhaled <- trip_set$stage_duration / 60 * trip_set$v_rate
  
  # PM inhaled (micro grams) = duration * ventilation rate * exposure rates * concentration
  trip_set$pm_inhaled <- trip_set$stage_duration / 60 * trip_set$v_rate * trip_set$e_rate * trip_set$conc_pm
  
  # Calculate total_travel_time_hrs of stage_duration
  trip_set <- trip_set %>% 
    group_by(participant_id, scenario) %>% 
    mutate(total_travel_time_hrs = sum(stage_duration) / 60) %>% 
    # total_typical_time_rem_hrs <- lt_sed_time_hrs + nd_time_hrs + other_time_hrs,
    # remaining_time_hrs = 24 - total_travel_time_hrs,
    # lt_sed_time_prop_hrs <- lt_sed_time_hrs / total_typical_time_rem_hrs * remaining_time_hrs,
    # nd_time_prop_hrs <- nd_time_hrs / total_typical_time_rem_hrs * remaining_time_hrs,
    # other_time_prop_hrs <- other_time_hrs / total_typical_time_rem_hrs * remaining_time_hrs) |> 
    ungroup()
  
  
  # Extract mets for sleep, rest, moderate and vigorous
  sleep_met <- met_df %>% filter(stage_mode == 'sleep') %>% 
    dplyr::select(met) %>% as.numeric()
  rest_met <- met_df %>% filter(stage_mode == 'rest') %>% 
    dplyr::select(met) %>% as.numeric()
  moderate_met <- met_df %>% filter(stage_mode == 'moderate') %>% 
    dplyr::select(met) %>% as.numeric()
  vigorous_met <- met_df %>% filter(stage_mode == 'vigorous') %>% 
    dplyr::select(met) %>% as.numeric()
  leisure_met <- met_df %>% filter(stage_mode == 'leisure') %>% 
    dplyr::select(met) %>% as.numeric()
  light_met <- met_df %>% filter(stage_mode == 'light_activities') %>% 
    dplyr::select(met) %>% as.numeric()
  
  #----
  # Dan: We are saying that time spent for any person is:
  # Sleep: 8.3 hours
  # Leisure sedentary screen time: 3.15 hours
  # Light activities: 10.75 hours
  # For a total of 22.2 hours
  # Since it is possible that the unknown time is less than 22.2, I am going
  # to use proportions for now to get values for these activities.
  # There are a couple of rules to consider:
  # - When sleep time is less than 6 hours, we have to assign it to 6 and split
  # proportionally the remaining time in leisure and light activities
  # - When the known time is greater than 18 hours, we assign it to 18 and sleep
  # to 6 hours. Leisure and light activities will be zero in this case
  
  sleep_hours <- 8.3
  leisure_hours <- 3.15
  light_hours <- 10.75
  
  # Calculate total air and pm inhaled in each person
  synth_pop <- trip_set  %>%  
    filter(participant_id != 0) %>% 
    group_by(participant_id, scenario) %>% 
    # reframe instead of summarise in the latest version
    reframe(
      total_travel_time_hrs = max(total_travel_time_hrs, na.rm = T),
      travel_air_inhaled = sum(air_inhaled, na.rm = T),
      travel_pm_inhaled = sum(pm_inhaled, na.rm = T),
      work_ltpa_marg_met = max(work_ltpa_marg_met, na.rm=T)
    ) %>% 
    distinct(participant_id, scenario, .keep_all = T) %>% 
    # Merge participant information to calculate ventilation rates for sleep and rest
    left_join(people_for_vent_rates, by = 'participant_id') %>% 
    # Merge scenario concentration
    left_join(conc_pm_df, by = 'scenario') %>% 
    # Calculate vent_rates and air inhaled in the same way as with travel modes
    rowwise() %>% 
    mutate(
      # Transforming work_ltpa_marg_met to a daily value
      daily_work_ltpa_marg_met = work_ltpa_marg_met / 7,
      # Calculate time spent in moderate and vigorous activities from work_ltpa_marg_met
      # I subtract 1 because the initial values are MET not Marginal MET
      time_moderate = (daily_work_ltpa_marg_met * MODERATE_PA_CONTRIBUTION) / (moderate_met - 1),
      time_vigorous = (daily_work_ltpa_marg_met * (1 - MODERATE_PA_CONTRIBUTION)) / (vigorous_met - 1),
      # Calculate known time
      known_time = total_travel_time_hrs + time_moderate + time_vigorous,
      known_time = ifelse(known_time > 18, 18, known_time), # Just to avoid outliers
      # Calculate unknown time
      unknown_time = 24 - known_time,
      
      # Calculate ventilation rate for each activity (durations are calculated from the unknown time)
      ## Sleep
      sleep_duration = (unknown_time * (sleep_hours/(sleep_hours + leisure_hours + light_hours))) * 60, # In minutes
      ### Conditional to check if sleep duration is less than 6 hours
      sleep_less_6h = ifelse(sleep_duration < (6 * 60), 1, 0),
      ### Assign 6 hours when sleep time is less than 6
      sleep_duration = ifelse(sleep_less_6h == 1, (6 * 60), sleep_duration),
      vo2_sleep = ecf * sleep_met * rmr,
      pct_vo2max_sleep = ifelse(sleep_duration < 5, 100,
                                ifelse(sleep_duration > 540, 33,
                                       121.2 - (14 * log(sleep_duration)))),
      upper_vo2max_sleep = vo2max * (pct_vo2max_sleep / 100),
      adj_vo2_sleep = ifelse(vo2_sleep > upper_vo2max_sleep, 
                             upper_vo2max_sleep, vo2_sleep),
      ### Draw sample from normal distribution taking into account the variance 
      ### between activities
      e_ijk_sleep = rnorm(1, 0, sd_test_level),
      log_vent_rate_sleep = intercept_a + (slope_b * log(adj_vo2_sleep/body_mass)) + d_k + e_ijk_sleep,
      vent_rate_sleep = exp(log_vent_rate_sleep) * body_mass,
      v_rate_sleep = vent_rate_sleep * 60 / 1000,
      ### Calculate air and pm inhaled
      sleep_air_inhaled = sleep_duration / 60 * v_rate_sleep,
      sleep_pm_inhaled = sleep_duration / 60 * v_rate_sleep * conc_pm,
      
      
      ## Leisure sedentary screen time
      #### Leisure has a correction if sleep time is less than 6 hours
      leisure_duration = ifelse(sleep_less_6h == 1,
                                (unknown_time - 6) * (leisure_hours/(leisure_hours + light_hours)) * 60,# In minutes
                                unknown_time * (leisure_hours/(sleep_hours + leisure_hours + light_hours))) * 60, # In minutes
      leisure_duration = ifelse(leisure_duration < 0, 0, leisure_duration), # This happened in 35 rows out of the 132k
      vo2_leisure = ecf * leisure_met * rmr,
      pct_vo2max_leisure = ifelse(leisure_duration < 5, 100,
                                  ifelse(leisure_duration > 540, 33,
                                         121.2 - (14 * log(leisure_duration)))),
      upper_vo2max_leisure = vo2max * (pct_vo2max_leisure / 100),
      adj_vo2_leisure = ifelse(vo2_leisure > upper_vo2max_leisure, 
                               upper_vo2max_leisure, vo2_leisure),
      ## Draw sample from normal distribution taking into account the variance 
      ## between activities
      e_ijk_leisure = rnorm(1, 0, sd_test_level),
      log_vent_rate_leisure = intercept_a + (slope_b * log(adj_vo2_leisure/body_mass)) + d_k + e_ijk_leisure,
      vent_rate_leisure = exp(log_vent_rate_leisure) * body_mass,
      v_rate_leisure = vent_rate_leisure * 60 / 1000,
      # Calculate air and pm inhaled
      leisure_air_inhaled = leisure_duration / 60 * v_rate_leisure,
      leisure_pm_inhaled = leisure_duration / 60 * v_rate_leisure * conc_pm,
      
      
      ## Light activities
      #### Light activities has a correction if sleep time is less than 6 hours
      light_duration = ifelse(sleep_less_6h == 1,
                              (unknown_time - 6) * (light_hours/(leisure_hours + light_hours)) * 60,# In minutes
                              unknown_time * (light_hours/(sleep_hours + leisure_hours + light_hours))) * 60, # In minutes
      light_duration = ifelse(light_duration < 0, 0, light_duration), # This happened in 35 rows out of the 132k
      vo2_light = ecf * light_met * rmr,
      pct_vo2max_light = ifelse(light_duration < 5, 100,
                                ifelse(light_duration > 540, 33,
                                       121.2 - (14 * log(light_duration)))),
      upper_vo2max_light = vo2max * (pct_vo2max_light / 100),
      adj_vo2_light = ifelse(vo2_light > upper_vo2max_light, 
                             upper_vo2max_light, vo2_light),
      ## Draw sample from normal distribution taking into account the variance 
      ## between activities
      e_ijk_light = rnorm(1, 0, sd_test_level),
      log_vent_rate_light = intercept_a + (slope_b * log(adj_vo2_light/body_mass)) + d_k + e_ijk_light,
      vent_rate_light = exp(log_vent_rate_light) * body_mass,
      v_rate_light = vent_rate_light * 60 / 1000,
      # Calculate air and pm inhaled
      light_air_inhaled = light_duration / 60 * v_rate_light,
      light_pm_inhaled = light_duration / 60 * v_rate_light * conc_pm,
      
      # Total air and pm inhaled
      total_air_inhaled = travel_air_inhaled + sleep_air_inhaled + leisure_air_inhaled + light_air_inhaled,
      total_pm_inhaled = travel_pm_inhaled + sleep_pm_inhaled + leisure_pm_inhaled + light_pm_inhaled,
      
      # Calculate pm / air ratio
      conc_pm_inhaled = total_pm_inhaled / total_air_inhaled
    ) 
  
  # Change to wide format
  synth_pop <- synth_pop %>% 
    dplyr::select(participant_id, scenario, conc_pm_inhaled) %>% 
    pivot_wider(names_from = 'scenario', values_from = 'conc_pm_inhaled') %>% 
    # Rename columns
    rename_at(vars(starts_with(c("base", "sc"))), 
              ~ paste0("pm_conc_", SCEN_SHORT_NAME))
  
  # Get all participants without any travel (in the travel survey)
  id_wo_travel <- SYNTHETIC_POPULATION |> 
    filter(!participant_id %in% trip_set$participant_id)
  
  # Assign all participants without travel baseline + scenario specific base concentration
  id_wo_travel <- cbind(id_wo_travel |> 
                          dplyr::select(-work_ltpa_marg_met), conc_pm_df |> 
                          pivot_wider(names_from = "scenario", values_from = "conc_pm"))
  # Rename columns
  id_wo_travel <- id_wo_travel |> 
    rename_at(vars(starts_with(c("base", "sc"))), 
              ~ paste0("pm_conc_", SCEN_SHORT_NAME))
  
  # Join demographics info from trip_set 
  synth_pop <- left_join(trip_set |> 
                           filter(participant_id != 0) |> 
                           dplyr::select(participant_id, age, sex, age_cat) |> 
                           distinct(), 
                         synth_pop, by = 'participant_id')
  
  # Combine people with and without trips
  synth_pop <- dplyr::bind_rows(synth_pop, id_wo_travel)
  
  # Convert data type to integer
  synth_pop$participant_id <- as.integer(synth_pop$participant_id)
  
  # Return list with concentration and per person PM2.5 exposure (unit: ug/m3)
  list(scenario_pm = conc_pm, pm_conc_pp = as.data.frame(synth_pop))
}