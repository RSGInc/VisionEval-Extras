#' Predict injuries
#'
#' Predict injuries for baseline and scenarios based on Poisson regression model fitted
#' on baseline fatality counts and distances
#'
#' This function uses the Poisson regression model built in the \code{\link{distances_for_injury_function()}} to predict fatality
#' counts for the Baseline and all the scenarios. It performs the following steps:
#'
#' \itemize{
#' \item create an injuries data frame containing all the distances travelled by mode, age, sex and scenario
#'
#' \item predict the fatalities for each strike and casualty mode combination, age and sex category
#'   and each scenario (incl Baseline). If the sample mode is set to 'constant' (and not 'sample'),
#'   we also predict upper and lower confidence interval boundaries
#'
#' \item create a whw_temp list containing the total predicted fatality counts for each casualty
#'   and strike mode pair for each scenario split into whw and nov matrices and, for the
#'   constant mode also give the upper and lower confidence interval limit predictions
#'
#' \item create an injuries2 data frame containing the total predicted fatality counts
#'   for each casualty mode by age and sex for each scenario. This dataframe also
#'   contains total death per age and sex category and, for the constant mode the
#'   upper and lower total death predictions of the confidence interval.
#' }
#'
#'
#'
#' @param true_distances data frame containing population distances for each scenario
#' @param injuries_list list of dataframes set up with scenario specific information to supply to regression model for prediction
#' @param reg_model Poisson injury regression model
#' @param constant_mode whether or not we are in constant (vs sampling) mode
#'
#' @return injuries2 - dataframe containing predicted fatality counts for each casualty mode by age and sex and for each scenario, plus confidence interval limits for constant mode
#' @return whw_temp - list containing the fatality predictions for each casualty and strike mode pair split into whw and nov matrices for each scenario. Upper and lower confidence interval predictions are also included for the constant mode
#'
#' @export


distances_for_injury_function <- function(journeys, dist) {
  browser()
  distances <- spread(journeys, mode, tot_dist, fill = 0)
  
  if ("walk_to_pt" %in% names(distances)) { # add walk to pt to pedestrian distance to have one measure for walking
    distances$walk <- distances$walk + distances$walk_to_pt
  }
  ## add all car related distances to car distance
  distances$auto <- rowSums(distances[, colnames(distances) %in% c("auto", "taxi", "shared_auto", "shared_taxi")])
  distances <- distances[, -which(names(distances) %in% c("taxi", "shared_auto", "shared_taxi", "walk_to_pt"))]
  
  true_distances_0 <- distances
  true_distances_0$sex_age <- paste0(true_distances_0$sex, "_", true_distances_0$age_cat) # add one sex and age column
  
  # add bus driver distances to total bus distances to cover all by people travelling by bus
  if (ADD_BUS_DRIVERS) true_distances_0$bus <- true_distances_0$bus + true_distances_0$bus_driver
  
  true_distances <- true_distances_0
  
  # find mode names
  mode_names <- names(true_distances)[!names(true_distances) %in% c("age_cat", "scenario", "sex_age", "sex")]
  
  
  injury_table <- INJURY_TABLE
  
  # add injury_reporting rate
  for (type in INJURY_TABLE_TYPES) {
    INJURY_TABLE[[type]]$injury_reporting_rate <- INJURY_REPORTING_RATE
  }
  
  INJURY_TABLE <<- INJURY_TABLE
  
  
  ## add distance columns for baseline to injury data
  injuries_for_model <- add_distance_columns(injury_table,
                                             mode_names,
                                             true_distances_0,
                                             dist,
                                             scenarios = SCEN[1]
  )
  
  # determine whether there are any age and sex combinations for which we don't have any distance information
  zero_dist <- list()
  zero_dist_pos_inj <- list()
  
  # flag to highlight if some age and gender categories have at least one fatality but zero distance for a strike and casualty mode combination
  zero_dist_flag <- F
  
  # finds cas and strike mode combinations for which there exist zero distances
  for (type in INJURY_TABLE_TYPES) {
    zero_dist[[type]] <- subset(injuries_for_model$baseline[[type]], strike_distance == 0 | cas_distance == 0)
    zero_dist_pos_inj[[type]] <- subset(zero_dist[[type]], count > 0)
    if (nrow(zero_dist_pos_inj[[type]]) > 0) {
      zero_dist_flag <- T
    }
  }
  
  
  # if there exists at least one age and sex category with at least one fatality but with
  # zero distance for a strike and casualty mode combination, aggregate by age and sex
  if (zero_dist_flag == T) {
    for (type in INJURY_TABLE_TYPES) {
      injuries_df <- injuries_for_model$baseline[[type]]
      setDT(injuries_df)
      injuries_for_model$baseline[[type]] <- as.data.frame(injuries_df[, .(
        count = sum(count), weight = mean(weight), strike_distance_sum = mean(strike_distance_sum),
        cas_distance_sum = mean(cas_distance_sum)
      ), by = c("cas_mode", "strike_mode")])
      injuries_for_model$baseline[[type]]$strike_distance <- injuries_for_model$baseline[[type]]$strike_distance_sum
      injuries_for_model$baseline[[type]]$cas_distance <- injuries_for_model$baseline[[type]]$cas_distance_sum
    }
  }
  
  # remove any injuries for which we don't have either casualty or strike mode distance
  for (type in INJURY_TABLE_TYPES) {
    injuries_for_model$baseline[[type]] <- subset(injuries_for_model$baseline[[type]], strike_distance > 0 & cas_distance > 0)
  }
  
  
  # create list where each element contains all age, sex, casualty and strike mode combination for
  # the baseline and all scenarios
  scenario_injury_table <- list()
  for (type in INJURY_TABLE_TYPES) {
    scenario_injury_table[[type]] <- expand.grid(
      age_cat = unique(DEMOGRAPHIC$age),
      cas_gender = unique(DEMOGRAPHIC$sex),
      cas_mode = unique(injuries_for_model[[1]][[type]]$cas_mode),
      strike_mode = unique(injuries_for_model[[1]][[type]]$strike_mode)
    )
  }
  
  
  # add distance information for the baseline and all scenarios
  injuries_list <- add_distance_columns(injury_table = scenario_injury_table, mode_names, true_distances_0, dist)
  
  
  for (n in 1:(NSCEN + 1)) { # loop through baseline and all scenarios
    for (type in INJURY_TABLE_TYPES) {
      # remove zero distances
      injuries_list[[n]][[type]] <- subset(
        injuries_list[[n]][[type]],
        strike_distance > 0 & cas_distance > 0
      )
      injuries_list[[n]][[type]]$injury_gen_age <- apply(
        cbind(
          as.character(injuries_list[[n]][[type]]$cas_gender),
          as.character(injuries_list[[n]][[type]]$age_cat)
        ), 1,
        function(x) paste(x, collapse = "_")
      ) # create an age sex column
      # remove strike and cas mode pairs where cas mode = strike mode as these are considered nov accidents
      injuries_list[[n]][[type]]$cas_strike_mode <- apply(
        cbind(
          as.character(injuries_list[[n]][[type]]$cas_mode),
          as.character(injuries_list[[n]][[type]]$strike_mode)
        ), 1,
        function(x) paste(x, collapse = "_")
      )
      injuries_list[[n]][[type]] <- injuries_list[[n]][[type]] %>%
        filter(cas_strike_mode != "car_car" & cas_strike_mode != "bus_bus" &
                 cas_strike_mode != "motorcycle_motorcycle" &
                 cas_strike_mode != "cycle_cycle" &
                 cas_strike_mode != "truck_truck")
    }
  }
  
  ## determine safety in number coefficients and add as columns to injuries_list and injuries_for_model
  if (CALL_INDIVIDUAL_SIN == F) { # if we have the same coefficients for all modes
    CAS_EXPONENT <<- SIN_EXPONENT_SUM * CASUALTY_EXPONENT_FRACTION
    STR_EXPONENT <<- SIN_EXPONENT_SUM - CAS_EXPONENT
    
    # when running in sampling mode, ensure that cas_exponent and str_exponent are always below 1
    if (CAS_EXPONENT > 1) CAS_EXPONENT <<- 1
    if (STR_EXPONENT > 1) STR_EXPONENT <<- 1
    
    for (type in INJURY_TABLE_TYPES) { # add the exponents to the the baseline table and also the scenario tables
      injuries_for_model$baseline[[type]]$cas_exponent_col <- CAS_EXPONENT
      injuries_for_model$baseline[[type]]$str_exponent_col <- STR_EXPONENT
      
      for (n in 1:(NSCEN + 1)) {
        injuries_list[[n]][[type]]$cas_exponent_col <- CAS_EXPONENT
        injuries_list[[n]][[type]]$str_exponent_col <- STR_EXPONENT
      }
    }
  }
  
  
  if (CALL_INDIVIDUAL_SIN == T) { # assign coefficients depending on strike and victim modes
    
    injuries_for_model$baseline$whw$cas_exponent_col <- SIN_EXPONENT_SUM_VEH * CASUALTY_EXPONENT_FRACTION_VEH
    injuries_for_model$baseline$whw$str_exponent_col <- SIN_EXPONENT_SUM_VEH - (SIN_EXPONENT_SUM_VEH * CASUALTY_EXPONENT_FRACTION_VEH)
    
    
    injuries_for_model$baseline$whw$cas_exponent_col[injuries_for_model$baseline$whw$cas_mode == "cycle"] <- SIN_EXPONENT_SUM_CYCLE * CASUALTY_EXPONENT_FRACTION_CYCLE
    injuries_for_model$baseline$whw$str_exponent_col[injuries_for_model$baseline$whw$cas_mode == "cycle"] <- SIN_EXPONENT_SUM_CYCLE - (SIN_EXPONENT_SUM_CYCLE * CASUALTY_EXPONENT_FRACTION_CYCLE)
    
    injuries_for_model$baseline$whw$cas_exponent_col[injuries_for_model$baseline$whw$cas_mode == "pedestrian"] <- SIN_EXPONENT_SUM_PED * CASUALTY_EXPONENT_FRACTION_PED
    injuries_for_model$baseline$whw$str_exponent_col[injuries_for_model$baseline$whw$cas_mode == "pedestrian"] <- SIN_EXPONENT_SUM_PED - (SIN_EXPONENT_SUM_PED * CASUALTY_EXPONENT_FRACTION_PED)
    
    injuries_for_model$baseline$nov$cas_exponent_col <- SIN_EXPONENT_SUM_NOV
    injuries_for_model$baseline$nov$str_exponent_col <- 1
    
    # when running in sampling mode, ensure that cas_exponent and str_exponent are always below 1
    injuries_for_model$baseline$whw$cas_exponent_col <- ifelse(injuries_for_model$baseline$whw$cas_exponent_col > 1, 1, injuries_for_model$baseline$whw$cas_exponent_col)
    injuries_for_model$baseline$whw$str_exponent_col <- ifelse(injuries_for_model$baseline$whw$str_exponent_col > 1, 1, injuries_for_model$baseline$whw$str_exponent_col)
    injuries_for_model$baseline$nov$cas_exponent_col <- ifelse(injuries_for_model$baseline$nov$cas_exponent_col > 1, 1, injuries_for_model$baseline$nov$cas_exponent_col)
    
    for (n in 1:(NSCEN + 1)) { # assign exponents to the scenario tables
      injuries_list[[n]]$whw$cas_exponent_col <- SIN_EXPONENT_SUM_VEH * CASUALTY_EXPONENT_FRACTION_VEH
      injuries_list[[n]]$whw$str_exponent_col <- SIN_EXPONENT_SUM_VEH - (SIN_EXPONENT_SUM_VEH * CASUALTY_EXPONENT_FRACTION_VEH)
      
      injuries_list[[n]]$whw$cas_exponent_col[injuries_list[[n]]$whw$cas_mode == "cycle"] <- SIN_EXPONENT_SUM_CYCLE * CASUALTY_EXPONENT_FRACTION_CYCLE
      injuries_list[[n]]$whw$str_exponent_col[injuries_list[[n]]$whw$cas_mode == "cycle"] <- SIN_EXPONENT_SUM_CYCLE - (SIN_EXPONENT_SUM_CYCLE * CASUALTY_EXPONENT_FRACTION_CYCLE)
      
      injuries_list[[n]]$whw$cas_exponent_col[injuries_list[[n]]$whw$cas_mode == "pedestrian"] <- SIN_EXPONENT_SUM_PED * CASUALTY_EXPONENT_FRACTION_PED
      injuries_list[[n]]$whw$str_exponent_col[injuries_list[[n]]$whw$cas_mode == "pedestrian"] <- SIN_EXPONENT_SUM_PED - (SIN_EXPONENT_SUM_PED * CASUALTY_EXPONENT_FRACTION_PED)
      
      # when running in sampling mode, ensure that cas_exponent and str_exponent are always below 1
      injuries_list[[n]]$whw$cas_exponent_col <- ifelse(injuries_list[[n]]$whw$cas_exponent_col > 1, 1, injuries_list[[n]]$whw$cas_exponent_col)
      injuries_list[[n]]$whw$str_exponent_col <- ifelse(injuries_list[[n]]$whw$str_exponent_col > 1, 1, injuries_list[[n]]$whw$str_exponent_col)
    }
    
    # when running in sampling mode, ensure that cas_exponent and str_exponent are always below 1
    injuries_for_model$baseline$nov$cas_exponent_col <- min(SIN_EXPONENT_SUM_NOV, 1)
    injuries_for_model$baseline$nov$str_exponent_col <- 1
    
    
    for (n in 1:(NSCEN + 1)) {
      # when running in sampling mode, ensure that cas_exponent and str_exponent are always below 1
      injuries_list[[n]]$nov$cas_exponent_col <- min(SIN_EXPONENT_SUM_NOV, 1)
      injuries_list[[n]]$nov$str_exponent_col <- 1
    }
  }
  
  
  ########################### Build regression model
  
  # run regression model on baseline data
  reg_model <- list()
  
  # define form of regression models for both whw and nov matrices
  forms <- list(
    whw = "count~cas_mode+strike_mode+offset(log(cas_distance)+(cas_exponent_col-1)*log(cas_distance_sum)+log(strike_distance)+
                    (str_exponent_col-1)*log(strike_distance_sum)+log(injury_reporting_rate)+log(weight))",
    nov = "count~cas_mode+offset(log(cas_distance)+(cas_exponent_col-1)*log(cas_distance_sum)+log(injury_reporting_rate)+log(weight))"
  )
  
  # add age and sex information if this information exists in the injuries_for_model data
  if ("age_cat" %in% names(injuries_for_model[[1]][[1]])) {
    for (type in INJURY_TABLE_TYPES) {
      forms[[type]] <- paste0(c(forms[[type]], "age_cat"), collapse = "+")
    }
  }
  if ("cas_gender" %in% names(injuries_for_model[[1]][[1]])) {
    for (type in INJURY_TABLE_TYPES) {
      forms[[type]] <- paste0(c(forms[[type]], "cas_gender"), collapse = "+")
    }
  }
  
  max_std <- list() # store largest standard errors
  max_std_def <- 10
  
  
  # try most sophisticated model first and then successively simplify model
  for (type in INJURY_TABLE_TYPES) {
    injuries_for_model[[1]][[type]]$injury_reporting_rate <- as.numeric(INJURY_REPORTING_RATE)
    
    # remove strike and cas mode pairs where cas mode = strike mode
    injuries_for_model[[1]][[type]]$cas_strike_mode <- apply(cbind(
      as.character(injuries_for_model[[1]][[type]]$cas_mode),
      as.character(injuries_for_model[[1]][[type]]$strike_mode)
    ), 1, function(x) paste(x, collapse = "_"))
    injuries_for_model[[1]][[type]] <- injuries_for_model[[1]][[type]] %>%
      filter(cas_strike_mode != "car_car" &
               cas_strike_mode != "bus_bus" &
               cas_strike_mode != "motorcycle_motorcycle" &
               cas_strike_mode != "cycle_cycle" &
               cas_strike_mode != "truck_truck")
    
    
    
    test <- "try-error"
    
    
    # try poisson distribution
    test <- try(glm(as.formula(forms[[type]]), data = injuries_for_model[[1]][[type]], family = "poisson"))
    
    if (length(test) != 1) {
      max_std[[type]] <- max(summary(test)$coefficients[, 2])
    } else {
      max_std[[type]] <- 1000
    }
    
    
    # if either standard errors are high or if no model could be built, try building model without age and gender categories if not already done
    if (length(test) == 1 | max_std[[type]] > max_std_def) {
      # remove age and gender if possible and fit model without age and gender categories
      if ("age_cat" %in% names(injuries_for_model[[1]][[1]])) {
        injuries_df <- injuries_for_model$baseline[[type]] # %>% dplyr::select(-c(age_cat, cas_gender))
        setDT(injuries_df)
        injuries_for_model$baseline[[type]] <- as.data.frame(injuries_df[, .(
          count = sum(count), weight = mean(weight), strike_distance_sum = mean(strike_distance_sum),
          cas_distance_sum = mean(cas_distance_sum),
          cas_exponent_col = mean(cas_exponent_col),
          str_exponent_col = mean(str_exponent_col)
        ), by = c("cas_mode", "strike_mode")])
        injuries_for_model$baseline[[type]]$strike_distance <- injuries_for_model$baseline[[type]]$strike_distance_sum
        injuries_for_model$baseline[[type]]$cas_distance <- injuries_for_model$baseline[[type]]$cas_distance_sum
        injuries_for_model[[1]][[type]]$injury_reporting_rate <- as.numeric(INJURY_REPORTING_RATE)
        
        # try model without age and gender categories
        forms_no_agecat <- list(
          whw = "count~cas_mode+strike_mode+offset(log(cas_distance)+(cas_exponent_col-1)*log(cas_distance_sum)
                                            +log(strike_distance)+(str_exponent_col-1)*log(strike_distance_sum)+
                                            log(injury_reporting_rate)+log(weight))",
          nov = "count~cas_mode+offset(log(cas_distance)+(cas_exponent_col-1)*log(cas_distance_sum)+log(injury_reporting_rate)+log(weight))"
        )
        
        # re-run regression model
        test <- try(glm(as.formula(forms_no_agecat[[type]]), data = injuries_for_model[[1]][[type]], family = "poisson"))
      }
    }
    
    # test whether standard error of newly built regression model is above cut-off distance and print
    # warning message if so
    if (length(test) != 1) {
      max_std[[type]] <- max(summary(test)$coefficients[, 2])
    } else {
      max_std[[type]] <- 0
    }
    
    if (max_std[[type]] > max_std_def) {
      print(paste0("!! The ", type, " injury model for ", city, " has large standard errors!!!"))
    }
    
    
    reg_model[[type]] <- test
    test <- NULL
  }
  
  browser()
  return(list(true_distances = true_distances, injuries_list = injuries_list, reg_model = reg_model, injuries_for_model = injuries_for_model))
}


dist_dur_tbls <- function(trips, scenarios) {
  
  #TODO: do we need to keep distances with trips? we can back out distances now
  WALK_SPEED <- 2.5
  BIKE_SPEED <- 7.2
  DRIVE_SPEED <- 13.8
  
  # fake it for now...
  trips$distance <- trips$minutes / 60 * DRIVE_SPEED
  trips[mode=="bike", distance := trips[mode=="bike", minutes] / 60 * BIKE_SPEED]  # in-place
  trips[mode=="walk", distance := trips[mode=="walk", minutes] / 60 * WALK_SPEED]
  trips[mode=="transit", distance := trips[mode=="transit", minutes]]
  
  # identify unique travel modes
  stage_modes <- unique(trips$mode)
  
  ## calculate all distances & durations for each scenario
  l_dist <- l_dur <- list()
  for (i in scenarios) { # loop through scenarios
    
    # get scenario trips
    local <- trips[trips$scenario == i]
    local_dist <- local[, .(sum_dist = sum(distance)), by = "mode"] # sum across distances by stage mode
    local_dur <- local[, .(sum_dur = sum(minutes)), by = "mode"] # sum across duration by stage mode
    
    # store results
    colnames(local_dist)[2] <- i
    l_dist[[i]] <- local_dist
    colnames(local_dur)[2] <- i
    l_dur[[i]] <- local_dur
  }
  
  ## join distances & durations
  for (i in names(l_dist)) {
    if (i == "baseline") {
      local_dist <- l_dist[[i]]
      local_dur <- l_dur[[i]]
    } else {
      local_dist <- local_dist[l_dist[[i]], on = "mode"] #<- left_join(local_dist, , by = "stage_mode")
      local_dur <- local_dur[l_dur[[i]], on = "mode"] # <- left_join(local_dur, l_dur[[i]], by = "stage_mode")
    }
  }
  
  # Remove walk to pt distances and durations as they were added to the pedestrian stages
  dist <- as.data.frame(local_dist[local_dist$mode != "walk_to_pt", ])
  dur <- as.data.frame(local_dur[local_dur$mode != "walk_to_pt", ])
  dist$stage_mode <- as.character(dist$mode)
  dur$stage_mode <- as.character(dur$mode)
  
  # return lsit of dist, duration dfs
  return(list(dist = dist, dur = dur, trips = trips))
}


get_all_distances <- function(trips, persons, scenarios) {
  
  # Generate distance and duration matrices taken from trip_scen_sets list as part of ithim_object
  # by calling the dist_dur_tbls.R function
  dist_dir_trips <- dist_dur_tbls(trips, scenarios)
  trips <- dist_dir_trips$trips
  keep <- c("PId", "scenario", "mode", "minutes", "distance")
  trips <- trips[, ..keep]
  #ithim_object$dist <- dist_and_dir$dist # create a distance matrix
  #ithim_object$dur <- dist_and_dir$dur # create a duration matrix
  #trip_scen_sets <- ithim_object$trip_scen_sets
  
  # grab baseline pop from persons table
  pop <- persons[scenario=="baseline"]
  
  # assign 5-year age bins
  pop$age_bin <- (pop$age %/% 5) * 5
  #pop <- pop %>% dplyr::rename(age_cat = age)
  pop_size <- nrow(pop)  # find the size of the total synthetic population
  
  # Recalculate distance by scaling to the total population (rather than the synthetic population)
  dist <- trips %>%
    group_by(mode, scenario) %>%
    summarise(
      ave_dist = sum(distance) / pop_size * nrow(pop)
    ) %>%
    spread(scenario, ave_dist)
  
  # add 'walk_to_pt' stage distances to 'pedestrian' stage distances
  if ("walk" %in% dist$mode && "walk_to_pt" %in% dist$mode) {
    dist[dist$mode == "walk", ][scenarios] <- dist[dist$mode == "pedestrian", ][scenarios] +
      dist[dist$mode == "walk_to_pt", ][scenarios]
    
    dist <- dist %>% filter(mode != "walk_to_pt")
  }
  
  ## Find population distances for each age and gender category by
  # - determining the proportion of age and sex specific mode distances
  #   to total mode distances in the synthetic population
  # - multiplying the total population distances by these proportions
  
  # join age bin, sex on trips
  join_cols <- c("PId", "age_bin", "sex")
  trips <- trips[pop[, ..join_cols], on="PId"]
  
  # find individual age / gender distances in synthetic population:
  trips_age_gender <- trips %>%
    group_by(age_bin, sex, mode, scenario) %>%
    summarise(dist_age = sum(distance))
  
  # find total trip distances by mode and scenario in the synthetic population
  trips_scen_mode <- trips %>%
    group_by(mode, scenario) %>%
    summarise(dist_total = sum(distance))
  
  # create data frame for each age, sex, mode and scenario combination containing the specific
  # distance travelled for this combination in the synthetic population and the total
  # distance travelled by the synthetic population for this specific mode and scenario
  trips_age_gender <- left_join(trips_age_gender, trips_scen_mode,
                                by = c("mode", "scenario")
  )
  
  # find proportion of total trip distance for each age and gender category distance in the synthetic population
  trips_age_gender$prop <- trips_age_gender$dist_age / trips_age_gender$dist_total
  
  # find total distance across entire population (and not just synthetic population)
  trips_age_gender$tot_dist <- trips_age_gender$dist_total / pop_size * nrow(pop)
  
  # scale total distance by trip proportion for each age and gender
  trips_age_gender$tot_dist <- trips_age_gender$tot_dist * trips_age_gender$prop
  journeys <- trips_age_gender %>% dplyr::select(-c(dist_age, dist_total, prop))
  
  # distances for injuries calculation but also parameterisation of Poisson injury regression model
  inj_distances <- distances_for_injury_function(
    journeys = journeys,
    dist = dist
  )
  
  return(list(dist, journeys, inj_distances))
}


set_injury_contingency <- function(trips, injuries) {
  
  # create list containing all modes in trip data set with speeds
  mode_names <- unique(trips$mode)
  
  injury_list <- list()
  injury_table_types <- c()
  
  # check whether there are any whw injuries given and create a whw matrix
  if (length(unique(injuries$strike_mode)) == 1 && !"nov" %in% injuries$strike_mode || length(unique(injuries$strike_mode)) > 1) {
    injury_list$whw <- subset(injuries, cas_mode %in% mode_names & strike_mode != "nov")
    injury_table_types <- c(injury_table_types, "whw")
  }
  # check if there are any nov injuries and create a nov matrix
  if ("nov" %in% injuries$strike_mode) {
    injury_list$nov <- subset(injuries, cas_mode %in% mode_names & strike_mode == "nov")
    injury_table_types <- c(injury_table_types, "nov")
  }
  
  injury_table <- list()
  
  # define column names to keep
  for (type in c(injury_table_types)) { # loop through 'whw' and 'nov'

    # Temporal injuries dataset
    temp_df <- injury_list[[type]]
    # Columns of interest
    group_columns <- c("cas_mode", "strike_mode", "age_bin", "sex")
    # Filter columns available in injuries dataset
    keep_names <- group_columns[group_columns %in% names(temp_df)]
    
    # summarise list of injuries by cas_mode, strike_mode, age_cat and cas_gender where this information exists
    # setDT(injury_list[[type]])
    injury_summary <- temp_df %>%
      group_by(across(keep_names)) %>%
      summarise(
        count = dplyr::n(),
        weight = mean(weight)
      ) %>%
      as.data.frame()
    
    # Conditional to restrict the number of injuries in dataset that don't have
    # cas_age nor cas_gender
    if (sum(names(injury_summary) %in% c("age_bin", "sex")) == 0) {
      # If no age nor gender exists, then each count is multiplied by the
      # proportion of deaths found in the GBD dataset
      injury_summary$count <- injury_summary$count * PROPORTION_INJURIES_AGERANGE
    }
    #
    
    ## create matrices for all cas and strike mode combinations (and all age and gender combinations) with 0 counts
    temp_table <- expand.grid(lapply(as.data.frame(temp_df)[, keep_names], unique))
    
    # Merge counts and weights to created matrix
    temp_table <- temp_table %>%
      left_join(injury_summary, by = group_columns)
    
    # Fill NAs because of no injuries in some cas and strike mode combinations
    temp_table$count[is.na(temp_table$count)] <- 0
    temp_table$weight[is.na(temp_table$weight)] <- mean(temp_df$weight)
    
    injury_table[[type]] <- temp_table
    rm(temp_df, temp_table)
  }
  
  # return injury tables in a list
  return(list(injury_table, injury_table_types))
}


prepareCrashData <- function(persons, trips, module_dir) {
  
  injuries <- fread(file.path(module_dir, "fars.csv"))
  
  # assign 5-year age bins
  injuries$age_bin <- (injuries$age %/% 5) * 5
  
  # if strike modeis null, replace with 'nov' (no vehicle)
  injuries[strike_mode=="", strike_mode := "nov"]  # in-place
  
  # Set weight as the unique number of years for which injury data exists
  injuries$weight <- length(unique(injuries$year))
  
  # Get all injuries where casualty and strike mode are identical for car, bus, motorcycle, cycle and truck
  # Treat bus_driver same as bus for strike mode
  same_cas_str_modes <- injuries %>% filter((cas_mode == "auto" & strike_mode == "auto") |
                                              (cas_mode == "bike" & strike_mode == "bike") |
                                              (cas_mode == "walk" & strike_mode == "walk"))
  
  # Filter all injuries where casualty equals strike mode
  # create dataset where casuality = strike mode fatality counts are removed
  injuries <- injuries %>% filter(!((cas_mode == "auto" & strike_mode == "auto") |
                                      (cas_mode == "bike" & strike_mode == "bike") |
                                      (cas_mode == "walk" & strike_mode == "walk")))
  
  # Where strike mode equals casualty mode, set the strike mode to 'nov'
  same_cas_str_modes <- same_cas_str_modes %>% mutate(strike_mode = "nov")
  
  # Join all injuries again
  injuries <- dplyr::bind_rows(injuries, same_cas_str_modes)
  
  # Call function to set tables for WHW and NOV
  injury_tables <- set_injury_contingency(trips, injuries)
  return(injury_tables)
}


injuries_function_2 <- function(true_distances, injuries_list, reg_model, constant_mode = F) {
  
  # create a list of all cas modes found within the whw and nov matrices
  cas_modes <- unique(as.character(injuries_list[[1]][[1]]$cas_mode))
  if (length(injuries_list[[1]]) == 2) {
    cas_modes <- unique(c(cas_modes, as.character(injuries_list[[1]]$nov$cas_mode)))
  }

  # create matrix containing all age categories for each sex
  demographic <- DEMOGRAPHIC
  demographic$dem_index <- 1:nrow(demographic)
  demographic <- demographic[, -which(colnames(demographic) == "population")]
  colnames(demographic)[which(colnames(demographic) == "age")] <- "age_cat"

  # join demographic with distance information to add demographic index
  injuries <- true_distances
  injuries <- dplyr::left_join(injuries, demographic, by = c("age_cat", "sex"))
  injuries$bus_driver <- 0
  injuries_lb <- injuries_ub <- injuries

  colnames(demographic)[which(colnames(demographic) == "sex")] <- "cas_gender"

  ############ predict fatality counts using the pre-defined Poisson regression model

  whw_temp <- list()
  for (scen in SCEN) { # loop through scenarios incl baseline
    whw_temp[[scen]] <- list()
    for (type in INJURY_TABLE_TYPES) { # loop through whw and nov matrices
      injuries_list[[scen]][[type]]$injury_reporting_rate <- 1

      # Set weight to 1 for prediction
      injuries_list[[scen]][[type]]$weight <- 1

      ## Use link type and use linkinv function to account for positive confidence interval for the poisson distr.
      # Ref: https://fromthebottomoftheheap.net/2018/12/10/confidence-intervals-for-glms/

      # grad the inverse link function
      ilink <- family(reg_model[[type]])$linkinv
      # add fit and se.fit on the **link** scale
      injuries_list[[scen]][[type]] <- bind_cols(
        injuries_list[[scen]][[type]],
        setNames(
          as_tibble(predict(reg_model[[type]],
            newdata = remove_missing_levels(
              reg_model[[type]],
              injuries_list[[scen]][[type]]
            ),
            type = "link", se.fit = TRUE
          )[1:2]),
          c("fit_link", "se_link")
        )
      )

      # create the 95% confidence interval using 2 times the standard error (se) and back-transform
      injuries_list[[scen]][[type]] <- mutate(injuries_list[[scen]][[type]],
        pred = ilink(fit_link),
        pred_ub = ilink(fit_link + (2 * se_link)),
        pred_lb = ilink(fit_link - (2 * se_link))
      )

      # for constant mode aggregate by strike mode and cas mode and create table with columns
      # containing the cas mode and rows for the strike modes for the whw and nov matrices
      if (constant_mode) {
        whw_temp[[scen]][[type]] <- sapply(unique(injuries_list[[scen]][[type]]$cas_mode), function(x) {
          sapply(
            unique(injuries_list[[scen]][[type]]$strike_mode),
            function(y) {
              sum(subset(
                injuries_list[[scen]][[type]],
                cas_mode == x & strike_mode == y
              )$pred, na.rm = T)
            }
          )
        })
        if (type == "whw") {
          colnames(whw_temp[[scen]][[type]]) <- unique(injuries_list[[scen]][[type]]$cas_mode)
          rownames(whw_temp[[scen]][[type]]) <- unique(injuries_list[[scen]][[type]]$strike_mode)
        } else {
          names(whw_temp[[scen]][[type]]) <- unique(injuries_list[[scen]][[type]]$cas_mode)
        }
      }

      # if constant_mode add additional tables where the predicted fatalities are aggregated by cas and strike mode for the
      # upper and lower 95% confidence interval limits
      if (constant_mode) {
        for (conf in c("ub", "lb")) { # loop through upper and lower confidence interval limits
          var_name <- paste0(type, "_", conf)
          whw_temp[[scen]][[var_name]] <- sapply(unique(injuries_list[[scen]][[type]]$cas_mode), function(x) {
            sapply(unique(injuries_list[[scen]][[type]]$strike_mode), function(y) {
              sum(subset(
                injuries_list[[scen]][[type]],
                cas_mode == x & strike_mode == y
              )[[paste0("pred_", conf)]], na.rm = T)
            })
          })
          if (type == "whw") {
            colnames(whw_temp[[scen]][[var_name]]) <- unique(injuries_list[[scen]][[type]]$cas_mode)
            rownames(whw_temp[[scen]][[var_name]]) <- unique(injuries_list[[scen]][[type]]$strike_mode)
          } else {
            names(whw_temp[[scen]][[var_name]]) <- unique(injuries_list[[scen]][[type]]$cas_mode)
          }
        }
      }

      # add demographic index to predicted fatality counts
      # dataframe and tibble, returns dataframe
      suppressWarnings(
        injuries_list[[scen]][[type]] <- dplyr::left_join(injuries_list[[scen]][[type]], demographic, by = c("age_cat", "cas_gender"))
      )
    }

    # set all columns values for the respective scenario and casualty mode to 0
    for (injured_mode in cas_modes) {
      for (index in unique(injuries$dem_index)) {
        injuries[
          injuries$scenario == scen & injuries$dem_index == index,
          match(injured_mode, colnames(injuries))
        ] <- 0
        injuries_lb[
          injuries_lb$scenario == scen & injuries_lb$dem_index == index,
          match(injured_mode, colnames(injuries_lb))
        ] <- 0
        injuries_ub[
          injuries_ub$scenario == scen & injuries_ub$dem_index == index,
          match(injured_mode, colnames(injuries_ub))
        ] <- 0
      }
    }

    # update values with new fatality predictions
    for (injured_mode in cas_modes) {
      for (index in unique(injuries$dem_index)) {
        for (type in INJURY_TABLE_TYPES) {
          injuries[injuries$scenario == scen & injuries$dem_index == index, match(injured_mode, colnames(injuries))] <-
            injuries[injuries$scenario == scen & injuries$dem_index == index, match(injured_mode, colnames(injuries))] +
            sum(injuries_list[[scen]][[type]][injuries_list[[scen]][[type]]$cas_mode == injured_mode &
              injuries_list[[scen]][[type]]$dem_index == index, ]$pred, na.rm = T) |> as.numeric()
        }
      }
    }

    # repeat for upper and lower confidence interval limits
    if (constant_mode) {
      for (injured_mode in cas_modes) {
        for (index in unique(injuries$dem_index)) {
          for (type in INJURY_TABLE_TYPES) {
            injuries_lb[injuries_lb$scenario == scen & injuries_lb$dem_index == index, match(injured_mode, colnames(injuries_lb))] <-
              injuries_lb[injuries_lb$scenario == scen & injuries_lb$dem_index == index, match(injured_mode, colnames(injuries_lb))] +
              sum(injuries_list[[scen]][[type]][injuries_list[[scen]][[type]]$cas_mode == injured_mode &
                injuries_list[[scen]][[type]]$dem_index == index, ]$pred_lb, na.rm = T) |> as.numeric()

            injuries_ub[injuries_ub$scenario == scen & injuries_ub$dem_index == index, match(injured_mode, colnames(injuries_ub))] <-
              injuries_ub[injuries_ub$scenario == scen & injuries_ub$dem_index == index, match(injured_mode, colnames(injuries_ub))] +
              sum(injuries_list[[scen]][[type]][injuries_list[[scen]][[type]]$cas_mode == injured_mode &
                injuries_list[[scen]][[type]]$dem_index == index, ]$pred_ub, na.rm = T) |> as.numeric()
          }
        }
      }
    }
  }


  # Create a total death count by summing across all casualty modes
  # Also remove NAs
  injuries <- injuries %>%
    ungroup() %>%
    mutate(Deaths = rowSums(dplyr::select(., cas_modes %>% as.character()), na.rm = T))

  # remove columns in injuries that still contain the original distances from the true_distances input
  # i.e remove columns where the mode is not in cas_mode
  to_keep <- c("age_cat", "sex", "scenario", "sex_age", "dem_index", levels(cas_modes), "Deaths")
  injuries2 <- injuries %>% dplyr::select(c(to_keep))

  # add lower and upper confidence interval death predictions
  if (constant_mode) {
    injuries_lb <- injuries_lb %>%
      ungroup() %>%
      mutate(Deaths_lb = rowSums(dplyr::select(., cas_modes %>% as.character()), na.rm = T))
    injuries_ub <- injuries_ub %>%
      ungroup() %>%
      mutate(Deaths_ub = rowSums(dplyr::select(., cas_modes %>% as.character()), na.rm = T))

    injuries2 <- dplyr::left_join(injuries2, injuries_lb %>% dplyr::select(age_cat, sex, dem_index, scenario, Deaths_lb),
      by = c("age_cat", "sex", "dem_index", "scenario")
    )
    injuries2 <- dplyr::left_join(injuries2, injuries_ub %>% dplyr::select(age_cat, sex, dem_index, scenario, Deaths_ub),
      by = c("age_cat", "sex", "dem_index", "scenario")
    )
  }

  list(injuries2, whw_temp)
}


# @title remove_missing_levels
#
# @description Accounts for missing factor levels present only in test data
# but not in train data by setting values to NA, i.e. if the data for which the predictions
# are made contains factor levels which do not appear in the baseline data used to
# parameterize the model, then we set the predictions for those factors to NA.
# Without this function, the entire model outputs would be NA if at least one factor level was unknown
#
# @import magrittr
# @importFrom gdata unmatrix
# @importFrom stringr str_split
#
# @param fit fitted model on training data, i.e. the baseline fatality counts in this case
#
# @param test_data data to make predictions for, i.e. injuries_list with cas, strike mode, age, sex and distance information for all scenarios
#
# @return data.frame with matching factor levels to fitted model
#
# @keywords internal
#
## !! temporary fix for missing (age) factors
#
# Adapted from  https://stackoverflow.com/a/39495480/4185785
#
#' @export


remove_missing_levels <- function(fit, test_data) {
  # drop empty factor levels in test data
  test_data <- as.data.frame(droplevels(test_data))

  # 'fit' object structure of 'lm' and 'glmmPQL' is different so we need to
  # account for it
  if (any(class(fit) == "glmmPQL")) {
    # Obtain factor predictors in the model and their levels
    factors <- (gsub(
      "[-^0-9]|as.factor|\\(|\\)", "",
      names(unlist(fit$contrasts))
    ))
    # do nothing if no factors are present
    if (length(factors) == 0) {
      return(test_data)
    }

    map(fit$contrasts, function(x) names(unmatrix(x))) %>%
      unlist() -> factor_levels
    factor_levels %>%
      str_split(":", simplify = TRUE) %>%
      extract(, 1) -> factor_levels

    model_factors <- as.data.frame(cbind(factors, factor_levels))
  } else {
    # Obtain factor predictors in the model and their levels
    factors <- (gsub(
      "[-^0-9]|as.factor|\\(|\\)", "",
      names(unlist(fit$xlevels))
    ))
    # do nothing if no factors are present
    if (length(factors) == 0) {
      return(test_data)
    }

    factor_levels <- unname(unlist(fit$xlevels))
    model_factors <- as.data.frame(cbind(factors, factor_levels))
  }

  # Select column names in test data that are factor predictors in
  # trained model

  predictors <- names(test_data[names(test_data) %in% factors])

  # For each factor predictor in your data, if the level is not in the model,
  # set the value to NA

  for (i in 1:length(predictors)) {
    found <- test_data[, predictors[i]] %in% model_factors[
      model_factors$factors == predictors[i],
    ]$factor_levels
    if (any(!found)) {
      # track which variable
      var <- predictors[i]
      # set to NA
      test_data[!found, predictors[i]] <- NA
      # drop empty factor levels in test data
      test_data <- droplevels(test_data)
      # issue warning to console
      message(sprintf(
        paste0(
          "Setting missing levels in '%s', only present",
          " in test data but missing in train data,",
          " to 'NA'."
        ),
        var
      ))
    }
  }
  return(test_data)
}
