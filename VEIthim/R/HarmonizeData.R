# this module prepares VE outputs for ITHIM

estimateVentilationRates <- function(persons, input_dir){
  
  # also need to estimate body mass for inhalation rates
  # we will draw from distributions supplied in core ITHIM codebase
  input_dir <- file.path(input_dir, "ithim")
  
  # Body mass
  body_mass <- fread(file.path(
    input_dir, "BodyMass.csv"),
    select = c("age", "sex", "gm", "gsd", "lower", "upper"))
  
  # Energy Conversion Factor
  ecf <- fread(file.path(
    input_dir, "ECF.csv"),
    select = c("age", "sex", "lower", "upper"))
  
  # Resting Metabolic Rate
  rmr <- fread(file.path(
    input_dir,"RMR.csv"),
    select = c("age", "sex", "slope", "interc", "se"))
  
  # maximum oxygen uptake rate (NVO2max)
  nvo2max <- fread(file.path(
    input_dir, "NVO2max.csv"),
    select = c("age", "sex", "mean", "sd", "lower", "upper"))
  
  # Ventilation from o2 Uptake 
  vent_from_oxygen <- fread(file.path(
    input_dir, "VentilationFromOxygenUptake.csv"),
    select = c("age", "sex", "intercept_a", "slope_b", "sd_person_level", "sd_test_level"))  
  
  # Dan: Draw from a log-normal distribution the body mass [kg] of each person
  # in the synthetic population
  body_mass <- persons[body_mass, on = c('age', 'sex'),
                       nomatch = NULL][order(PId)] %>%
    rowwise() %>% # To apply an operation in each row
    mutate(
      # Draw sample from distribution
      sample = rlnorm(1, log(gm), log(gsd)),
      # Check maximum and minimum values
      body_mass = ifelse(sample < lower, lower,
                         ifelse(sample > upper, upper, sample))) %>% 
    dplyr::select(PId, age, sex, body_mass)
  
  # Dan: Draw from a uniform distribution the Energy Conversion Factor [lt/kcal]
  # of each person in the synthetic population
  body_mass <- body_mass %>%
    left_join(ecf, by = c("age", "sex")) %>%
    rowwise() %>% # To apply an operation in each row
    mutate(
      # Draw sample from distribution
      ecf = runif(1, min = lower, max = upper),
    ) %>%
    dplyr::select(PId, age, sex, body_mass, ecf)
  
  # Dan: For each person in the synthetic population, calculate the Resting
  # Metabolic Rate [kcal/min] from a given regression formula with a normally
  # distributed error
  body_mass <- body_mass %>%
    left_join(rmr, by = c("age", "sex")) %>%
    rowwise() %>% # To apply an operation in each row
    mutate(
      # Draw sample from distribution
      error = rnorm(1, 0, se),
      # Calculate RMR
      rmr = 0.166 * (interc + (slope * body_mass) + error)
    ) %>%
    dplyr::select(PId, age, sex, body_mass, ecf, rmr)
  
  # Dan: Draw from a normalized maximum oxygen uptake rate (NVO2max) [lt/(min*kg)]
  # of each person in the synthetic population
  body_mass <- body_mass %>%
    left_join(nvo2max, by = c("age", "sex")) %>%
    rowwise() %>% # To apply an operation in each row
    mutate(
      # Draw sample from distribution
      sample = rnorm(1, mean, sd),
      # Check maximum and minimum values
      sample_check = ifelse(sample < lower, lower,
                            ifelse(sample > upper, upper, sample)
      ),
      # Transform [ml/(min*kg)] to [lt/(min*kg)]
      nvo2max = sample_check / 1000
    ) %>%
    dplyr::select(PId, age, sex, body_mass, ecf, rmr, nvo2max)
  
  # Dan: Calculate the maximum oxygen uptake rate [lt/min] of each person in the
  # synthetic population
  body_mass <- body_mass %>% mutate(vo2max = nvo2max * body_mass)
  
  # Dan: Get the parameters for the empirical equation to calculate the
  # Ventilation Rate from Oxygen Uptake Rate in each person in the
  # synthetic population
  body_mass <- body_mass %>%
    left_join(vent_from_oxygen, by = c("age", "sex")) %>%
    rowwise() %>% # To apply an operation in each row
    mutate(
      # Draw sample from normal distribution taking into account the variance
      # between persons
      d_k = rnorm(1, 0, sd_person_level)
    ) %>%
    dplyr::select(
      PId, body_mass, ecf, rmr, nvo2max, vo2max,
      intercept_a, slope_b, sd_person_level, sd_test_level, d_k
    )
  
  return(body_mass)
}


hhsToPersons <- function(Hhs, genderRatio, input_dir){
  
  WALK_SPEED <- 2.5
  BIKE_SPEED <- 7.2
  DRIVE_SPEED <- 13.8
  
  # first, rename mode cols
  setnames(Hhs,
           c("Dvmt", "TransitTrips", "WalkPMT", "BikePMT"),
           c("auto", "transit", "walk", "bike"))
  
  # and convert distances to duration
  Hhs$auto <- Hhs$auto / DRIVE_SPEED * 60
  Hhs$walk <- Hhs$walk / WALK_SPEED * 60
  Hhs$bike <- Hhs$bike / BIKE_SPEED * 60
  
  # melt Hh table so that each age bin is a row
  # replace 0 with NA onthe fly to so we can use na.rm in melt()
  persons <- melt(
    replace(Hhs, Hhs==0, NA),  # so we can use na.rm = TRUE
    id.vars=c("HhId", "Bzone"),
    measure = patterns("Age"),
    variable.name = "ageBin",
    value.name = "nPersons",
    na.rm = TRUE)
  
  # we need to join a few items back to persons from hh table
  join_cols <- c("HhId", "HhSize", "auto", "transit", "walk", "bike")
  persons <- persons[Hhs[, ..join_cols], on = .(HhId)]
  
  # now, each person becomes a row
  persons <- data.table(
    as.data.frame(lapply(persons, rep, persons$nPersons)))
  persons[,nPersons:=NULL]  # in-place
  
  # randomly assign each person an age in the relevant bin
  # first need to parse min, max age within each age bin
  persons$ageBin <- substring(persons$ageBin, 4)
  persons[, c("minAge", "maxAge") := tstrsplit(ageBin, "to", fixed=TRUE)]
  persons[minAge == "65Plus", minAge := "65"]  # in-place
  persons[is.na(maxAge), maxAge := "85"]  # in-place
  
  # TODO: this is a little slow
  persons$age <- mapply(
    function(x, y) sample(seq(x, y), 1), 
    persons$minAge, persons$maxAge)
  
  # next, assign gender to each person based on global gender ratio
  # for now, we're not going to worry about gender distribution within Hhs
  # in the future, we should do this more carefully
  values = c('male', 'female')
  weights = c(genderRatio, 1 - genderRatio)
  persons$sex <- sample(
    values,
    size=nrow(persons),
    replace = TRUE,
    prob = weights)
  
  # assign person uid so that we can keep track of persons across scenarios
  persons$PId <- 1:nrow(persons)
  
  # estimate ventilation rates and join to persons
  body_mass <- estimateVentilationRates(persons, input_dir)
  persons <- persons[setDT(body_mass), on = c('PId'), nomatch = NULL][order(PId)]
  
  # finally, decompose Hh travel totals to individuals
  # for now, just distribute equally across Hh members
  # in the future, we should do this more carefully
  travel_cols = c("auto", "transit", "walk", "bike")
  for (col in travel_cols) {
    replacement_values <- persons[, ..col] / persons$HhSize
    persons[, (col) := replacement_values]
  }
  
  # clean up schema and return persons df
  scrub <- c("ageBin", "HhSize", "minAge", "maxAge")
  persons[ ,(scrub) := NULL]  # in-place
  return(persons)
}


appendZoneAttributesToPersons <- function(Persons, zonePM, zonePA, year){
  
  # join background PM2.5 to Persons table
  # from Bzones
  join_cols <- c("Bzone", "PM25")
  Persons <- Persons[zonePM[zonePM$Year == year][, ..join_cols],
                     on = .(Bzone),
                     nomatch = NULL]
  
  # join background PA to Persons table
  # from Bzones
  join_cols <- c("Bzone", "pa0", "pa1", "pa2", "pa3")
  Persons <- Persons[zonePA[zonePA$Year == year][, ..join_cols],
                     on = .(Bzone),
                     nomatch = NULL]
  
  # background PA is a distribution, so assign a value to each person
  prob_samp <- function(x) sample.int(n=4, size=1, prob=x) # from jay.sf
  Persons[, pa_group := apply(.SD, 1, prob_samp), .SDcols=c("pa0", "pa1", "pa2", "pa3")]
  
  # now that we have a bin for each person, assign min/max PA values
  Persons$minPA <- 0
  Persons$maxPA <- 39
  Persons[pa_group == 2, c("minPA", "maxPA") := list(40, 74)]  # in-place
  Persons[pa_group == 3, c("minPA", "maxPA") := list(75, 149)]
  Persons[pa_group == 4, c("minPA", "maxPA") := list(150, 225)]

  # TODO: this is a little slow
  Persons$leisurePA <- mapply(
    function(x, y) sample(seq(x, y), 1), 
    Persons$minPA, Persons$maxPA)
  
  # clean up schema and return persons df
  scrub <- c("pa0", "pa1", "pa2", "pa3", "pa_group", "minPA", "maxPA")
  Persons[ ,(scrub) := NULL]  # in-place
  return(Persons)
}
  
