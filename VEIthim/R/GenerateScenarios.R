# this module appends counterfactual scenario to persons table

appendCounterfactualScenario <- function(persons){
  
  # establish baseline scenario
  persons$scenario <- "baseline"
  
  # for now, we'll just make up a counterfactual
  # that everyone drives more, walks/bikes/takes transit less
  # TODO: this will eventually reference other VE outputs
  counterfactual <- copy(persons)
  counterfactual$scenario <- "test_scenario"
  counterfactual$auto <- counterfactual$auto * 1.2
  counterfactual$transit <- counterfactual$transit * 0.8
  counterfactual$walk <- counterfactual$walk * 0.8
  counterfactual$bike <- counterfactual$bike * 0.8
  
  # bind two scenarios
  persons <- rbind(persons, counterfactual)
  
  # easiest way to integrate with ITHIM is to treat each mode as a single trip
  # so will convert persons into a quasi trips table
  trip_cols <- c("PId", "auto", "transit", "walk", "bike", "leisurePA",
                 "body_mass", "ecf", "rmr", "vo2max",
                 "intercept_a", "slope_b", "sd_test_level", "d_k"
                 )
  trips <- melt(persons[, ..trip_cols],
                measure.vars = c("auto", "transit", "walk", "bike"),
                variable.name = "mode", value.name = "minutes")[order(PId)]
  
  return(list(trips, persons))
}

  
