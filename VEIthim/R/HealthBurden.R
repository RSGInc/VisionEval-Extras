#' Compute health burden
#'
#' Compute health burden for population in scenarios given relative risks for diseases
#'
#' This function performs the following steps:
#'
#' \itemize{
#'  \item get the demographic and disease burden (subset of Global Burden of Disease dataset) data
#'    into the correct formats and join the two datasets
#'
#'  \item scale the burden data by the CHRONIC_DISEASE_SCALAR to account for bias in the data
#'
#'  \item split the above dataframe into two dataframes, one for deaths and one
#'    for years of life lost (YLLs)
#'
#'  \item add a demographic index (by age and sex category) to the dataframe containing the
#'    individual relative risk for different diseases
#'
#'  \item set the reference and the other scenarios
#'
#'  \item iterate over all disease outcomes:
#'    \itemize{
#'    \item define column names
#'
#'    \item loop either over 1 or 2 pathways depending on whether both PA and AP are affecting
#'      the disease and whether the AP and PA pathways are combined or not:
#'      \itemize{
#'      \item extract the relevant burden of disease for the specific scenario both for YLLs and deaths
#'
#'      \item find the sum of the relative risks (RR) for the specific disease for each age
#'        and sex category for the reference scenario
#'
#'      \item loop through the non-reference scenarios:
#'        \itemize{
#'        \item define column names
#'
#'        \item find the sum of the relative risks (RR) for the specific disease for each age
#'          and sex category for the non-reference scenario
#'
#'        \item calculate the PIF (potential impact fraction), i.e the proportional change
#'          in the sum of relative risks between the reference and the non-reference
#'          scenario for each age and sex category
#'
#'        \item calculate the health burden (deaths and ylls) for the non-reference scenario
#'          compared to the reference scenario by multiplying the current burden of disease
#'          by the PIF (combine_health_and_pif.R)
#'          }
#'
#'      \item if confidence intervals are required, loop through the upper and lower confidence
#'        interval limits and calculate the health burden for deaths and YLLs using
#'        the upper and lower confidence relative risks. If no upper and lower relative risk
#'        values exist, use the median value instead
#'        }
#'    }
#' }
#'
#'
#' @param ind_ap_pa dataframe of all individuals' relative risks for diseases
#' @param conf_int=F logic: whether to include confidence interval from dose response relationships or not
#' @param combined_AP_PA=T logic: whether to combine the two exposure pathways (AP and PA) or to compute them independently
#'
#' @return list of dataframes: one for deaths per disease per demographic group and scenario, and likewise for YLLs
#'
#' @export


health_burden <- function(ind_ap_pa, module_dir, input_dir, scenarios, conf_int = F, combined_AP_PA = T) {
  
  # load all relative risk values
  # these are stored in global env in ITHIM, which is not ideal....
  diseases <- fread(file.path(module_dir, "disease_outcomes_lookup.csv"))
  
  # load gbd data for region
  gbd_data <- fread(file.path(input_dir, "ithim", "gbd_oregon.csv"))
  
  # convert gbd data to rates
  gbd_data$deathrate <- gbd_data$deaths / gbd_data$ref_pop
  gbd_data$yllrate <- gbd_data$ylls / gbd_data$ref_pop
  
  # initiate dfs to store yll, deaths
  ind_ap_pa <- ind_ap_pa[order(PId)]
  cols <- c("PId")
  ylls <- deaths <- ind_ap_pa[, ..cols]
  
  ### iterate over all disease outcomes
  for (j in 1:nrow(diseases)) {
    
    # Disease acronym and full name
    ac <- as.character(diseases$acronym[j])
    gbd_dn <- as.character(diseases$GBD_name[j])
    
    # calculate health outcome, or independent pathways?
    pathways_to_calculate <- ifelse(combined_AP_PA, 1,
                                    diseases$physical_activity[j] + diseases$air_pollution[j]
    )
    
    # loop through relevant pathways depending on whether disease affects AP or PA and
    # whether AP and PA are combined or not
    for (path in 1:pathways_to_calculate) { # pathways_to_calculate is either 1 or 2
      # set up column names
      if (combined_AP_PA) { # if combined, find whether only pa or ap or both
        middle_bit <-
          paste0(
            ifelse(diseases$physical_activity[j] == 1, "pa_", ""),
            ifelse(diseases$air_pollution[j] == 1, "ap_", "")
          )
      } else {
        # if independent, choose which one
        middle_bit <- c("pa_", "ap_")[which(c(diseases$physical_activity[j], diseases$air_pollution[j]) == 1)[path]]
      }
      
      # set column names
      base_var <- paste0("RR_", middle_bit, scenarios[1], "_", ac)
      scen_var <- paste0("RR_", middle_bit, scenarios[2], "_", ac)
      
      # subset gbd data for deaths and ylls by the given disease
      gbd_disease <- gbd_data[cause == ac]
      
      # isolate rr for pathway
      person_cols <- c("PId", "age_bin", "sex", base_var, scen_var)
      disease_risk <- ind_ap_pa[, ..person_cols]
      
      # join the population, gbd data
      join_cols <- c("age_bin", "sex", "deathrate", "yllrate")
      disease_risk <- disease_risk[gbd_disease[, ..join_cols],
                                   on = .(age_bin, sex),
                                   nomatch = NULL]
      disease_risk <- disease_risk[order(PId)]
      
      ## set up pif tables
      # extract the relative baseline risk and demographic index
      #pif_table <- setDT(ind_ap_pa[, colnames(ind_ap_pa) %in% c(base_var, "dem_index")])
      # set the relative risk column name to 'outcome'
      #setnames(pif_table, base_var, "outcome")
      # sum the outcomes for the synthetic population by age and sex category
      #pif_ref <- pif_table[, .(sum(outcome)), by = "dem_index"]
      ## sort pif_ref
      #setorder(pif_ref, dem_index) # order by age and sex category, i.e dem_index
        
      # set up naming conventions
      scen <- scenarios[2]
      yll_name <- paste0(scen, "_ylls_", middle_bit, ac)
      deaths_name <- paste0(scen, "_deaths_", middle_bit, ac)
      
      # calculate PIF (i.e. change in RR compared to reference scenario) for this scenario and convert to vector
      pif_scen <- (disease_risk[, ..base_var] - disease_risk[, ..scen_var]) / disease_risk[, ..base_var]
      disease_risk$pif <- (
        (disease_risk[, ..base_var] - disease_risk[, ..scen_var]) / disease_risk[, ..base_var])
        
      # Calculate the difference in ylls between the non-reference and the reference scenario
      # by multiplying the current burden of disease for particular disease by the PIF value,
      # i.e the expected change between the reference scenario and the given scenario
      # for each age and sex category
      disease_ylls <- disease_risk$yllrate * disease_risk$pif
      ylls[[yll_name]] <- disease_ylls
        
      # Calculate the difference in deaths between the non-reference and the reference scenario
      # by multiplying current burden of disease for particular disease by the PIF value,
      # i.e the expected change between the reference scenario and the given scenario
      # for each age and sex category
      disease_deaths <- disease_risk$deathrate * disease_risk$pif
      deaths[[deaths_name]] <- disease_deaths
      
    } # end of pathways loop
  } # end of disease loop
  
  # return list of (deaths, ylls)
  return(list(deaths = deaths, ylls = ylls))
}



#' Join disease health burden and injury data
#'
#' Join the two data frames for health burden: that from disease, and that from road-traffic injury
#'
#' This function performs the following steps:
#'
#' \itemize{
#'  \item extract the yll and deaths data from the AP and PA pathways
#'  \item extract the yll and deaths data from the injury data
#'  \item create one dataframe for yll and one for deaths containing all the AP, PA and injury data
#' }
#'
#' @param ind_ap_pa list (deaths, YLLs) of data frames of all demographic groups' burdens of diseases
#' @param inj list (deaths, YLLs) of data frames of all demographic groups' burdens for road-traffic injury
#'
#' @return list of dataframes: one for deaths per cause per demographic group, and likewise for YLLs
#'
#' @export


join_hb_and_injury <- function(ind_ap_pa, inj) {
  deaths <- ind_ap_pa$deaths
  ylls <- ind_ap_pa$ylls
  
  # Select deaths columns from injury data
  inj_deaths <- dplyr::select(inj, c(age_cat, sex, contains("deaths")))
  
  # Select yll columns from injury data
  inj_ylls <- dplyr::select(inj, c(age_cat, sex, contains("yll")))
  
  # Join injuries data to global deaths and yll datasets
  deaths <- dplyr::left_join(deaths, inj_deaths, by = c("age_cat", "sex"))
  ylls <- dplyr::left_join(ylls, inj_ylls, by = c("age_cat", "sex"))
  list(deaths = deaths, ylls = ylls)
}


combine_health_and_pif <- function(pif_values, hc = DISEASE_BURDEN) {
  setorder(hc, dem_index)
  hm_cn_values <- hc$burden
  return_values <- hm_cn_values * pif_values # multiply burden of disease times PIF
  round(as.vector(return_values), 5)
}
