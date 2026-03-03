# Model
# Record of Revision:
#   23-06-2025 Kelsi
#   03-07-2025 Kelsi
#   04-07-2025 Kelsi
#   07-07-2025 Kelsi
#   10-07-2025 Kelsi 
#   23-07-2025 Tara 
#   24-07-2025 Kelsi
#   25-07-2025 Kelsi
#   27-11-2025 Tara - DOC correction 
#   17-12-2025 Tara - set seed

##### Packages & Data ---------------------------------------------------------

#setwd("C:/Users/P079516/OneDrive - Amsterdam UMC/One View/Documenten/PATTERN2/Study 4- PATTERN II model/Model/")# for tara
setwd("C:/Users/P077952/Amsterdam UMC/Boute, T.C. (Tara) - Study 4- PATTERN II model/Model") # for kelsi

# packages
library(tidyverse)
library(readxl)
library(ggplot2)
library(survival)
library(flexsurv)
library(flexsurvcure)
library(Rlab)

##### The Model -------------------
# ***REMEMBER*** cured_status: 0 = cured, 1 = non-cured



fup_model_fast <- function(data, sensitivities, specificities = NULL,
                           ct_freq, cea_freq_min, cea_freq_max, colo_freq, 
                           ct_freq2, cea_freq2_min, cea_freq2_max, colo_freq2, 
                           ct_switch, cea_switch, colo_switch, 
                           ct_stop=60, cea_stop=60, colo_stop=60,
                           adh_ct, adh_cea, adh_colo, 
                           DOD_shape, DOD_scale){
  # set.seed(06092025) # for reproducible results (so 90 day mortality stays the same)
  
  #source("2. Informing parameters.R") # do not run everytime, quicker to just load final models that we need
  # Loads 4 models:
  # 1- mor90D_1: model for 90-day mortality (gompertz model)
  # 2- DOC1: model for death other causes (gompertz model)
  # 3- TTR_cure_model: model for time to recurrence based on topo and pt/pn high-risk (lognormal cure model)
  # 4- symp_cure_model: model for symptomatic detection time / delayed TTR
  load("model_inputs.RData")
  
  n_patients <- nrow(data)
  
  ###### Event times: generate 6 times (always using months except for the 90 day mortality) -------------------------------------------------
  #   1 - 90 day mortality time
  #   2 - death other causes (DOC) time
  #   3 - time detectable / time until recurrence - based on risk profile 
  #   4 - symptomatic time 
  #   5 - follow-up (FUP) detection time
  #   6 - death of disease (DOD) time
  
  # 90-day mortality time (can add seed = XXX) to keep it reproducible  ###, seed = 1234
  set.seed(1234)
  data$mor90d_time <- simulate(mor90D_1, newdata=data.frame(age_cat2=data$age_cat))$time/30.4167 # sample event times from Gompertz distribution for 90D mortality (convert from days to months)
  data$mor90d_time[data$mor90d_time > 3] <- Inf
  
  # Death other causes time
  #set.seed(1234)
  data$DOC_time <- simulate(DOC1, newdata=data.frame(age_cat2=data$age_cat))$time + 3  # sample event times from Gompertz distribution for death other causes 
  # + 3 months correction 90DM 
  
  # Detectable time 
  #set.seed(1234)
  mean_lnorm <- 3.3110  -0.2118*(data$topo == 2) - 0.5136*(data$high_risk == 1) # using parameter from cure model with effects for topo/high_risk on both theta and lognormal distribution mean
  recurrence_time <- rlnorm(n_patients, mean = mean_lnorm, sd = 1.0745) # time until recurrence based on topo and pt/pn high-risk 
  data$recurrence_time <- recurrence_time
  
  data$recurrence_time[data$recurrence_time > 60] <- Inf
  data$cured_status[data$recurrence_time > 60] <- 0
  
  data$detectable_time <- ifelse(data$cured_status == 1, # if cured status is not cured (0=cured, 1=not-cured)
                                 recurrence_time - rnorm(n_patients, 6, 1), # then find time when detectable-- this will change to be dependent on disease specific things 
                                 Inf) # otherwise set at Inf months (i.e., cured patients never get a recurrence so never detectable again)
  data$detectable_time[data$detectable_time < 0] <- 0 
  
  
  
  # Simulate symptomatic detection time / delayed TTR
  #set.seed(1234)
  log_HR_samples <- rnorm(n_patients, mean = symp_cure_model$res["FUP_pool1", "est"], sd = symp_cure_model$res["FUP_pool1", "se"])  # Introduce uncertainty: generate 1000 samples of log-HR using SE directly
  lambda_samples <- 1 / exp(log_HR_samples) # Convert to hazard ratio (HR) and Compute lambda for each simulated HR
  lambda_patient <- sample(lambda_samples, size = n_patients, replace = TRUE) # Trek per patiënt een unieke lambda uit de gesimuleerde distributie
  data$lambda_patient <- lambda_patient
  data$symprec_times <- ifelse(data$cured_status == 1, # simulate symptomatic time for patients that were not cured (0=cured, 1=not-cured)
                               lambda_patient*recurrence_time, # scale the recurrence time by the lambda 
                               Inf) # otherwise set at Inf months (i.e., cured patients never get a recurrence so never have a symptomatic time)
  
  data$symprec_times[data$mor90d_time < data$symprec_times | data$DOC_time < data$symprec_times] <- Inf # patients that died before symptom detected time cannot get a symptomatic detection, so set to Inf 
 # data$symprec_times[data$symprec_times > 120] <- Inf
  
  
  ######  Follow-up schedule:  -------------------------------------------------
  # simulate test visits and results according to frequency, adherence and sensitivity of the test
  simulate_test <- function(freq, freq2, switch, adherence, sens, data, stoptime, spec) {
    ## Precompute test schedules for each unique frequency
    unique_freqs1 <- unique(freq) # see how many unique values of freq1 there are (before switch)
    unique_freqs2 <- unique(freq2) # do the same for freq2 (after switch)
    unique_freqs <- expand.grid(unique_freqs1,unique_freqs2) # make all combinations of freq before and after switch
    
    sched_map  <- lapply(1:nrow(unique_freqs), function(i) c(seq(unique_freqs[i,1], switch, by = unique_freqs[i,1]),
                                                             seq(switch + unique_freqs[i,2], stoptime, by = unique_freqs[i,2])))
    names(sched_map) <- paste0(unique_freqs[,1], ",", unique_freqs[,2])
    
    schedule_list <- sched_map[ paste0(freq, ",", freq2)]
    
    if (length(freq)==1 & length(freq2)==1){
      freq_new <- rep(as.character(freq), n_patients)
      freq2_new <- rep(as.character(freq2), n_patients)
      schedule_list <- sched_map[paste0(freq_new, ",", freq2_new)]
    }
    
    test_list <- vector("list", n_patients)
    result_list <- vector("list", n_patients)
    sens_vec <- sens[data$loc_rec] # loc_num is coded 0-4 so we add 1 so that it takes position 1-5 from the sensitivites  
    
    if (sum(as.numeric(freq))==0) return(list(tests = test_list, results = result_list))
    
    for (i in seq_len(n_patients)) {
      
      test_months <- schedule_list[[i]] # months that they COULD have a visit 
      
      if (length(test_months) == 0) { # no schedule at all
        test_list[[i]]   <- integer(0)
        result_list[[i]] <- integer(0)
        next
      }
      
      attend <- rbinom(length(test_months), 1, adherence) # person attends ct scan with probability adher http://127.0.0.1:22397/graphics/757192f1-57ee-4f4b-af6c-407c6379996d.pngence, 1=attend, otherwise 0
      test_list[[i]] <- test_months[attend == 1] # select months that the patient actually attended based on above probability 
      
      if (data$cured_status[i] == 0) {  #0 = cured, 1 = non-cured
        # if they were cured but they showed up for a test then that test will always be negative
        result_list[[i]] <- rep(0, length(test_list[[i]]))
        next
      } 
      
      test_visits <- test_list[[i]] # times patient showed up for a test visit
      if (length(test_visits) == 0) {
        result_list[[i]] <- integer(0) # if they never showed up then they never got a positive test
        next
      }
      
      res <- rbinom(length(test_visits), 1, sens_vec[i]) # generate results at visit i
      
      detect_time <- data$detectable_time[i] # detectable time for this patient 
      if (is.null(spec)){
        res[test_visits < detect_time] <- 0      # force zero (i.e. negative test results) before detectable time
        
      }else{# code false positives with an 8
        # before detectability time it is possible that a false positive is recorded, then an 8 will be placed in the results 
        n_before_detect <- sum(test_visits < detect_time)
        res[test_visits < detect_time] <- ifelse(rbinom(n_before_detect, 1, 1 - spec[data$loc_rec][i])==1, 8, 0)
      }
      
      # apply stopping after first positive (remove visits and results after first positive)
      pos_idx <- which(res == 1)
      if (length(pos_idx) > 0) {
        first_month <- test_visits[min(pos_idx)]
        keep <- test_visits <= first_month
        test_list[[i]]   <- test_visits[keep]
        result_list[[i]] <- res[keep]
      } else {
        result_list[[i]] <- res
      }
      
      # check if patients has already died from death other causes (DOC) or 90-day mortality, and then remove future tests: 
      # --> if a patient dies (also with DOC), a patient does not receive any tests anymore
      death_month <- floor(min(data$mor90d_time[i], data$DOC_time[i])) # month of death (earliest of DOC or 90d death)
      
      # only do something if death occurred before the maximum possible time
      if (death_month < 60) {
        keep_idx <- test_list[[i]] <= death_month
        test_list[[i]]   <- test_list[[i]][keep_idx]
        result_list[[i]] <- result_list[[i]][keep_idx]
      }
    }
    return( list(tests = test_list, results = result_list))
  }
  
  # CEA tests can have a range so for each person randomly pick a number of visits 
  if (cea_freq_min == cea_freq_max) {
    cea_freq <-  12/cea_freq_max
  }else{
    cea_freq <- as.numeric(12/sample(seq(cea_freq_min, cea_freq_max, 1), n_patients, replace=T))
  }

  if (cea_freq2_min == cea_freq2_max) { #after the switch
    cea_freq2 <-  12/cea_freq2_max
  }else{
    cea_freq2 <- as.numeric(12/sample(seq(cea_freq2_min, cea_freq2_max, 1), n_patients, replace=T))

  }

  if (is.null(specificities)){
    # make matrices for test visits of results at each visit and timing of each visit 
    ct_res <- simulate_test(ct_freq, ct_freq2, ct_switch, adh_ct, sensitivities$CT, data, ct_stop)
    cea_res <- simulate_test(cea_freq, cea_freq2, cea_switch, adh_cea, sensitivities$CEA, data, cea_stop)
    colo_res <- simulate_test(colo_freq, colo_freq2, colo_switch, adh_colo, sensitivities$Coloscopie, data, colo_stop)
    
  }else{
    # specificity incorporation
    # make matrices for test visits of results at each visit and timing of each visit 
    ct_res <- simulate_test(ct_freq, ct_freq2, ct_switch, adh_ct, sensitivities$CT, data, ct_stop, specificities$CT)
    cea_res <- simulate_test(cea_freq, cea_freq2, cea_switch, adh_cea, sensitivities$CEA, data, cea_stop, specificities$CEA)
    colo_res <- simulate_test(colo_freq, colo_freq2, colo_switch, adh_colo, sensitivities$Coloscopie, data, colo_stop, specificities$Coloscopie)
    
  }
  ### home based testing: error with test 
  
  # make new columns and fill with NA to start with
  data$detection_method <- data$detection_time <- data$FUP_detection_time <- NA 
  
  # find earliest time between 3 types tests and between symptomatic detection
  # create vectors CT, CEA, Colonoscopy detection times
  ct <- rep(Inf, n_patients)
  cea <- rep(Inf, n_patients)
  colo <- rep(Inf, n_patients)
  
  # For non-cured patients, we need to compute test times
  # Find the first positive index & Assign the corresponding time
  for (i in seq_len(n_patients)) {
    if (data$cured_status[i] == 1) {  # Non-cured patients (cured=0, recurrence=1)
      # CT test times
      ct_res_i <- ct_res$results[[i]]
      pos_ct <- match(1, ct_res_i)
      if (!is.na(pos_ct)) ct[i] <- ct_res$tests[[i]][pos_ct]
      
      # CEA test times
      cea_res_i <- cea_res$results[[i]]
      pos_cea <- match(1, cea_res_i)
      if (!is.na(pos_cea)) cea[i] <- cea_res$tests[[i]][pos_cea]
      
      # Colonoscopy test times
      colo_res_i <- colo_res$results[[i]]
      pos_colo <- match(1, colo_res_i)
      if (!is.na(pos_colo)) colo[i] <- colo_res$tests[[i]][pos_colo]
    }
  }
  
  
  data$FUP_detection_time <- pmin(ct, cea, colo, na.rm = TRUE) # find first follow-up detection time for each patient (min of CT, CEA, Colonoscopy)
  
  all_times_matrix <- cbind(ct, cea, colo,  symp=data$symprec_times) # Create a matrix with test times and symptom times
  data$detection_time <- apply(all_times_matrix, 1, min, na.rm = TRUE)# Calculate overall detection time (min of all times for each patient - including symptom)
  
  # Identify the detection method based on the min detection time
  # Apply condition for detection method
  data$detection_method <- apply(all_times_matrix, 1, function(times) {
    min_time <- min(times, na.rm = TRUE)
    if (!is.finite(min_time)) return('')  # No detection
    # Check if symptomatic is the earliest (if Symp matches min time)
    if (times[4] == min_time) return("Symp")
    # Otherwise, return the tests that have tied for the earliest time
    tied_tests <- names(times)[times == min_time & !is.na(times)]
    paste(tied_tests, collapse = "+")
  })
  
  # for patients who never had a recurrence detected (time=infinity), fill entry with blank space (easier to read)
  data$detection_method[data$detection_time=='Inf'] <- ""
  
  ###### For each patient, only keep test results before the first positive result -------------------------------------------------
  # (first positive result is same as detection_time)
  # i.e., If a CEA detects recurrence at 6 months, then no other tests are also performed (e.g., CT and colo, same for CT detecting, then no more CEA, 6 months is just example)
  for (i in 1:nrow(data)) {
    if (data$detection_time[i] != Inf) {
      stop_month <- floor(as.numeric(data$detection_time[i]))
      keep_ct <- ct_res$tests[[i]] <= stop_month
      ct_res$tests[[i]]   <- ct_res$tests[[i]][keep_ct]
      ct_res$results[[i]] <- ct_res$results[[i]][keep_ct]
      
      keep_cea <- cea_res$tests[[i]] <= stop_month
      cea_res$tests[[i]]   <- cea_res$tests[[i]][keep_cea]
      cea_res$results[[i]] <- cea_res$results[[i]][keep_cea]
      
      keep_colo <- colo_res$tests[[i]] <= stop_month
      colo_res$tests[[i]]   <- colo_res$tests[[i]][keep_colo]
      colo_res$results[[i]] <- colo_res$results[[i]][keep_colo]
      
    }
  }
  
  ###### Curative treatment: probability of receiving curative treatment -------------------------------------------------
  #                        group proportion lower_CI upper_CI
  # 1        FUP- LRR & CUR      0.600    0.434    0.747
  # 2     FUP- LRR & NONCUR      0.400    0.253    0.566
  # 3         FUP- DR & CUR      0.557    0.484    0.628
  # 4      FUP- DR & NONCUR      0.443    0.372    0.516
  # 5       FUP- 2REC & CUR      0.180    0.109    0.279
  # 6   FUP- 2REC & NONCURe      0.820    0.721    0.891
  # 7       SYMP- LRR & CUR      0.333    0.155    0.569
  # 8    SYMP- LRR & NONCUR      0.667    0.431    0.845
  # 9        SYMP- DR & CUR      0.172    0.090    0.299
  # 10    SYMP- DR & NONCUR      0.828    0.701    0.910
  # 11     SYMP- 2REC & CUR      0.067    0.012    0.235
  # 12 SYMP- 2REC & NONCURe      0.933    0.765    0.988
  
  data$detection_type <- ifelse(data$detection_method == 'Symp', "Symp", ifelse(data$detection_method != '', "FUP", ""))

  # store the probability of being cured for each combination of locnum and detection type 
  cured_prob <- case_when(
    data$detection_type == 'Symp' & data$locnum == 0 ~ 0.33,
    data$detection_type == 'Symp' & data$locnum == 1 ~ 0.172,
    data$detection_type == 'Symp' & data$locnum == 2 ~ 0.067, # 2/30 = 0.067
    data$detection_type == 'FUP' & data$locnum == 0 ~ 0.6,
    data$detection_type == 'FUP' & data$locnum == 1 ~ 0.557,
    data$detection_type == 'FUP' & data$locnum == 2 ~ 0.180,
    data$detection_type == '' ~ 0
  )
  # draw from binomial distributions based on combinations of detection type and locnum using the specific cured_prob for each patient 
  curtreat <- rbinom(n_patients, 1, cured_prob) 
  ### add later ### ---> can we also make it flexible so that we can incorporate the 95%-CI in a PSA?

  data$curtreat <- ifelse(curtreat == 1, "curative",ifelse(data$detection_type!='', "non-curative", ''))
  
  # DOD_simulate <- function(shape, scale, eta, n_sim) { # beta changes for curative or non-curative
  #   U <- runif(n_sim)
  #   scale_adj <- scale * exp(eta)  # adjusted scale for number of recurrences
  #   t <- scale_adj * (-log(1 - U))^(1 / shape)
  #   return(t)
  # }
  
  combos <- expand.grid( # make all combinations of locnum, detection_type, and curtreat 
    locnum = c(0, 1, 2),
    detection_type = c('FUP', 'Symp'),
    curtreat = c('curative', 'non-curative'),
    stringsAsFactors = FALSE
  )
  
  combos$cur_beta <- case_when( # Add corresponding beta values for each combination
    combos$curtreat == 'curative' & combos$detection_type == 'FUP' ~ beta_cure_fup,
    combos$curtreat == 'non-curative' & combos$detection_type == 'FUP' ~ beta_noncure_fup,
    combos$curtreat == 'curative' & combos$detection_type == 'Symp' ~ beta_cure_symp,
    combos$curtreat == 'non-curative' & combos$detection_type == 'Symp' ~ beta_noncure_symp
  )
  
  combos$locnum_beta <- case_when(
    combos$locnum == 0 ~  0,
    combos$locnum == 1 ~ (beta_1meta),
    combos$locnum == 2 ~ (beta_2meta)
  )
  
  combos$overall_beta <- combos$cur_beta + combos$locnum_beta
  
  DOD_time <- rep(Inf, n_patients)# create empty vector for DOD times only for those that had a detection 
  
  # Loop over each combination of locnum, detection_type, and curtreat (12 combos)
  for (i in 1:nrow(combos)) {
    # Identify rows that match the current combination
    matching_rows <- data$locnum == combos$locnum[i] &
      data$detection_type == combos$detection_type[i] &
      data$curtreat == combos$curtreat[i]
  
    DOD_time[matching_rows]  <- rweibull(n = sum(matching_rows), shape = DOD_shape, scale= DOD_scale*exp(combos$overall_beta[i]))
    # simulate DOD times for the matching rows --> that combination of locnum, detection_type, and curtreat (12 combos)
  }
  
  data$DOD_time_raw <- DOD_time 
  data$DOD_time_raw[data$DOD_time_raw >  60] <- Inf # if death 5 years after detection then set at infinity  
  
  # add DOD time to time of detection (so that DOD time is not before recurrence time)
  # because DOD simulates from detection time not resection time
  data$DOD_time_corrected <- data$DOD_time_raw + data$detection_time 
  
  ###### final time of death (death_time): -------------------------------------------------
  # check which was earlier: DOC (death other causes), mor90d_time (90 day mortality), DOD (death of disease)
  # save this time as final death time, and save reason for death 
  death_times <- cbind(data$mor90d_time, data$DOC_time, data$DOD_time_corrected)
  data$death_time <- do.call(pmin, as.data.frame(death_times)) # find minimum of the three times 
  data$death_cause <- c("90day", "DOC", "DOD")[max.col(-death_times, ties.method='first')] # get corresponding cause of the minimum
  data$death_cause[data$mor90d_time ==Inf & data$DOC_time==Inf &  data$DOD_time_corrected==Inf] <- 'DOC'
  
  data$CRC_death <- ifelse(data$death_cause == 'DOD', 1, 0)
  
  # count number of each type of test performed before first detection
  data$num_ct <- sapply(ct_res$tests, function(x) length(x))
  data$num_cea <-  sapply(cea_res$tests, function(x) length(x))
  data$num_colo <-  sapply(colo_res$tests, function(x) length(x))
  
  data$time60month <- pmin(data$death_time, 60)
  data$status60month <- ifelse(data$death_time < 60, 1, 0)
  
  data$rec_to_death <- data$DOD_time_raw
  
  # count number of false positives per test type
  data$num_ct_FP <- sapply(ct_res$results, function(x) sum(x == 8))
  data$num_cea_FP <- sapply(cea_res$results, function(x)  sum(x == 8))
  data$num_colo_FP <- sapply(colo_res$results, function(x)  sum(x == 8))
  
  return(list(data,  cea_res))
}
