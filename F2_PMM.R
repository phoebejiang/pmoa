#####
# PROJECT. Precision Medicine for Knee OA
# AUTHOR. Xiaotong Jiang
#####
# USAGE. Fit precision medicine model (PMM) using jackknife value function estimator. Run main.R for PMOA.
#####



F2_PMM <- function(n.sim = 1, n, data.type, models, dat, test.dat, estimator, num_rep, num_folds, n.try, n.pick){
  
  #' Fit precision medicine model with the jackknife method
  #' 
  #' @param n.sim integer, number of simulations (set to be 1 for example code; was 100 in the manuscript)
  #' @param n integer, sample size, i.e. number of rows
  #' @param data.type string, type of data, {"circle", "steps", "line", "quadratic", "null"}; "null" means one treatment is better than the other for all X values 
  #' @param models a vector of PMM models of user's choice
  #' @param dat data.frame, training data, where covariate names start with "x", treatment name is "a" with values 0:2, and outcome name is "y"; used for all four estimators
  #' @param test.dat data.frame, test data, where covariate names start with "x", treatment name is "a" with values 0:2, and outcome name is "y"; used for jackknife+test and empirical+test estimators 
  #' @param estimator string, which estimator to use for the value function, originally {"empirical", "jackknife", "empirical+test", "jackknife+test"}; for the purpose of this example code, estimator is set to be "jackknife"
  #' @param num_rep integer, number of repetitions for cross validation in the LIST + super learning model
  #' @param num_folds integer, number of cross validation folds in the LIST+superlearning model
  #' @param n.try integer, randomly generate n.try alpha (base learner weight)
  #' @param n.pick integer, pick the best n.pick alpha (base learner weight)
  #' 
  #' @return ds, data.frame, long version of value.object, a list of PMM results, (components: sim.param, data.param, U, W, A)
  #' @return mean_var, data.frame of n.models x 7 (cols = model, estimator, variance, mean, V, n, n.sim)
  
  
  ##########################
  #### Define constants ####
  ##########################
  
  size.rule <- c(2L, 3L, 5L, 10L)       # number of list nodes 
  error.sd <- 1                         # standard deviation for error term
  num_trt = 3                           # number of treatments
  propensity <- rep(1/num_trt, num_trt) # propensity score of the treatments, assume equal probability
  
  n.models <- length(models)
  
  ## Create empty shells
  tmp.list = array(NA, c(n.sim, n, length(models), 3),
                   dimnames = list(paste0("sim", 1:n.sim),
                                   paste0("test.",c(paste0("jk", 1:n))),
                                   models,
                                   c("Ui", "Wi", "d.hat")))
    # n.sim x n x n.models x 3 (U or W or d.hat)
  
  value.object <- list(sim.param = list(n.sim = n.sim, sample.size = n, data.type = data.type),
                       data.param = list(propensity = propensity, estimator = estimator,
                                         num_trt = num_trt, size.rule = size.rule, 
                                         seed.fn = seed.fn, error.sd = error.sd),
                       U = NA,
                       W = NA,
                       A = NA)
  
  tt(1)
  
  ds <- NULL
  
  ##################################
  #### Apply PMM with Jackknife ####
  ##################################
  
  #### Fit PMM and get estimated value functions (Ui, Wi, Zi) ####
  for (i in 1:n.sim){ 
    if (class(dat) == "list"){
      data = dat[[i]]
      test.data = test.dat[[i]]
    } else if (class(dat) == "data.frame"){
      data = dat
      test.data = test.dat
    }
    cat("Simulation: No.", i, "Estimator:", estimator, "\n")
    
    for (j in 1:length(models)){
      model.fullname <- models[j]
      cat("Model:", model.fullname, "\n")
      
      model <- ifelse(grepl("^list[0-9]+\\+.*", model.fullname), 
                      sub("[0-9]+\\+.*", "", model.fullname), # if list+rf
                      sub("[\\.|-].*$", "", model.fullname)) # if not list+rf
      
      node <- ifelse(grepl("^list[0-9]+\\+.*$", model.fullname), 
                     model.fullname %>% str_match_all("[0-9]+") %>% unlist(),  # if list+rf
                     str_split(model.fullname, "\\.")[[1]][2]) # if not list+rf
      
      Q.model <- ifelse(grepl("^list[0-9]+\\+", model.fullname), 
                        sub("^list[0-9]+\\+", "", model.fullname), 
                        NA)
      
      sl.bundle <- NULL # only used for super learning model
      if ((Q.model == "sl") & !is.na(Q.model)){sl.bundle = list(data, num_rep, num_folds, n.try, n.pick, NA)} # remove_ids initialized to be NA for right now
      
      which_wl <- str_split(model.fullname, "-")[[1]][2]
        
      if (grepl("test", estimator)){test.data = test.data} else {test.data = NULL}
      
      combined <- Get_UW_by_JK(data = data, model = model,
                               size = ifelse(!is.na(node), as.integer(node), NA),
                               empirical = grepl("empirical", estimator),
                               test = grepl("test", estimator),
                               test.data = test.data,
                               kernel = ifelse(!is.na(which_wl), which_wl, NA),
                               Q.model = Q.model, sl.bundle = sl.bundle)
      # print(table(combined$A))
      
      
      tmp.list[i,,j,1] <- combined$U
      tmp.list[i,,j,2] <- combined$W
      tmp.list[i,,j,3] <- combined$A

      #### Summarize results and calculate value function ####
      tmpU <- tmp.list[i,,j,1] %>% reshape2::melt() %>% 
        rownames_to_column(var = "test.n") %>% 
        mutate(sim = i, model = model.fullname, estimator = estimator, which.Z = "Ui") %>% 
        select(sim, test.n, model, value, estimator, which.Z)
      tmpW <- tmp.list[i,,j,2] %>% reshape2::melt() %>% 
        rownames_to_column(var = "test.n") %>% 
        mutate(sim = i, model = model.fullname, estimator = estimator, which.Z = "Wi") %>% 
        select(sim, test.n, model, value, estimator, which.Z)
      each.ds <- rbind(tmpU, tmpW)
      
      ds <- rbind(ds, each.ds)
    }#j loop for models
  }#i loop for partitioned simulations
  
  # Update the empty value object shell
  value.object$U = tmp.list[,,,1] # n.sim x n x n.models
  value.object$W = tmp.list[,,,2] # n.sim x n x n.models
  value.object$A = tmp.list[,,,3] # n.sim x n x n.models
  
  tt(2) %>% print()
  
  #### Variance of the jackknife estimator (across J)
  ds$test.n <- as.integer(sub(pattern = "test.jk", replacement = "", x = ds$test.n))
  stopifnot(n.sim == length(unique(ds$sim)))
  cat("\nn.models", n.models, "n", n, "ds.shape", nrow(ds), ncol(ds), "\n")
  stopifnot(nrow(ds) == n.models * n * 2 * n.sim)
  n.sims.obv <- length(unique(ds$sim))
  
  # First take the mean across J, i.e. n samples and simulations, for each model, estimator, & Z component
  mu <- ds %>% 
    filter(model %in% models) %>%
    group_by(model, estimator, which.Z) %>% 
    summarise(mean = mean(value, na.rm = T)) %>% 
    spread(which.Z, mean) %>% rename(Ubar = Ui, Wbar = Wi) 
  # dim: n.models x 4
  
  # Convert ds from long (U,W separate rows) to wide (U,W same row)
  ds.wide <- ds %>% 
    arrange(sim, test.n, model, which.Z) %>% 
    mutate(id = rep(1:(n()/2), each = 2)) %>% 
    spread(which.Z, value) %>% 
    select(-id) %>% 
    rename(Fbar.U = Ui, Fbar.W = Wi)
  # dim: nrow(ds)/2 x 6
  
  # Calculate F double bar across test.n and simulation (J = mK, where K = n)
  # F bar is Ui, Wi (i.e. Zi) for jackknife
  preps <- ds.wide %>% 
    group_by(model, estimator) %>% 
    summarise(Fbarbar.U = mean(Fbar.U, na.rm = T), Fbarbar.W = mean(Fbar.W, na.rm = T)) %>% 
    right_join(ds.wide, by = c("model", "estimator")) %>% # Match the mean with each value
    mutate(U.diff = Fbar.U - Fbarbar.U, W.diff = Fbar.W - Fbarbar.W, 
           S1 = U.diff^2, S2 = U.diff * W.diff, S4 = W.diff^2) # Calculate the MSE and variance matrix
  # dim: nrow(ds)/2 x 13
  
  # Calculate covariance by entry
  preps2 <- preps %>%  
    group_by(model, estimator) %>% 
    summarise(S1.sum = sum(S1, na.rm = T) / (n*(n*n.sims.obv - 1)), S2.sum = sum(S2, na.rm = T) / (n*(n*n.sims.obv - 1)), S4.sum = sum(S4, na.rm = T) / (n*(n*n.sims.obv - 1))) %>%  # S2 == S3
    left_join(mu, by = c("model", "estimator")) %>% 
    mutate(Q1 = 1/Wbar, Q2 = -Ubar/(Wbar^2)) %>% 
    ungroup() 
  # n = sample size, i.e. 100; n.sims = number of total simulations, i.e. 10*10 = 100.
  # dim: n.models x 9
  
  # Variance is Q^T * Sn * Q
  var <- cbind(preps2[,1:2], 
               variance = apply(preps2[,-(1:2)], 1, function(x) {t(c(x['Q1'], x['Q2'])) %*% matrix(c(x['S1.sum'], x['S2.sum'], x['S2.sum'], x['S4.sum']), nrow = 2) %*% c(x['Q1'], x['Q2'])} ))
  
  #### Mean of the jackknife estimator (across J)
  # Mean is Q^T Zi
  preps3 <- preps %>% 
    left_join(preps2, by = c("model", "estimator")) %>% 
    mutate(Ri = Q1 * Fbar.U + Q2 * Fbar.W) %>% 
    group_by(model, estimator) %>% 
    summarise(mean = mean(Ri, na.rm = T), V = sum(Fbar.U, na.rm = T) / sum(Fbar.W, na.rm = T))
  
  mean_var = left_join(var, preps3 %>% select(model, estimator, mean, V), by = c("model", "estimator")) %>% 
    mutate(n = n, sims = n.sim, sims.obv = n.sims.obv, data.type = data.type) 
  # means should be 0 or very close to 0
  # V is value function summed across all J; 
  # values are summarized by simulations first; 
  # close but a bit different
  
  return(list(ds, mean_var))
}

