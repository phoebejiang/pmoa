#####
# PROJECT. Precision Medicine for Knee OA
# AUTHOR. Xiaotong Jiang
#####
# USAGE. Fit zero-order model (ZOM). Run main.R for PMOA.
#####

F3_ZOM <- function(n.sim = 1, n, data.type, dat, test.dat, estimator){
  
  #' Zero-order model, i.e., all subjects receiving the same treatment
  #' 
  #' @param n.sim number of simulations (set to be 1 for example code; was 100 in the manuscript)
  #' @param n sample size, i.e. number of rows
  #' @param data.type type of data, {"circle", "steps", "line", "quadratic", "null"}; "null" means one treatment is better than the other for all X values 
  #' @param dat training data, where covariate names start with "x", treatment name is "a" with values 0:2, and outcome name is "y"; used for all four estimators
  #' @param test.dat test data, where covariate names start with "x", treatment name is "a" with values 0:2, and outcome name is "y"; used for jackknife+test and empirical+test estimators 
  #' @param estimator which estimator to use for the value function, originally {"empirical", "jackknife", "empirical+test", "jackknife+test"}; for the purpose of this example code, estimator is set to be "jackknife"
  #' 
  #' @return ds, data.frame, long version of tmp.list, an array of ZOM results, dimension n.sim x n x n.models x 2 (U,W)
  #' @return mean_var, data.frame of n.models x 7 (cols = model, estimator, variance, mean, V, n, n.sim)
 
  models <- paste0("zero_order", unique(dat$a))
  n.models = length(models) 
  
  #### Create empty shells ####
  # ZOM U's and W's
  tmp.list = array(NA, c(n.sim, n, length(models), 2),
                   dimnames = list(paste0("sim", 1:n.sim),
                                   paste0("test.",c(paste0("jk", 1:n))),
                                   models,
                                   c("Ui", "Wi")))
  # n.sim x n x n.models x 2 (U or W)
  
  #### Get estimated value functions (Ui, Wi, Zi) ####
  ## Loop through all n.sims in each partition
  for (i in 1:n.sim){ 
    if (class(dat) == "list"){
      data = dat[[i]]
      test.data = test.dat[[i]]
    } else if (class(dat) == "data.frame"){
      data = dat
      test.data = test.dat
    }
    cat("Estimator:", estimator, "\n")
    
    ### Loop through all ZOMs
    for (m in 1:length(models)){
      
      model <- models[m]
      if (grepl("test", estimator)){test.data = test.data} else {test.data = NULL}
      combined <- Get_UW_by_JK_zero_order(data = data,
                                          which_zero = model %>% str_match("[0-9]+") %>% unlist() %>% as.numeric(),
                                          empirical = grepl("empirical", estimator),
                                          test = grepl("test", estimator),
                                          test.data = test.data)
      
      tmp.list[i,,m,1] <- combined$U
      tmp.list[i,,m,2] <- combined$W
      
    }#m loop through ZOMs
  }#i loop though n.sim
  
  #### Summarize results and estimator value function ####
  tmpU <- tmp.list[,,,1] %>% reshape2::melt() %>% mutate(estimator = "jackknife", which.Z = "Ui") 
  tmpW <- tmp.list[,,,2] %>% reshape2::melt() %>% mutate(estimator = "jackknife", which.Z = "Wi")
  if (ncol(tmpU) == 5){
    tmpU %<>% mutate(sim = 1) %>% select(sim, everything())
    tmpW %<>% mutate(sim = 1) %>% select(sim, everything())
  } 
  colnames(tmpU) <- c("sim", "test.n", "model", "value", "estimator", "which.Z")
  colnames(tmpW) <- c("sim", "test.n", "model", "value", "estimator", "which.Z")
  ds <- rbind(tmpU, tmpW)
  ds$test.n <- as.integer(sub(pattern = "test.jk", replacement = "", x = ds$test.n))
  stopifnot(n.sim == length(unique(ds$sim)))
  
  #### Variance of the jackknife estimator (across J) ####
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
    summarise(S1.sum = sum(S1, na.rm = T) / (n*(n*n.sim - 1)), S2.sum = sum(S2, na.rm = T) / (n*(n*n.sim - 1)), S4.sum = sum(S4, na.rm = T) / (n*(n*n.sim - 1))) %>%  # S2 == S3
    left_join(mu, by = c("model", "estimator")) %>%
    mutate(Q1 = 1/Wbar, Q2 = -Ubar/(Wbar^2)) %>%
    ungroup()
  # n = sample size, i.e. 100; n.sim = number of total simulations, i.e. 10*10 = 100.
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
    mutate(n = n, n.sim = n.sim, data.type = data.type)
  # means should be 0 or very close to 0
  # V is value function summed across all J;
  # values in data3 are summarized by simulations first;
  # close but a bit different
  
  return(list(ds, mean_var))
}

  