#####
# PROJECT. Precision Medicine for Knee OA
# AUTHOR. Xiaotong Jiang
#####
# USAGE. Generate simulated data. Run main.R for PMOA.
#####
# Simulation Data: 
# 4 data types {circle, steps, line, quadratic} 
# 3 X's, X3 is nuisance
# 3 treatments (independent of X), 0 is considered standard of care and would be the control group
#####

F1_genSimData <- function(n.sim, n, data.type){
  
  #' Generate simulated dataset according to supplemental material
  #' 
  #' @param n.sim number of simulations (set to be 1 for example code; was 100 in the manuscript)
  #' @param n sample size, i.e. number of rows
  #' @param data.type type of data, {"circle", "steps", "line", "quadratic", "null"}; "null" means one treatment is better than the other for all X values, not reported in the manuscript
  #' 
  #' @return simdat.list -- a list of training data and independent test data of same size; if n.sim != 1, training data and test data are lists of n.sim simulation copies
  
  library(tidyverse)
  
  ##########################
  #### Define constants ####
  ##########################
  
  size.rule <- c(2L, 3L, 5L, 10L)   # number of list nodes for LIST models
  error.sd <- 1                         # standard deviation for error term
  num_x = 3                             # number of covariates; at least 2; x1 and x2 determine the decision boundary, rest are nuisance variables
  num_trt = 3                           # number of treatments
  propensity <- rep(1/num_trt, num_trt) # propensity score of the treatments, assume equal probability
  
  ##################################
  #### Generate simulation data ####
  ##################################
  
  for (i in 1:n.sim){
    set.seed(seed.fn(1000, i))
    
    if ((i == 1) & (n.sim != 1)){
      dat.list = vector("list", length = n.sim)
      test.dat.list = vector("list", length = n.sim)
    }
    
    dat <- data.generator(n_ = n, propensity_ = propensity, num_trt_ = num_trt, num_x_ = num_x, 
                          type_ = data.type, error.sd_ = error.sd) 
    test.dat <- data.generator(n_ = n, propensity_ = propensity, num_trt_ = num_trt, num_x_ = num_x, 
                               type_ = data.type, error.sd_ = error.sd)
    
    if (n.sim != 1){
      dat.list[[i]] <- dat
      test.dat.list[[i]] <- test.dat
    }
  }
  
  if (n.sim == 1){output = list(dat, test.dat)} else {output = list(dat.list, test.dat.list)}

  return(output)
}

