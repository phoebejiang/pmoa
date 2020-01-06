#####
# PROJECT. Precision Medicine for Knee OA
# AUTHOR. Xiaotong Jiang
#####
# USAGE. To provide utility functions for all F functions (not including the super learning part). Run main.R for PMOA.
#####

library(tidyverse)
library(magrittr)
library(randomForest)
library(listdtr)
library(varhandle) # unfactor()
library(data.table) # %like%
library(DynTxRegime)
library(glmnet)
library(BART)
library(RLT)

### 1. Define sample and population data
## sample data (for training)
data.generator <- function(n_ = 1000, propensity_ = 0.5, num_trt_ = 3, num_x_ = 3, type_ = "circle", error.sd_) {
  # set.seed(seed_)
  x.nus <- NULL
  stopifnot(num_x_ >= 2)
  
  # Create feature variables
  x1 = runif(n_)*4 - 2                 # Xj ~ U(-2, 2) i.i.d.
  x2 = runif(n_)*4 - 2
  for (j in 1:(num_x_ - 2)){
    x.nu = runif(n_)*4 - 2                # nuisance variable
    x.nus = cbind(x.nus, x.nu)
  }#j for nuisance variables
  
  sim.data <- data.frame(x1, x2, x.nus)
  colnames(sim.data) <- paste0("x", 1:num_x_)
  
  if (num_trt_ == 2){
    # Treatment and error
    a = rbinom(n_, 1, propensity_)      # A ~ Binom(n, 0.5)
    error = rnorm(n = n_, mean = 0, sd = error.sd_)
    
    # Outcomes
    if (type_ == "circle") {
      EY1 = (x1 + x2) + 3 - x1^2 - x2^2   # decision boundary: Treat if (X1^2 + X2^2 - 2 < 0)
      EY0 = (x1 + x2)                     # value = pi*(int_0^3 (3-y) dy = 9/2) /16 = 9pi/32 = 0.8836
    } else if (type_ == "steps") {
      EY1 = (x1 + x2) + 1*(x2 <= ceiling(x1)) - 1*(x2  > ceiling(x1)) 
      # decision boundary: Treat if (X2 - ceiling(X1) <0)
      EY0 = (x1 + x2)                     # value = 10/16 = 5/8 = 0.6250
    } else if (type_ == "line") {
      EY1 = (x1 + x2) + (1 - x1 - x2)
      EY0 = (x1 + x2)
    } else if (type_ == "quadratic") {
      EY1 = (x1 + x2) - x1^2 + 1 + x2
      EY0 = (x1 + x2)
    } else if (type_ == "null") {
      EY1 = (x1 + x2) + 0.5
      EY2 = (x1 + x2) 
    } else stop("No valid type is chosen.")
    
    y = ifelse(a, EY1, EY0) + error
    
    sim.data <- cbind(sim.data, a, y, EY1, EY0, error)
  }# if (num_trt == 2)
  
  if (num_trt_ == 3){
    # Treatment and error
    if (length(propensity_) != 3){cat("Wrong dimension of propensity score!!")}
    
    a = c(0,1,2) %*% rmultinom(n_, 1, propensity_) %>% as.vector() # A ~ Multinom(n_, c(1/3, 1/3, 1/3))
    error = rnorm(n = n_, mean = 0, sd = error.sd_)
    
    # Outcomes
    if (type_ == "circle") {
      EY0 = (x1 + x2) 
      EY1 = EY0 + (3 - x1^2 - x2^2)*(x1^2 + x2^2 - 1)   
      EY2 = EY0 + (1 - x1^2 - x2^2)
    } else if (type_ == "steps") {
      EY0 = (x1 + x2)
      EY1 = EY0 + 1 * (x2 <= ceiling(x1)) - 1 * (x2  > ceiling(x1)) 
      EY2 = EY0 + 1 * (x2 <= ceiling(x1 - 2)) - 1 * (x2 > ceiling(x1 - 2))
    } else if (type_ == "line") {
      EY0 = (x1 + x2)
      EY1 = EY0 + (1 - x1 - x2) * (x1 + x2 + 1)
      EY2 = EY0 + (x1 + x2 - 1)
    } else if (type_ == "quadratic") {
      EY0 = (x1 + x2)
      EY1 = EY0 + (x2 - x1^2 + 2) * (x1^2 - x2)
      EY2 = EY0 + (x2 - x1^2)
    } else if (type_ == "null") {
      EY0 = (x1 + x2)
      EY1 = (x1 + x2) + 0.5
      EY2 = (x1 + x2) 
    } else stop("No valid type is chosen.")
    
    y = ifelse(a == 0, EY0, ifelse(a == 1, EY1, EY2)) + error
    
    sim.data <- cbind(sim.data, a, y, EY2, EY1, EY0, error)
  }# if (num_trt == 3)
  
  return(sim.data)
}


## population data

# Grid Method
# expand.grid(seq(-2,2, by=0.01),seq(-2,2, by=0.01)) %>% 
#   transmute(x1 = Var1, x2 = Var2, a=NA) -> pop

# Random Method
population.generator <- function(n_, propensity_, num_trt_, num_x_, type_, seed_){
  set.seed(seed_)
  x.nus <- NULL
  stopifnot(num_x_ >= 2)
  
  # Create feature variables
  x1 = runif(n_)*4 - 2                 # Xj ~ U(-2, 2) i.i.d.
  x2 = runif(n_)*4 - 2
  
  for (j in 1:(num_x_ - 2)){
    x.nu = runif(n_)*4 - 2                # nuisance variable
    x.nus = cbind(x.nus, x.nu)
  }#j for nuisance variables
  
  sim.data <- data.frame(x1, x2, x.nus)
  colnames(sim.data) <- paste0("x", 1:num_x_)
  
  if (num_trt_ == 2){
    pop <- cbind(sim.data, A = rbinom(n = n_, size = 1, p = propensity_))
    
    if (type_ == "circle"){
      ## circle data
      pop.1 <- pop %>% mutate(a=1, Y1 = (x1 + x2) + 3 - x1^2 - x2^2) 
      pop.0 <- pop %>% mutate(a=0, Y0 = (x1 + x2)) 
      #0.8791 /0.8836, random 0.8862
    } else if (type_ == "steps"){
      ## steps data
      pop.1 <- pop %>% mutate(a=1, Y1 = (x1 + x2) + 1*(x2 <= ceiling(x1)) - 1*(x2  > ceiling(x1)))
      pop.0 <- pop %>% mutate(a=0, Y0 = (x1 + x2))
      #0.6244 /0.6250, random 0.6292
    } else if (type_ == "line") {
      ## line data
      pop.1 = pop %>% mutate(a=1, Y1 = (x1 + x2) + (1 - x1 - x2))
      pop.0 = pop %>% mutate(a=0, Y0 = (x1 + x2))
    } else if (type_ == "quadratic") {
      ## quadratic data
      pop.1 = pop %>% mutate(a=1, Y1 = (x1 + x2) - x1^2 + 1 + x2)
      pop.0 = pop %>% mutate(a=0, Y0 = (x1 + x2))
    } else if (type_ == "null") {
      pop.1 = (x1 + x2) + 0.5
      pop.0 = (x1 + x2) 
    }else {print("ERROR: wrong data.type!")}
    return(list(pop.0 = pop.0, pop.1 = pop.1))
  }
  
  if (num_trt_ == 3){
    if (length(propensity_) != 3){cat("Wrong dimension of propensity score!!")}
    
    pop <- cbind(sim.data, 
                 A = c(0,1,2) %*% rmultinom(n_, 1, propensity_) %>% as.vector()) 
    # A ~ Multinom(n_, c(1/3, 1/3, 1/3)))
    
    if (type_ == "circle"){
      ## circle data
      pop.2 <- pop %>% mutate(a=2, Y2 = (x1 + x2) + 1 - x1^2 - x2^2)
      pop.1 <- pop %>% mutate(a=1, Y1 = (x1 + x2) + (3 - x1^2 - x2^2)*(x1^2 + x2^2 - 1))
      pop.0 <- pop %>% mutate(a=0, Y0 = (x1 + x2)) 
    } else if (type_ == "steps"){
      ## steps data
      pop.1 <- pop %>% mutate(a=1, Y1 = (x1 + x2) + 1*(x2 <= ceiling(x1)) - 1*(x2  > ceiling(x1)))
      pop.2 <- pop %>% mutate(a=2, Y2 = (x1 + x2) + 1*(x2 <= ceiling(x1 - 2)) - 1*(x2 > ceiling(x1 - 2)))
      pop.0 <- pop %>% mutate(a=0, Y0 = (x1 + x2))
    } else if (type_ == "line") {
      ## line data
      pop.1 = pop %>% mutate(a=1, Y1 = (x1 + x2) + (1 - x1 - x2) * (x1 + x2 + 1))
      pop.2 = pop %>% mutate(a=2, Y2 = (x1 + x2) + (x1 + x2 - 1))
      pop.0 = pop %>% mutate(a=0, Y0 = (x1 + x2))
    } else if (type_ == "quadratic") {
      ## quadratic data
      pop.2 = pop %>% mutate(a=2, Y2 = (x1 + x2) + (x2 - x1^2))
      pop.1 = pop %>% mutate(a=1, Y1 = (x1 + x2) + (x2 - x1^2 + 2)*(x1^2 - x2) )
      pop.0 = pop %>% mutate(a=0, Y0 = (x1 + x2))
    } else if (type_ == "null") {
      pop.2 = pop %>% mutate(a=2, Y2 = (x1 + x2))
      pop.1 = pop %>% mutate(a=1, Y1 = (x1 + x2) + 0.5)
      pop.0 = pop %>% mutate(a=0, Y0 = (x1 + x2))
    } else {print("ERROR: wrong data.type!")}
    return(list(pop.0 = pop.0, pop.1 = pop.1, pop.2 = pop.2))
  }# if (num_trt == 3)
}


### 2. Define main functions

## Estimate propensity function (in the denominator of value)
Get_pi.hat <- function(A, prop_data, num_trt){
  stopifnot(unique(A) %in% 0:num_trt)
  list <- prop_data[prop_data$a == A, "count"]
  # list <- ifelse(A == 0, prop_data[prop_data$a == 0, "count"], prop_data[prop_data$a == 1, "count"])
  return(unlist(list))
}

## Get the estimated U and W for jackknife algorithm (only 1 person per test set)
Get_Each_UW <- function(R, A, d.hat, which_one, prop_data, num_trt){
  W = ifelse(d.hat == A, 1, 0) / Get_pi.hat(A, prop_data, num_trt)
  U = R * W
  if (which_one == "U"){return(U)} else 
    if (which_one == "W"){return(W)} else {print("ERROR: wrong value of which_one!")}
}

## Define function for RWL, used in Get_UW_by_JK(), only for covariates that start with "x"
runRWL <- function(data, x.train, outcome_name, group_name, vars, kernel){
  moPropen <- buildModelObj(model = ~1,
                            solver.method = 'glm',
                            solver.args = list('family'='binomial'),
                            predict.method = 'predict.glm',
                            predict.args = list(type='response'))
  
  RWL.formula <- as.formula(paste0("~", paste(colnames(x.train), collapse="+")))
  
  # Create modelObj object for main effect component
  moMain <- buildModelObj(model = RWL.formula, solver.method = 'lm')
  
  if (kernel == "poly2"){
    rwlRes <- rwl(moPropen = moPropen, moMain = moMain,
                  data = data, reward = outcome_name, txName = group_name,
                  kernel = "poly", cvFolds = 0L, kparam = 2, lambdas = 2,
                  regime = RWL.formula, verbose = F)
  } 
  if (kernel == "linear"){
    rwlRes <- rwl(moPropen = moPropen, moMain = moMain,
                  data = data, reward = outcome_name, txName = group_name,
                  regime = RWL.formula, kernel = "linear", verbose = F)
  }
  
  return(rwlRes)
}


## Modified utility function for the LISTDTR model in listdtr package
Modified_listdtr <- function (y, a, x, stage.x, Q.model, maxlen = node.size, seed = seed_, kfolds = 5L, fold = NULL, zeta.choices = NULL, eta.choices = NULL, sl.bundle = sl.bundle) {
  if (!is.matrix(y)) { 
    y <- as.matrix(y)
  }
  if (!is.data.frame(a)) {
    a <- as.data.frame(a)
  }
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  if (is.null(colnames(x))) {
    colnames(x) <- paste0("x", seq_len(ncol(x)))
  }
  stopifnot(nrow(y) == nrow(a) && nrow(y) == nrow(x))
  n <- nrow(y)
  stopifnot(ncol(y) == ncol(a))
  n.stage <- ncol(y)
  dtr <- vector("list", n.stage)
  future.y <- double(n)
  if (is.null(colnames(a))) {
    colnames(a) <- paste0("a", 1L:n.stage)
  }
  a.mm <- lapply(1L:n.stage, function(j) model.matrix(as.formula(sprintf("~ -1 + %s", colnames(a)[j])), a))
  stage.a.mm <- rep.int(1L:n.stage, sapply(a.mm, ncol))
  a.mm <- do.call("cbind", a.mm)
  if (is.null(fold)) {
    if (!is.null(seed)) {
      set.seed(seed)
    }
    fold <- rep_len(1L:kfolds, n)[sample.int(n)]
  }
  for (i.stage in n.stage:1L) {
    current.x <- cbind(x[, which(stage.x <= i.stage), drop = FALSE], 
                       a.mm[, which(stage.a.mm < i.stage), drop = FALSE], 
                       y[, seq_len(i.stage - 1L), drop = FALSE])
    if (ncol(current.x) < 2L) {
      current.x <- cbind(x = current.x, dummy_ = 0)
    }
    current.a <- a[, i.stage]
    current.y <- y[, i.stage] + future.y
    outcomes <- NULL
    
    # Modified to accept different Q.models (more than just "krr")
    if (Q.model == "krr"){
      
      krr.mod <- krr(current.x, current.y, current.a)
      options <- krr.mod$options
      outcomes <- predict(krr.mod, current.x)
      
    } else if (Q.model == "rf"){
      
      ds <- cbind(current.y, current.a, current.x) %>% as.data.frame()
      colnames(ds)[1:2] <- c("y", "a") 
      
      rf.mod <- randomForest(
        formula = as.formula(paste0("y ~ ", paste0(colnames(current.x), collapse = "+"), "+ a")), data = ds, proximity = T) 
      
      pred0 <- predict(rf.mod, ds %>% mutate(a = 0))
      pred1 <- predict(rf.mod, ds %>% mutate(a = 1))
      pred2 <- predict(rf.mod, ds %>% mutate(a = 2))
      
      outcomes <- cbind(pred0, pred1, pred2)
      a.vals <- unique(as.matrix(current.a))
      colnames(outcomes) <- sort(a.vals)
      
    } else if (Q.model == "sl"){
      
      dat <- sl.bundle[[1]] # the original dataset
      num_rep <- sl.bundle[[2]] # number of repetition
      num_folds <- sl.bundle[[3]] # number of cross validation
      num_trt <- 3 
      n.try <- sl.bundle[[4]] %>% as.numeric()
      n.pick <- sl.bundle[[5]] %>% as.numeric()
      remove_ids <- sl.bundle[[6]]
      
      outcomes <- Super_Learning("min_mse", dat, num_rep, num_folds, num_trt, remove_ids, 
                                 seed = seed, n.try = n.try, n.pick = n.pick, plot = FALSE) %>% as.matrix()
      
    } else if (Q.model == "en"){
      
      #  Convert a to factor with 0 as reference group
      data.all <- cbind(current.x, a = current.a) %>% as.data.frame()
      data.all$a = as.factor(data.all$a) %>% relevel(ref = "0")
      
      # Create interactions
      f <- as.formula(~ x1*a + x2*a + x3*a)
      x.train <- model.matrix(f, data.all) 
      y.train <- as.matrix(current.y, ncol=1)
      
      # Training
      en.mod <- glmnet(x.train, y.train, alpha = 0.5, family = "gaussian")
      cv = cv.glmnet(x.train, y.train, type.measure = "mse", nfolds = 5)
      
      # Testing
      r <- nrow(data.all)
      
      x.test <- suppressWarnings(data.all %>% slice(rep(row_number(), 3)) %>% 
        mutate(a = as.factor(rep(0:2, each = r))) %>% 
        model.matrix(f, .))
      
      preds <- predict(en.mod, x.test, s = cv$lambda.min)
      outcomes <- cbind(preds[1:r,], preds[(r+1):(2*r),], preds[(2*r + 1): (3*r),])
  
      a.vals <- unique(as.matrix(current.a))
      colnames(outcomes) <- sort(a.vals)

    } else {print("ERROR: Invalid Q.model!")}
    
    regrets <- get.regrets(outcomes)
    obj <- build.rule.cv(current.x, regrets, kfolds, fold, maxlen, zeta.choices, eta.choices)
    dtr[[i.stage]] <- obj
    future.y <- outcomes[cbind(1L:n, obj$action)]
    
  }# end of stage loop
  class(dtr) <- "listdtr"
  
  return(dtr)
}

## Utility function from the listdtr package
get.regrets <- function(outcomes){
  regrets <- .Call("R_get_regrets_from_outcomes", outcomes)
  colnames(regrets) <- colnames(outcomes)
  
  return(regrets)
}

## Modified function for the function in listdtr package
# Added cv.loss = NULL for when simple.loss is small
build.rule.cv <- function (x, y, kfolds = 5L, fold = NULL, maxlen = 10L, zeta.choices = NULL, 
                           eta.choices = NULL, cv.only = FALSE) 
{
  if (!is.matrix(x) || ncol(x) < 2) {
    x <- cbind(x = x, dummy_ = 0)
  }
  y <- as.matrix(y)
  stopifnot(ncol(y) > 1)
  stopifnot(nrow(x) == nrow(y))
  n <- nrow(x)
  simple.loss <- min(colMeans(y))
  if (simple.loss < 1e-08) {
    zeta.selected <- simple.loss * n
    eta.selected <- simple.loss * n
    cv.loss <- NULL # added by XJ
  } else {
    if (is.null(zeta.choices) || is.null(eta.choices)) {
      zeta.grid <- simple.loss * c(2, 0.75, 0.3, 0.12, 
                                   0.05)
      eta.grid <- simple.loss * n * c(0.3, 0.1, 0.03)
      zeta.choices <- rep(zeta.grid, times = 3L)
      eta.choices <- rep(eta.grid, each = 5L)
    }
    if (is.null(fold)) {
      kfolds <- as.integer(kfolds)
      fold <- rep_len(1L:kfolds, n)[sample.int(n)]
    } else {
      fold <- as.integer(fold)
      kfolds <- max(fold)
      if (any(tabulate(fold, kfolds) <= 5L)) {
        stop("Some fold(s) have too few observations.")
      }
    }
    fold <- fold - 1L
    cv.loss <- .Call("R_cv_tune_rule", x, y, zeta.choices, 
                     eta.choices, maxlen, fold, kfolds)
    min.cv.loss <- min(cv.loss)
    if (min.cv.loss > simple.loss - 1e-08) {
      zeta.selected <- simple.loss * n
      eta.selected <- simple.loss * n
    }
    else {
      index <- which(cv.loss - min.cv.loss - 1e-08 <= 0)[1L]
      zeta.selected <- zeta.choices[index]
      eta.selected <- eta.choices[index]
    }
  }
  cv <- list(zeta.selected = zeta.selected, eta.selected = eta.selected, 
             metrics = data.frame(zeta.choices = zeta.choices, eta.choices = eta.choices, 
                                  cv.loss = cv.loss))
  if (cv.only) {
    object <- list(cv = cv)
  } else {
    object <- build.rule(x, y, maxlen, zeta.selected, eta.selected)
    object$cv <- cv
  }
  object
}



## Training (for all models)
Training <- function(data.train, model, size, kernel, Q.model, sl.bundle){
  num_trt <- length(unique(data.train$a))
  
  if (model == "rf"){
    x.train <- data.train %>% select(starts_with("x"))
    
    original.mtry <- floor((ncol(x.train) + 1)/3)
    if (original.mtry <= 1){new.mtry = 2}
    
    rf.formula <- as.formula(
      paste0("y~ a+", paste0(x.train %>% colnames(), collapse = "+")))
    
    mod.train <- randomForest(rf.formula, data = data.train %>% as.matrix(), mtry = new.mtry)
  } else if (model == "rlt") {
    mod.train <- RLT::RLT(x = data.train %>% select(starts_with("x"), a), y = data.train[["y"]], reinforcement = T, use.cores = 1, combsplit = 2, ntrees = 50)
  } else if ((model == "list") & is.na(Q.model)){
    mod.train = tryCatch(
      {listdtr(y = data.train$y, a = data.train$a, 
               x = data.train %>% dplyr::select(starts_with("x")) %>% as.matrix, 
               stage.x = rep(1, ncol(data.train %>% dplyr::select(starts_with("x")))), maxlen = size)
      }, error = function(error_message){
        message("\nHere's the original error message:")
        message(error_message)
        return(NA)
      } 
    )
  } else if ((model == "list") & (!is.na(Q.model))){
    mod.train = tryCatch(
      {Modified_listdtr(y = data.train$y, a = data.train$a, 
                        x = data.train %>% select(starts_with("x")) %>% as.matrix(), 
                        stage.x = rep(1, ncol(data.train %>% select(starts_with("x")))), 
                        Q.model = Q.model, maxlen = size, seed = NULL, sl.bundle = sl.bundle)
      }, error = function(error_message){
        message("\nHere's the original error message:")
        message(error_message)
        return(NA)
      } 
    )
  } else if (model == "krr"){
    mod.train <- krr(x = data.train %>% select(starts_with("x")) %>% as.matrix(), y = data.train$y, group = data.train$a)
  } else if (model == "lasso"){
    data.train$a = as.factor(data.train$a) %>% relevel(ref = "0")
    
    # Create interactions
    f <- as.formula(~ x1*a + x2*a + x3*a)
    x.train <- model.matrix(f, data.train)
    y.train <- as.matrix(data.train$y, ncol=1)
    
    # Training
    mod.train <- glmnet(x.train, y.train, alpha = 1, family = "gaussian")
  } else if (model == "ridge"){
    data.train$a = as.factor(data.train$a) %>% relevel(ref = "0")
    
    # Create interactions
    f <- as.formula(~ x1*a + x2*a + x3*a)
    x.train <- model.matrix(f, data.train)
    y.train <- as.matrix(data.train$y, ncol=1)
    
    # Training
    mod.train <- glmnet(x.train, y.train, alpha = 0, family = "gaussian")
  } else if (model == "elastic_net"){
    data.train$a = as.factor(data.train$a) %>% relevel(ref = "0")
    
    # Create interactions
    f <- as.formula(~ x1*a + x2*a + x3*a)
    x.train <- model.matrix(f, data.train)
    y.train <- as.matrix(data.train$y, ncol=1)
    
    # Training
    mod.train <- glmnet(x.train, y.train, alpha = 0.5, family = "gaussian")
  } else if (model == "rwl"){
    
    x.train <- data.train %>% select(starts_with("x"))
    
    if (num_trt == 2){
      data.train$a = ifelse(data.train$a == 0, -1, data.train$a) # if 0 then -1, if 1 then 1
      mod.train=tryCatch(
        {runRWL(data.train, x.train, data.train$y, "a", x.train %>% colnames(), kernel)
        }, error = function(error_message){
          message("\nHere's the original error message:")
          message(error_message)
          return(NA)
        } 
      )
    } 
    if (num_trt == 3){
      # Recode treatment group for two stages
      training <- data.train %>% mutate(group1 = ifelse(a == 0, 0, 1), 
                                        group2 = ifelse(a == 1, 0, ifelse(a == 2, 1, NA)))
      rwl.s1 <- runRWL(training, x.train, training$y, "group1", x.train %>% colnames(), kernel)
      mod.train = rwl.s1 # model from the first stage
    }
    
  } else if (model == "bart"){
    # BART model needs both training and testing data, so both training and testing can be found in Testing() function 
    # mod.train <- wbart(x.train = data.train %>% select(starts_with("x")), 
    #                   y.train = data.train$y,
    #                   nskip = 500L, ndpost = 5000L, ntree = 500L) 
    mod.train <- NULL
  } else {cat("\nError: Invalid model name!")}
  return(mod.train)
}
  

## Estimate value function (U, W) by jackknife training and testing samples
Get_UW_by_JK <- function(data, model, size, empirical, test, test.data, kernel, Q.model, sl.bundle){
  
  # empirical -- whether or not d.hat depends on all n data or jackknife (n-1) samples for training
  # test -- whether test data is the ith jackknife sample or a new independent copy
  # sl.bundle -- extra parameters needed by the super learner model; NULL for other models
  
  Ui.vec <- NULL
  Wi.vec <- NULL
  A.hat.vec <- NULL
  
  num_trt <- length(unique(data$a))

  #### Training (empirical) ####
  # To save computation time, split jackknife and empirical training
  # Fit model once for empirical estimator because training data do not change by i
  if (empirical == T){
    data.train = data
    mod.train <- Training(data.train, model, size, kernel, Q.model, sl.bundle)
  }
  
  #### Start the loop for each jackknife sample ####
  # Find out which dataset is the test set
  if (class(test.data) == "list") {
    which_test = test.data$pop.1
    which_test %<>% rename(A = a, a = A) 
  } else if (is.null(test.data)) {
    which_test = data
  } else {which_test = test.data}
  # if test is population; if test is null, no honesty; if test is another independent sample, honesty
  
  if (!is.data.frame(which_test)){which_test = as.data.frame(which_test)}
  
  # s is the number of rows a test set has to go through in a loop
  s <- nrow(which_test)
  
  # Generate propensity scores 
  prop_scores <- which_test %>% dplyr::group_by(a) %>% dplyr::summarise(count = n() / nrow(which_test)) 
  
  prop_scores$a <- ifelse(prop_scores$a == -1, 0, prop_scores$a) # convert back to 0 and 1 for propensity scores
  
  # Saving optTx results for stage-1 RWL models
  optTx.s1.all <- NULL 
  
  # Loop for each testing sample
  for (i in 1:s){
    #### Training (jackknife) ####
    # Fit model (should be the same code as the empirical == T but the data.train is different)
    if (empirical == F){
      data.train = data[-i,] # n-1 x 9
      sl.bundle[[6]] <- i
      mod.train <- Training(data.train, model, size, kernel, Q.model, sl.bundle)
    } 
    
    # Define test sets, # 1 x 9 
    if (test == F){
      data.test = data[i,]
      population = F
    } 
    if (test == T){
      if (class(test.data) == "list"){
        population = T # Indicator that the test data is population data
      }
      if (class(test.data) == "data.frame"){
        population = F 
        data.test = test.data[i,]
      }
    } 
    
    #### Testing #####
    # Get estimated decision rule with test data and UW
    if (population == F){
      if (model == "rf"){
        rf.preds <- NULL
        for (l in 0:(num_trt - 1)){
          rf.pred <- predict(mod.train, data.test %>% mutate(a = l) %>% as.matrix())
          rf.preds <- c(rf.preds, rf.pred)
        }
        rf.preds <- rf.preds %>% t() %>% as.data.frame()
        colnames(rf.preds) <- 0:(num_trt-1)
        if (length(which(rf.preds == max(rf.preds))) == 1){
          A.hat <- which.max(rf.preds) %>% names() %>% as.numeric()
          Ui <- Get_Each_UW(data.test$y, data.test$a, A.hat, "U", prop_scores, num_trt)
          Wi <- Get_Each_UW(data.test$y, data.test$a, A.hat, "W", prop_scores, num_trt)
        } else {Ui = NA; Wi = NA; A.hat = NA; print("Tied rf.preds!");print(rf.preds)}
      } else if (model == "rlt") {
        rlt.preds <- NULL
        for (l in 0:(num_trt - 1)){
          rlt.pred <- predict(mod.train, data.test %>% select(starts_with("x")) %>% mutate(a = l))$Prediction
          rlt.preds <- c(rlt.preds, rlt.pred)
        }
        rlt.preds <- rlt.preds %>% t() %>% as.data.frame()
        colnames(rlt.preds) <- 0:(num_trt - 1)
        A.hat <- which.max(rlt.preds) %>% names() %>% as.numeric()
        Ui <- Get_Each_UW(data.test$y, data.test$a, A.hat, "U", prop_scores, num_trt)
        Wi <- Get_Each_UW(data.test$y, data.test$a, A.hat, "W", prop_scores, num_trt)
      } else if (model == "list"){
        UW = tryCatch(
          {A.hat <- predict(mod.train, stage = 1, 
                            xnew = data.test %>% dplyr::select(starts_with("x")) %>% as.matrix()) %>% unfactor()
          Ui <- Get_Each_UW(data.test$y, data.test$a, A.hat, "U", prop_scores, num_trt)
          Wi <- Get_Each_UW(data.test$y, data.test$a, A.hat, "W", prop_scores, num_trt)
          UW <- list(Ui, Wi, A.hat)
          }, error = function(error_message){
            message("\nHere's the original error message:")
            message(error_message)
            return(list(NA, NA, NA))}
        )
        Ui = UW[[1]]
        Wi = UW[[2]]
        A.hat = UW[[3]]
      } else if (model == "krr"){
        krr.pred <- predict(mod.train, data.test %>% as.matrix())
        A.hat <- colnames(krr.pred)[which.max(krr.pred)] %>% as.numeric() 
        Ui <- Get_Each_UW(data.test$y, data.test$a, A.hat, "U", prop_scores, num_trt)
        Wi <- Get_Each_UW(data.test$y, data.test$a, A.hat, "W", prop_scores, num_trt)
      } else if (model == "lasso"){
        f <- as.formula(~ x1*a + x2*a + x3*a)
        
        # Create three test sets
        data.all <- rbind(data.test, data.train)
        data.all$a = as.factor(data.all$a) %>% relevel(ref = "0")
        new.train = data.all[-1,]
        x.train <- model.matrix(f, new.train)
        y.train <- as.matrix(new.train$y, ncol=1)
        
        x.test = suppressWarnings(data.frame(data.all[1,], freq = 1:num_trt) %>% 
                                    mutate(a = as.factor(0:(num_trt-1))) %>% 
                                    model.matrix(f, .))
        
        # Testing
        cv = cv.glmnet(x.train, y.train, type.measure = "mse", nfolds = 5)
        lasso.preds <- predict(mod.train, x.test, s = cv$lambda.min)
        lasso.preds %<>% t() %>% as.data.frame()
        colnames(lasso.preds) <- 0:(num_trt - 1)
        A.hat <- which.max(lasso.preds) %>% names() %>% as.numeric()
        Ui <- Get_Each_UW(data.test$y, data.test$a, A.hat, "U", prop_scores, num_trt)
        Wi <- Get_Each_UW(data.test$y, data.test$a, A.hat, "W", prop_scores, num_trt)
      } else if (model == "ridge"){
        f <- as.formula(~ x1*a + x2*a + x3*a)
        
        # Create three test sets
        data.all <- rbind(data.test, data.train)
        data.all$a = as.factor(data.all$a) %>% relevel(ref = "0")
        new.train = data.all[-1,]
        x.train <- model.matrix(f, new.train)
        y.train <- as.matrix(new.train$y, ncol=1)
        
        x.test = suppressWarnings(data.frame(data.all[1,], freq = 1:num_trt) %>% 
                                    mutate(a = as.factor(0:(num_trt-1))) %>% 
                                    model.matrix(f, .))
        
        # Testing
        cv = cv.glmnet(x.train, y.train, type.measure = "mse", nfolds = 5)
        ridge.preds <- predict(mod.train, x.test, s = cv$lambda.min)
        ridge.preds %<>% t() %>% as.data.frame()
        colnames(ridge.preds) <- 0:(num_trt - 1)
        A.hat <- which.max(ridge.preds) %>% names() %>% as.numeric()
        Ui <- Get_Each_UW(data.test$y, data.test$a, A.hat, "U", prop_scores, num_trt)
        Wi <- Get_Each_UW(data.test$y, data.test$a, A.hat, "W", prop_scores, num_trt)
      } else if (model == "elastic_net"){
        f <- as.formula(~ x1*a + x2*a + x3*a)
        
        # Create three test sets
        data.all <- rbind(data.test, data.train)
        data.all$a = as.factor(data.all$a) %>% relevel(ref = "0")
        new.train = data.all[-1,]
        x.train <- model.matrix(f, new.train)
        y.train <- as.matrix(new.train$y, ncol=1)
        
        x.test = suppressWarnings(data.frame(data.all[1,], freq = 1:num_trt) %>% 
                                    mutate(a = as.factor(0:(num_trt-1))) %>% 
                                    model.matrix(f, .))
        
        # Testing
        cv = cv.glmnet(x.train, y.train, type.measure = "mse", nfolds = 5)
        en.preds <- predict(mod.train, x.test, s = cv$lambda.min)
        en.preds %<>% t() %>% as.data.frame()
        colnames(en.preds) <- 0:(num_trt - 1)
        A.hat <- which.max(en.preds) %>% names() %>% as.numeric()
        Ui <- Get_Each_UW(data.test$y, data.test$a, A.hat, "U", prop_scores, num_trt)
        Wi <- Get_Each_UW(data.test$y, data.test$a, A.hat, "W", prop_scores, num_trt)
      } else if (model == "bart"){
        # Make and combine x.test data
        test1 = data.test %>% mutate(a = 0)
        test2 = data.test %>% mutate(a = 1)
        test3 = data.test %>% mutate(a = 2)
        x.test.new = rbind(test1, test2, test3)
        # BART for continuous outcomes, posterior
        bart.mod <- wbart(x.train =  data.train %>% select(starts_with("x")), 
                          y.train = data.train$y,
                          x.test = x.test.new, 
                          nskip = 500L, ndpost = 5000L, ntree = 500L) 
        # nskip is number of MCMC draws to burn in, ndpost is number of MCMC samples to keep, 
        # default settings gave always models, so trying different options (increased nskip, ndpost, ntree)
        
        # Find the BART ITR for test set
        N <- nrow(data.test)
        pred <- bart.mod$yhat.test.mean 
        if (nrow(x.test.new) == 1){ # jackknife
          bart.preds = pred %>% as.matrix() %>% t()
        } else {
          bart.preds = cbind(pred[(1:N)], pred[N+(1:N)], pred[2*N+(1:N)])
        }
        colnames(bart.preds) = 0:(num_trt - 1)
        itr.pick <- integer(N)
        for(i in 1:N) itr.pick[i] <- which.max(bart.preds[i,]) %>% names() %>% as.numeric()
        A.hat <- itr.pick
        Ui <- Get_Each_UW(data.test$y, data.test$a, A.hat, "U", prop_scores, num_trt)
        Wi <- Get_Each_UW(data.test$y, data.test$a, A.hat, "W", prop_scores, num_trt)
        
      } else if (model == "rwl"){
        
        if (num_trt == 2){
          UW = tryCatch(
            {rwl.pred.ori <- optTx(mod.train, data.test %>% dplyr::select(starts_with("x")))$optimalTx # should be -1 or 1
            A.hat <- ifelse(rwl.pred.ori == -1, 0, 1) # should be 0 or 1
            Ui <- Get_Each_UW(data.test$y, data.test$a, A.hat, "U", prop_scores, num_trt)
            Wi <- Get_Each_UW(data.test$y, data.test$a, A.hat, "W", prop_scores, num_trt)
            UW <- list(Ui, Wi, A.hat)
            }, error = function(error_message){
              message("\nHere's the original error message:")
              message(error_message)
              return(list(NA, NA, NA))}
          )
        }
        if (num_trt == 3){
          pred.s1 <- optTx(mod.train, data.test %>% dplyr::select(starts_with("x")))$optimalTx # should be -1 or 1
          optTx.s1 <- data.frame(pid = i, a = data.test$a, opt_group = pred.s1, obs_outcome = data.test$y)
          optTx.s1.all <- rbind(optTx.s1.all, optTx.s1)
          UW <- list(NA, NA, NA)
        }
        
        Ui = UW[[1]]
        Wi = UW[[2]]
        A.hat = UW[[3]]
        
      } else if (grepl("^zero_order[0-9]$", model)) {
        which_zero <- model %>% str_match("[0-9]+") %>% unlist() %>% as.numeric()
        A.hat = which_zero # zero_order d.hat
        Ui <- Get_Each_UW(data.test$y, data.test$a, A.hat, "U", prop_scores, num_trt)
        Wi <- Get_Each_UW(data.test$y, data.test$a, A.hat, "W", prop_scores, num_trt)
      } else {cat("\nError: Invalid model name!")}
      
      Ui.vec <- c(Ui.vec, Ui)
      Wi.vec <- c(Wi.vec, Wi)
      A.hat.vec <- c(A.hat.vec, A.hat)
    }#end of if population == F 
  }#i as removal of the ith subject
  
  ## Stage 2 of RWL models
  if (num_trt == 3 & model == "rwl"){
    
    # Part 1, for those whose optimal group and observed group are both the combined group (1+2)
    selected.p1 <- optTx.s1.all %>% filter(a %in% c(1,2), opt_group == 1)
    optTx.s2.p1.all <- NULL
    if (nrow(selected.p1) != 0){
      dat.s2.p1 <- data %>% mutate(pid = 1:nrow(data)) %>% filter(pid %in% selected.p1$pid)
      # Loop through each test set of the subset of those whose optimal group and observed group are both the combined group (1+2)
      for (i in 1:nrow(dat.s2.p1)){
        dat.train.s2.p1 = dat.s2.p1[-i,]
        dat.test.s2.p1 = dat.s2.p1[i,]
        x.train.s2.p1 = dat.train.s2.p1 %>% dplyr::select(contains("x"))
        x.test.s2.p1 = dat.test.s2.p1 %>% dplyr::select(contains("x"))
        
        optTx.s2.p1 = tryCatch(
          {
            # Recode treatment group for two stages
            training <- dat.train.s2.p1 %>% mutate(group1 = ifelse(a == 0, 0, 1), 
                                              group2 = ifelse(a == 1, 0, ifelse(a == 2, 1, NA)))
            rwl.s2.p1 <- runRWL(training, x.train.s2.p1, training$y, "group2", x.train.s2.p1 %>% colnames(), kernel)
            pred.s2.p1 <- optTx(rwl.s2.p1, x.test.s2.p1 %>% dplyr::select(starts_with("x")))$optimalTx
            optTx.s2.p1 <- data.frame(pid = dat.test.s2.p1$pid, a = dat.test.s2.p1$a, opt_group = pred.s2.p1, obs_outcome = dat.test.s2.p1$y)
            optTx.s2.p1
          }, error = function(error_message){
            message("\nHere's the original error message:")
            message(error_message)
            return(data.frame(pid = dat.test.s2.p1$pid, a = dat.test.s2.p1$a, opt_group = NA, obs_outcome = dat.test.s2.p1$y, row.names = "y"))
          }
        )
        optTx.s2.p1.all <- rbind(optTx.s2.p1.all, optTx.s2.p1)
      }#i for each jackknife test set in RWL Stage 2, Part 1
    }# ifelse nrow(selected.p1) != 0
    
    # Part 2, for those whose optimal group is 1+2 but observed group is 0 (test set)
    # No jackknife for this because nothing to learn from 1,2 when don't observe 1,2
    # Just a lump sum prediction on all test set using s2.p1 data
    selected.p2 <- optTx.s1.all %>% filter(a == 0, opt_group == 1)
    optTx.s2.p2 <- NULL
    
    if (nrow(selected.p2) != 0){
      dat.s2.p2 <- data %>% mutate(pid = 1:nrow(data)) %>% filter(pid %in% selected.p2$pid)
      x.train.s2.p2 <- dat.s2.p1 %>% select(contains("x"))
      x.test.s2.p2 <- dat.s2.p2 %>% select(contains("x"))
      
      optTx.s2.p2 = tryCatch(
        { # Recode treatment group for two stages
          training <- dat.s2.p1 %>% mutate(group1 = ifelse(a == 0, 0, 1), 
                                                 group2 = ifelse(a == 1, 0, ifelse(a == 2, 1, NA)))
          rwl.s2.p2 <- runRWL(training, x.train.s2.p2, training$y, "group2", x.train.s2.p2 %>% colnames(), kernel)
          pred.s2.p2 <- optTx(rwl.s2.p2, x.test.s2.p2 %>% dplyr::select(starts_with("x")))$optimalTx 
          optTx.s2.p2 <- data.frame(pid = dat.s2.p2$pid, a = dat.s2.p2$a, opt_group = pred.s2.p2, obs_outcome = dat.s2.p2$y)
          optTx.s2.p2
        }, error = function(error_message){
          message("\nHere's the original error message:")
          message(error_message)
          return(data.frame(pid = dat.s2.p2$pid, a = dat.s2.p2$a, opt_group = NA, obs_outcome = dat.s2.p2$y))
        }
      )
    }# ifelse nrow(selected.p2) != 0
    
    # Put all optimal group predictions together
    optTx.s1.all$opt_group <- ifelse(optTx.s1.all$opt_group == 0, 0, 
                                     ifelse(optTx.s1.all$opt_group == 1, 12, NA)) # group 1
    optTx.s2.p1.all$opt_group <- ifelse(optTx.s2.p1.all$opt_group == 0, 1,
                                        ifelse(optTx.s2.p1.all$opt_group == 1, 2, NA)) # group 2
    optTx.s2.p2$opt_group <- ifelse(optTx.s2.p2$opt_group == 0, 1,
                                    ifelse(optTx.s2.p2$opt_group == 1, 2, NA)) # group 2
    optTx_result <- rbind(optTx.s1.all %>% filter(opt_group == 0),
                          optTx.s2.p1.all,
                          optTx.s2.p2)
    optTx_results.all <- optTx_result %>% select(pid, a, opt_group, obs_outcome) %>% arrange(pid)
    stopifnot(optTx_results.all$pid == 1:nrow(data))
    stopifnot(optTx_results.all$obs_outcome == data$y)
    
    Ui.vec <- apply(optTx_results.all, 1, function(x){Get_Each_UW(x['obs_outcome'], x['a'], x['opt_group'], "U", prop_scores, num_trt)})
    Wi.vec <- apply(optTx_results.all, 1, function(x){Get_Each_UW(x['obs_outcome'], x['a'], x['opt_group'], "W", prop_scores, num_trt)})
    A.hat.vec <- optTx_results.all$opt_group
  }
  
  # When test set is population data
  if (population == T){
    if (model == "rf"){
      rf.preds <- NULL
      for (l in 1:num_trt){
        data.pop.each <- test.data[[l]]
        rf.pred <- predict(mod.train, data.pop.each %>% select(starts_with("x"), a) %>% as.matrix()) # no need to mutate a to all 0 or all 1 because the a in population is already so
        rf.preds <- cbind(rf.preds, rf.pred)
      }#l 
      colnames(rf.preds) <- sub("pop.", "", names(test.data))
      A.hat.vec <- colnames(rf.preds)[apply(rf.preds, 1, which.max)] %>% as.numeric()
      stopifnot(A.hat.vec %in% 0:(num_trt-1))
      Ui.vec <- ifelse(A.hat.vec == 1, test.data[[2]]$Y1, ifelse(A.hat.vec == 2, test.data[[3]]$Y2, test.data[[1]]$Y0))
      Wi.vec <- rep(NA, s)
    } else if (model == "rlt") {
      rlt.preds <- NULL
      for (l in 1:num_trt){
        data.pop.each <- test.data[[l]]
        rlt.pred <- predict(mod.train, data.pop.each %>% select(starts_with("x"), a) %>% as.matrix())$Prediction
        rlt.preds <- cbind(rlt.preds, rlt.pred)
      }#l
      colnames(rlt.preds) <- sub("pop.", "", names(test.data))
      A.hat.vec <- colnames(rlt.preds)[apply(rlt.preds, 1, which.max)] %>% as.numeric()
      stopifnot(A.hat.vec %in% 0:(num_trt-1))
      Ui.vec <- ifelse(A.hat.vec == 1, test.data[[2]]$Y1, ifelse(A.hat.vec == 2, test.data[[3]]$Y2, test.data[[1]]$Y0))
      Wi.vec <- rep(NA, s)
    }  else if (model == "list"){
      UW = tryCatch(
        {A.hat.vec <- predict(mod.train, stage = 1, xnew = test.data[[1]] %>% dplyr::select(starts_with("x")) %>% as.matrix()) %>% unfactor()
        stopifnot(A.hat.vec %in% 0:(num_trt-1))
        Ui.vec <- ifelse(A.hat.vec == 1, test.data[[2]]$Y1, ifelse(A.hat.vec == 2, test.data[[3]]$Y2, test.data[[1]]$Y0))
        Wi.vec <- rep(NA,s)
        UW <- list(Ui.vec, Wi.vec, A.hat.vec)
        }, error = function(error_message){
          message("\nHere's the original error message:")
          message(error_message)
          return(list(rep(NA, s), rep(NA, s), rep(NA, s)))}
      )
      Ui.vec = UW[[1]]
      Wi.vec = UW[[2]]
      A.hat.vec = UW[[3]]
    } else if (model == "krr"){
      krr.pred <- predict(mod.train, test.data[[1]] %>% as.matrix())
      A.hat.vec <- colnames(krr.pred)[apply(krr.pred, 1, which.max)] %>% as.numeric()
      stopifnot(A.hat.vec %in% 0:(num_trt-1))
      Ui.vec <- ifelse(A.hat.vec == 1, test.data[[2]]$Y1, ifelse(A.hat.vec == 2, test.data[[3]]$Y2, test.data[[1]]$Y0))
      Wi.vec <- rep(NA, s)
    } else if (model == "lasso"){
      f <- as.formula(~ x1*a + x2*a + x3*a)
      
      # Create three test sets
      new.train <- data.train
      new.train$a = as.factor(new.train$a) %>% relevel(ref = "0")
      x.train <- model.matrix(f, new.train)
      y.train <- as.matrix(new.train$y, ncol=1)
      
      r <- nrow(test.data[[1]])
      
      x.test = suppressWarnings(rbind(test.data[[1]] %>% select(starts_with("x")) %>% mutate(a = 0), 
                                      test.data[[2]] %>% select(starts_with("x")) %>% mutate(a = 1), 
                                      test.data[[3]] %>% select(starts_with("x")) %>% mutate(a = 2)) %>% 
                                  as.data.frame() %>% 
                                  mutate(a = as.factor(a)) %>% 
                                  model.matrix(f, .))
      
      # Testing
      cv = cv.glmnet(x.train, y.train, type.measure = "mse", nfolds = 5)
      outcomes <- predict(mod.train, x.test, s = cv$lambda.min)
      lasso.preds <- cbind(outcomes[1:r,], outcomes[(r+1):(2*r),], outcomes[(2*r + 1): (3*r),])
      colnames(lasso.preds) <- sub("pop.", "", names(test.data))
      A.hat.vec <- colnames(lasso.preds)[apply(lasso.preds, 1, which.max)] %>% as.numeric()
      stopifnot(A.hat.vec %in% 0:(num_trt-1))
      Ui.vec <- ifelse(A.hat.vec == 1, test.data[[2]]$Y1, ifelse(A.hat.vec == 2, test.data[[3]]$Y2, test.data[[1]]$Y0))
      Wi.vec <- rep(NA, s)
    } else if (model == "ridge"){
      f <- as.formula(~ x1*a + x2*a + x3*a)
      
      # Create three test sets
      new.train <- data.train
      new.train$a = as.factor(new.train$a) %>% relevel(ref = "0")
      x.train <- model.matrix(f, new.train)
      y.train <- as.matrix(new.train$y, ncol=1)
      
      r <- nrow(test.data[[1]])
      
      x.test = suppressWarnings(rbind(test.data[[1]] %>% select(starts_with("x")) %>% mutate(a = 0), 
                                      test.data[[2]] %>% select(starts_with("x")) %>% mutate(a = 1), 
                                      test.data[[3]] %>% select(starts_with("x")) %>% mutate(a = 2)) %>% 
                                  as.data.frame() %>% 
                                  mutate(a = as.factor(a)) %>% 
                                  model.matrix(f, .))
      
      # Testing
      cv = cv.glmnet(x.train, y.train, type.measure = "mse", nfolds = 5)
      outcomes <- predict(mod.train, x.test, s = cv$lambda.min)
      ridge.preds <- cbind(outcomes[1:r,], outcomes[(r+1):(2*r),], outcomes[(2*r + 1): (3*r),])
      colnames(ridge.preds) <- sub("pop.", "", names(test.data))
      A.hat.vec <- colnames(ridge.preds)[apply(ridge.preds, 1, which.max)] %>% as.numeric()
      stopifnot(A.hat.vec %in% 0:(num_trt-1))
      Ui.vec <- ifelse(A.hat.vec == 1, test.data[[2]]$Y1, ifelse(A.hat.vec == 2, test.data[[3]]$Y2, test.data[[1]]$Y0))
      Wi.vec <- rep(NA, s)
      
    } else if (model == "elastic_net"){
      f <- as.formula(~ x1*a + x2*a + x3*a)
      
      # Create three test sets
      new.train <- data.train
      new.train$a = as.factor(new.train$a) %>% relevel(ref = "0")
      x.train <- model.matrix(f, new.train)
      y.train <- as.matrix(new.train$y, ncol=1)
      
      r <- nrow(test.data[[1]])
      
      x.test = suppressWarnings(rbind(test.data[[1]] %>% select(starts_with("x")) %>% mutate(a = 0), 
                                      test.data[[2]] %>% select(starts_with("x")) %>% mutate(a = 1), 
                                      test.data[[3]] %>% select(starts_with("x")) %>% mutate(a = 2)) %>% 
                                  as.data.frame() %>% 
                                  mutate(a = as.factor(a)) %>% 
                                  model.matrix(f, .))
      
      # Testing
      cv = cv.glmnet(x.train, y.train, type.measure = "mse", nfolds = 5)
      outcomes <- predict(mod.train, x.test, s = cv$lambda.min)
      en.preds <- cbind(outcomes[1:r,], outcomes[(r+1):(2*r),], outcomes[(2*r + 1): (3*r),])
      colnames(en.preds) <- sub("pop.", "", names(test.data))
      A.hat.vec <- colnames(en.preds)[apply(en.preds, 1, which.max)] %>% as.numeric()
      stopifnot(A.hat.vec %in% 0:(num_trt-1))
      Ui.vec <- ifelse(A.hat.vec == 1, test.data[[2]]$Y1, ifelse(A.hat.vec == 2, test.data[[3]]$Y2, test.data[[1]]$Y0))
      Wi.vec <- rep(NA, s)
      
    } else if (model == "rwl"){
      
      
      UW = tryCatch(
        {rwl.pred.ori <- optTx(mod.train, test.data[[1]] %>% dplyr::select(starts_with("x")))$optimalTx # should be -1 or 1
        A.hat.vec <- ifelse(rwl.pred.ori == -1, 0, 1) # should be 0 or 1
        stopifnot(A.hat.vec %in% 0:(num_trt-1))
        Ui.vec <- ifelse(A.hat.vec == 1, test.data[[2]]$Y1, ifelse(A.hat.vec == 2, test.data[[3]]$Y2, test.data[[1]]$Y0))
        Wi.vec <- rep(NA, s)
        UW <- list(Ui.vec, Wi.vec, A.hat.vec)
        }, error = function(error_message){
          message("\nHere's the original error message:")
          message(error_message)
          return(list(rep(NA,s), rep(NA,s), rep(NA,s)))}
      )
      Ui.vec = UW[[1]]
      Wi.vec = UW[[2]]
      A.hat.vec = UW[[3]]
    } else if (grepl("^zero_order[0-9]$", model)){
      A.hat.vec <- rep(which_zero, s)
      Ui.vec <- ifelse(A.hat.vec == 1, test.data[[2]]$Y1, ifelse(A.hat.vec == 2, test.data[[3]]$Y2, test.data[[1]]$Y0))
      Wi.vec <- rep(NA,s)
    } else {cat("\nError: Invalid model name!")}
  }# end of if population == T
  
  return(list(U = Ui.vec, W = Wi.vec, A = A.hat.vec))
}

## F3_ZOM.R
Get_UW_by_JK_zero_order <- function(data, which_zero, empirical, test, test.data){
  # "empirical" means whether or not d.hat depends on all n data or jackknife (n-1) samples for training
  # "test" means whether test data is the ith jackknife sample or a new independent copy
  Ui.vec <- NULL
  Wi.vec <- NULL
  A.hat.vec <- NULL
  num_trt <- length(unique(data$a))
    
  #### Start the loop for each jackknife sample ####
  # Find out which dataset is the test set
  if (class(test.data) == "list") {
    which_test = test.data$pop.1
    which_test %<>% rename(A = a, a = A) 
  } else if (is.null(test.data)) {
    which_test = data
  } else {which_test = test.data}
  # if test is population; if test is null, no honesty; if test is another independent sample, honesty
  
  if (!is.data.frame(which_test)){which_test = as.data.frame(which_test)}
  
  # s is the number of rows a test set has to go through in a loop
  s <- nrow(which_test)
  
  # Generate propensity scores 
  prop_scores <- which_test %>% dplyr::group_by(a) %>% dplyr::summarise(count = n() / nrow(which_test)) 
  
  prop_scores$a <- ifelse(prop_scores$a == -1, 0, prop_scores$a) # convert back to 0 and 1 for propensity scores
  
  # Define test sets 
  if (test == F){
    data.test = data
    population = F
  } 
  if (test == T){
    if (class(test.data) == "list"){
      population = T # Indicator that the test data is population data
    }
    if (class(test.data) == "data.frame"){
      population = F 
      data.test = test.data
    }
  } 
  
  A.hat.vec = rep(which_zero, s) # zero_order d.hat
  
  # Get estimated decision rule with test data and UW
  if (population == F){
    stopifnot(unique(data.test$a) %in% 0:num_trt)
    
    data.test %<>% mutate(indicator = ifelse(which_zero == a, 1, 0)) %>% 
      left_join(prop_scores, by = "a")
    
    Wi.vec <- data.test$indicator / data.test$count
    Ui.vec <- data.test$y * data.test$indicator / data.test$count
    
    # Ui.vec <- Get_Each_UW(data.test$y, data.test$a, A.hat.vec, "U", prop_scores, num_trt)
    # Wi.vec <- Get_Each_UW(data.test$y, data.test$a, A.hat.vec, "W", prop_scores, num_trt)
  }
  
  if (population == T){
    Ui.vec <- ifelse(A.hat.vec == 1, test.data[[2]]$Y1, ifelse(A.hat.vec == 2, test.data[[3]]$Y2, test.data[[1]]$Y0))
    Wi.vec <- rep(NA,s)
  }
  
  return(list(U = Ui.vec, W = Wi.vec, A = A.hat.vec))
}


## Get estimated policy, d.hat(X.pop), with d.hat trained on all X.train; output would be 1 x n.pop
Get_Z.hat <- function(data, model, size, test, propensity){
  if (model == "rf"){
    rf.formula <- as.formula(paste0("y~a+", paste0(colnames(data)[colnames(data) %like% "x"], collapse = "+")))
    rf.mod <- randomForest(rf.formula, data = data)
    y1.hat <- predict(rf.mod, test %>% mutate(a=1))
    y0.hat <- predict(rf.mod, test %>% mutate(a=0))
    d.hat <- (y1.hat >= y0.hat) # n.test x 1
  } else if (model == "krr"){
    krr.mod <- krr(x = data %>% dplyr::select(starts_with("x")) %>% as.matrix(), y = data$y, group = data$a)
    y.hat <- predict(krr.mod, test %>% dplyr::select(starts_with("x")) %>% as.matrix())
    d.hat <- (y.hat[,2] >= y.hat[,1]) # n.test x 1
  } else if (model == "list"){
    d.hat = tryCatch(
      {list.mod <- listdtr(y = data$y, a = data$a, x = data %>% dplyr::select(starts_with("x")) %>% as.matrix, stage.x = rep(1,ncol(data %>% dplyr::select(starts_with("x")))), maxlen = size)
      d.hat <- predict(list.mod, stage = 1, xnew = test %>% dplyr::select(starts_with("x")) %>% as.matrix()) %>% unfactor()
      }, error = function(error_message){
        message(error_message)
        d.hat <- NA
      }
    )
  } else {print("ERROR: wrong model name!")}
  Z.hat <- mean(test$y * ifelse(test$a == d.hat, 1, 0)) / propensity # scalar, sum over n.test
  return(Z.hat)
}

## Get estimated policy, d.tilde.hat(X.pop), using d.tilde.hat trained on all X except X.tilde.i; output would be n x n.pop
Get_Z.tilde.hat <- function(data, model, size, test, seed_, error.sd_, propensity){
  data.tilde <- data.generator(n_ = n, propensity_ = propensity, type_ = data.type, error.sd_ = error.sd_, seed_ = seed_) 
  d.tilde.hat <- NULL
  
  for (i in 1:nrow(data)){
    data.new <- data
    data.new[i, grepl( "x", colnames(data) )] <- data.tilde[i, grepl("x", colnames(data))] # replace the ith patient's all dimensions of X with the independent replicate X tilde
    if (model == "rf"){
      rf.formula <- as.formula(paste0("y~a+", paste0(colnames(data.new)[colnames(data.new) %like% "x"], collapse = "+")))
      rf.mod <- randomForest(rf.formula, data = data.new)
      y1.hat <- predict(rf.mod, test %>% mutate(a=1))
      y0.hat <- predict(rf.mod, test %>% mutate(a=0))
      d.tilde.hat.i <- (y1.hat >= y0.hat) # n.test x 1
    } else if (model == "krr"){
      krr.mod <- krr(x = data.new %>% dplyr::select(starts_with("x")) %>% as.matrix(), y = data.new$y, group = data.new$a)
      y.hat <- predict(krr.mod, test %>% dplyr::select(starts_with("x")) %>% as.matrix())
      d.tilde.hat.i <- (y.hat[,2] >= y.hat[,1]) # n.test x 1
    } else if (model == "list"){
      d.tilde.hat.i = tryCatch(
        {list.mod <- listdtr(y = data.new$y, a = data.new$a, x = data.new %>% dplyr::select(starts_with("x")) %>% as.matrix, stage.x = rep(1,ncol(data.new %>% dplyr::select(starts_with("x")))), maxlen = size)
        d.hat <- predict(list.mod, stage = 1, xnew = pop.1 %>% dplyr::select(starts_with("x")) %>% as.matrix()) %>% unfactor()
        }, error = function(error_message){
          message(error_message)
          d.tilde.hat.i <- NA
        }
      )
    } else {print("ERROR: wrong model name!")}
    d.tilde.hat <- rbind(d.tilde.hat, d.tilde.hat.i) # n.train x n.test
    
  }#i loop for each subject
  
  Z.tilde.hat <- apply(d.tilde.hat, 1, function(x) {mean(test$y * ifelse(test$a == x, 1, 0)) / propensity}) %>% mean() # scalar, sum over n.train and n.test
  return(Z.tilde.hat)
}


### 3. Define utility functions
# randomization seed
seed.fn = function(i, j) {i + j}

# time calculator
tt <- function(s){
  if (s==1) {time.tmp <<- Sys.time() # record time
  } else if (s==2) { # calculate time
    return(data.frame(begin = time.tmp, end = Sys.time(), elapsed = Sys.time() - time.tmp))
  }
}

# Calculate value function 
Get_Value_CV <- function(UW_list, num_folds){
  
  Umk <- UW_list[["Umk"]]
  Wmk <- UW_list[["Wmk"]]
  Ubar <- mean(Umk)
  Wbar <- mean(Wmk)
  n <- length(Umk)  # n = mk = 500
  stopifnot(length(Umk) == length(Wmk))
  
  mean_value <- Ubar / Wbar # or equivalently, sum(Umk) / sum(Wmk)
  
  Rmk <- Umk / Wbar - Wmk * Ubar / (Wbar^2) 
  se_value <- sqrt(sum(Rmk^2) / (num_folds*(n-1)))
  
  return(list(mean_value = mean_value, se_value = se_value))
}

Get_Rmk_CV <- function(UW_list){
  Umk <- UW_list[["Umk"]]
  Wmk <- UW_list[["Wmk"]]
  
  Ubar <- mean(Umk, na.rm = T)
  Wbar <- mean(Wmk, na.rm = T)
  
  Rmk <- Umk / Wbar - Wmk * Ubar / (Wbar^2) 
  
  return(Rmk)
}