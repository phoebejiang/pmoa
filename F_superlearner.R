#####
# PROJECT. Precision Medicine for Knee OA
# AUTHOR. Xiaotong Jiang
#####
# USAGE. To provide utility functions for the LIST super learning model. Run main.R for PMOA.
#####

# Generate simplex of length n
get.simplex = function(n){
  alpha = rep(0, n)
  sum = 0
  for(i in 1:(n-1)){
    alpha[i] = runif(1, 0, 1 - sum)
    sum = sum + alpha[i]
    # cat("alpha: ", alpha[i], "sum", sum, "\n")
  }
  alpha[n] = 1 - sum
  
  return (sample(alpha))
}


## Base learners

# Model 1: RF
RF <- function(outcome, x.train, x.test, dat.train, dat.test){
  # x.train and x.test, only X variables (not including pid but including a)
  # dat.test has pid that test or x.test does not have
  # train and test have X variables and the outcome needed
  
  train <- cbind(x.train, dat.train[[outcome]]) 
  colnames(train)[ncol(train)] <- outcome
  
  rf.formula <- as.formula(paste0(outcome, "~ ."))
  rf.mod <- randomForest(rf.formula, data = train, importance=TRUE, proximity = TRUE) 
  
  pred1 <- predict(rf.mod, newdata = x.test %>% mutate(a = 0))
  pred2 <- predict(rf.mod, newdata = x.test %>% mutate(a = 1))
  pred3 <- predict(rf.mod, newdata = x.test %>% mutate(a = 2))
  
  pred_result <- data.frame(
    pred1 = pred1, pred2 = pred2, pred3 = pred3, 
    rand = dat.test$a, 
    id = dat.test$pid) 
  
  return(pred_result)
}

# Model 2: RLT (added 110818)
rlt <- function(outcome, x.train, x.test, dat.train, dat.test){
  # x.train and x.test, only X variables (not including rand and pid)
  # dat.test has pid that test or x.test does not have
  # train and test have X variables, rand, and the outcome needed
  
  rlt.mod <- RLT::RLT(x = x.train, y = dat.train[[outcome]], reinforcement = T, use.cores = 1, combsplit = 2, ntrees = 50)
  pred1 <- predict(rlt.mod, x.test %>% mutate(a = 0))$Prediction
  pred2 <- predict(rlt.mod, x.test %>% mutate(a = 1))$Prediction
  pred3 <- predict(rlt.mod, x.test %>% mutate(a = 2))$Prediction
  
  pred_result <- data.frame(
    pred1 = pred1, pred2 = pred2, pred3 = pred3, 
    rand = dat.test$a, 
    id = dat.test$pid)
  
  return(pred_result)
}

# Model 3: Penalized
PENAL <- function(outcome, x.train, x.test, dat.train, dat.test, alpha){
  allX <- colnames(x.train)
  otherX <- allX[!grepl("^a", allX)]
  rand <- setdiff(allX, otherX)
  interactions <- apply(expand.grid(otherX, rand), 1, paste, collapse=":")
  
  f <- as.formula(paste0("~", paste0(rand, collapse = "+"), "+", 
                         paste0(otherX, collapse = "+"), "+",
                         paste0(interactions, collapse = "+")))
  
  new.x.train <- model.matrix(f, x.train) 
  new.y.train <- as.matrix(dat.train[[outcome]], ncol=1)
  
  # Training
  mod = glmnet(new.x.train, new.y.train, alpha = alpha, family = "gaussian")
  cv = cv.glmnet(new.x.train, new.y.train, type.measure = "mse", nfolds = 5)
  
  # 3 possible testings where rand is hypothesized to be 3 treatment groups
  new.test <- NULL
  nrow(x.test) -> r
  for (m in 1:length(unique(c(dat.train$a, dat.test$a)))){
    tmp <- x.test %>% mutate(a = m - 1)
    new.test <- rbind(new.test, tmp)
  }
  
  new.x.test <- model.matrix(f, new.test)
  pred = predict(mod, new.x.test, s = cv$lambda.min) 
  pred.all <- cbind(pred[1:r,], pred[(r+1):(2*r),], pred[(2*r + 1): (3*r),])
  colnames(pred.all) <- c("pred1", "pred2", "pred3")
  
  pred_result <- pred.all %>% 
    as.data.frame() %>% 
    mutate(rand = dat.test$a,
           id = dat.test$pid)
  
  return(pred_result)
}

# Model 4: KRR
KRR <- function(outcome, x.train, x.test, dat.train, dat.test){
  
  if (all(class(x.train) != "matrix")){
    x.train %<>% data.matrix()
    x.test %<>% data.matrix()
  }
  
  krr.mod <- krr(x = x.train[, colnames(x.train) != "a"], y = dat.train[[outcome]] %>% as.numeric, group = dat.train$a)
  
  pred <- predict(krr.mod, x.test[, colnames(x.test) != "a"]) %>% as.data.frame()
  colnames(pred) <- c("pred1", "pred2", "pred3")
  
  pred_result <- pred %>% 
    mutate(rand = dat.test$a, 
           id = dat.test$pid) 
  
  return(pred_result)
}


# Model 5: BART 
BART <- function(outcome, x.train, x.test, dat.train, dat.test){
  
  # Make and combine x.test data
  test1 = x.test %>% mutate(a = 0)
  test2 = x.test %>% mutate(a = 1)
  test3 = x.test %>% mutate(a = 2)
  x.test.new = rbind(test1, test2, test3)
  
  # BART for continuous outcomes, posterior
  bart.mod <- wbart(x.train = x.train, 
                    y.train = dat.train[[outcome]],
                    x.test = x.test.new, 
                    nskip = 500L, ndpost = 5000L, ntree = 500L) 
  # nskip is number of MCMC draws to burn in, ndpost is number of MCMC samples to keep, 
  # default settings sometimes gave always-2 and always-3 models, so trying different options (increased nskip, ndpost, ntree)
  
  # Find the BART ITR for test set
  N <- nrow(x.test)
  pred <- bart.mod$yhat.test.mean 
  if (nrow(x.test.new) == 1){ # jackknife
    pred = pred %>% as.matrix() %>% t()
  } else {
    pred = cbind(pred[(1:N)], pred[N+(1:N)], pred[2*N+(1:N)])
  }
  # itr.pick <- integer(N)
  # for(i in 1:N) itr.pick[i] <- which.max(pred[i,])
  colnames(pred) = c("pred1", "pred2", "pred3")
  
  pred_result <- pred %>% 
    as.data.frame() %>% 
    mutate(rand = dat.test$a, 
           id = dat.test$pid)
  
  return(pred_result)
}



# Objective function for finding alpha's (coefficient of each base learner)
obj.function = function(type, alpha, outcome, pred.files, k, prop_scores){
  
  #' @param type If "max_value", objective function is value function; If "min_mse", objective function is mse loss of prediction
  #' @param pred.files two layers of list for the CV outputs 
  
  if (type == "max_value"){
    J = length(pred.files)
    
    U = 0; W = 0; values = 0
    
    # IDs for each patient (same order)
    ids = pred.files[[1]][[1]]$id
    
    ## For each super learner repetition
    for (v in 1:k){
      # calculate the weighted average of for all methods
      O = subset(pred.files[[1]][[v]], select = -id) * alpha[1]     
      for (j in 2:J){
        O = O + subset(pred.files[[j]][[v]], select = -id) * alpha[j]
      }
      
      # optimal treatment group accoding to O
      opt_group = apply(O, 1, which.max)
      
      # for each individual, calculate value
      for(i in 1:length(ids)){
        this_person = dat %>% filter(pid == ids[i])
        indicator = this_person$rand == opt_group[i]
        weight = prop_scores %>% filter(rand == this_person$rand) 
        U.decision = indicator * this_person[[outcome]] / weight$count
        W.decision = indicator / weight$count
        U = U + U.decision
        W = W + W.decision
      }#i, individual id
      
      value = U / W
      
      if (v == 1){
        values = value
      } else {
        values = c(values, value)
      }
    }#v, fold
    
    # To maximize f.value is to minimize -f.values
    neg.values = - mean(values) 
    
    return(neg.values)
  }#if max_value
  
  if (type == "min_mse"){
    J = length(pred.files)
    
    loss = 0
    
    # # IDs for each patient (same order)
    # ids = pred.files[[1]][[1]]$id
    
    ## For each super learner repetition
    for (v in 1:k){
      # calculate the weighted average of for all methods
      O = pred.files[[1]][[v]]$pred_outcome * alpha[1]     
      for (j in 2:J){
        O = O + pred.files[[j]][[v]]$pred_outcome * alpha[j]
      }
      
      # Calculate MSE loss
      mse = mean((pred.files[[1]][[1]]$obs_outcome - O)^2)
      
      if (v == 1){
        loss = mse
      } else {
        loss = c(loss, mse)
      }
    }#v, fold
    
    return(mean(loss)) # take mean across all CV repetitions
  }#if min_mse
}


# super learning function to get preds using SL (min_mse), to be fed into modified_list_dtr
Super_Learning = function(type, dat, num_reps, num_folds, num_trt, remove_ids, seed = 2019, n.try, n.pick, plot = FALSE){
  
  #' @param type if "max_value", objective function is value function; if "min_mse", objective function is mse loss of prediction
  #' @param num_reps number of repetition of cross validations
  #' @param remove_ids if not using the whole input data (i.e., one id if jk, a list of ids if cv)
  #' @param dat original dataset
  #' @param n.try randomly generate n.try alpha
  #' @param n.pick pick the best n.pick alpha
  #' @param plot whether to plot the f.record that shows how alphas are updated to improve objective function
  
  library(caret)
  
  # Base learners (can modify this to make it user specified)
  base.learners <- c("bart", "elastic_net", "krr", "lasso", "rf", "ridge", "rlt")
   
  ## Calculate cross validated value functions from all base learners ##
  dat %<>% mutate(pid = 1:nrow(dat))
  cvpred.files <- list()
  
  # Loop through each model to get CV predictions
  for (j in 1:length(base.learners)){ 
    base.learner <- base.learners[j]
    cat("\nBase Learner:", base.learner, "\n")
    
    # Prediction data
    pred_result.all <- NULL
    for (m in 1:num_reps){
      cat("|| Repetition ", m, "||\n")
      set.seed(seed.fn(2,m))
      folds <- createFolds(as.factor(dat$a), k = num_folds) # stratified CV
      
      # loop through each of the CV folds
      for (i in 1:num_folds){ 
        pred_result <- NULL
        
        dat.train = dat[-folds[[i]],] # K-1 folds 
        dat.test = dat[folds[[i]],] # 1 fold
        x.train <- dat.train %>% dplyr::select(contains("x"), a)  
        x.test <- dat.test %>% dplyr::select(contains("x"), a) 
        
        ## Fit training and testing model
        outcome = "y"
        if (base.learner == "rf"){
          pred_result = RF(outcome, x.train, x.test, dat.train, dat.test)
        } else if (base.learner == "rlt"){
          pred_result = rlt(outcome, x.train, x.test, dat.train, dat.test)
        } else if (base.learner == "krr"){
          pred_result = KRR(outcome, x.train, x.test, dat.train, dat.test)
        } else if (base.learner == "lasso"){
          pred_result = PENAL(outcome, x.train, x.test, dat.train, dat.test, 1)
        } else if (base.learner == "ridge"){
          pred_result = PENAL(outcome, x.train, x.test, dat.train, dat.test, 0)
        } else if (base.learner == "elastic_net"){
          pred_result = PENAL(outcome, x.train, x.test, dat.train, dat.test, 0.5)
        } else if (base.learner == "bart"){
          pred_result = BART(outcome, x.train, x.test, dat.train, dat.test)
        } else{
          print("ERROR: Invalid model name!")
        }
        
        pred_result.all <- rbind(pred_result.all, 
                                 pred_result %>% 
                                   dplyr::select(pred1, pred2, pred3, rand, id) %>% 
                                   mutate(rep = m, fold = i))
        
      }#i loop for num_folds
      
    }#m, reptitions of CVs
    
    # Organize optTx_result.all for each outcome
    colnames(pred_result.all) <- c("pred1", "pred2", "pred3", "rand", "id", "rep", "fold")
    cvpreddata <- pred_result.all %>% select(id, rep, fold, rand, everything()) %>% arrange(id, rep, fold) 
    
    cvpredoptdata <- cvpreddata %>% 
      left_join(dat %>% select(pid, y, a), by = c("id" = "pid", "rand" = "a")) %>% 
      rename(obs_outcome = y) %>% 
      mutate(pred_outcome = ifelse(rand == 1, pred1, ifelse(rand == 2, pred2, pred3))) %>% 
      select(-pred1, -pred2, -pred3, -contains("opt_group"))
    
    # remove the test set person, just training set people (unless remove_ids = NA)
    cvpredoptdata %<>% filter(!(id %in% remove_ids)) 
    
    # split prediction data into lists by rep 
    # (if by unique fold, use preddata, fold, and k = num_reps * num_folds)
    if (type == "max_value"){
      cvpred.files[[j]] <- split(cvpreddata[which(names(cvpreddata) %in% c(paste0("pred", 1:num_trts), "id", "rand"))], f = cvpreddata$rep)
    } 
    if (type == "min_mse"){
      cvpred.files[[j]] <- split(cvpredoptdata, f = cvpredoptdata$rep) 
    }
    
  }#j, end of loop for base.learners

  # Training IDs
  ids <- setdiff(unique(dat$pid), remove_ids)
  J <- length(base.learners)
  stopifnot(nrow(dat) == length(unique(cvpreddata$id))) # number of subjects
  stopifnot(num_folds == length(unique(cvpreddata$fold))) # number of cv folds in the above files 
  stopifnot(num_reps == length(unique(cvpreddata$rep))) # number of repetitions
  
  # Get propensity scores
  prop_scores = dat %>% group_by(a) %>% summarise(count = n() / nrow(dat)) 
  
  ## Simulated annealing to find the weight optimization ##
  # Yunshu Zhang at NCSU helped providing code for simulated annealing 
  # Inspired by manuscript "High dimensional precision medicine from patient-derived xenografts" (Rashid et al. 2019 submitted to JASA)
  cat("\nSimulated annealing")
  # randomly generate n.try alpha
  # and pick the best n.pick alpha
  # n.try = 1000, 2000, 3000
  # n.pick = 10, 20, 30
  alpha.random = matrix(data = 0, nrow = n.try, ncol = J)
  f.values = rep(0, n.try)
  
  set.seed(seed)
  for(i in 1:n.try){
    alpha.random[i,] = get.simplex(J)
    f.values[i] = obj.function(type, alpha.random[i,], outcome, cvpred.files, num_reps, prop_scores)
  }
  
  # pick the best n.pick alpha and randomly choose from the tie case
  alpha.start = alpha.random[which(f.values %in% sort(f.values)[1:n.pick]),][1:n.pick,]
  
  n.stage = 60
  temp = rep(5, n.stage)
  for(i in 2:n.stage){
    temp[i] = 0.9 * temp[i - 1] # decaying temperature
  }
  stage.length = rep(c(60,120,220), each = 20)
  e = 0.5
  
  print(paste0("n.stage: ", n.stage, ", ", "total stage length: ", sum(stage.length)))
  
  # record the final alpha and function value
  alpha.final = alpha.start
  f.final = rep(0, n.pick)
  
  f.record = rep(0, sum(stage.length)) # for plotting
  
  for(start in 1:n.pick){
    # cat("n.pick = ", start, "\n")
    alpha.old = alpha.start[start,]
    f.old = obj.function(type, alpha.old, outcome, cvpred.files, num_reps, prop_scores)
    
    step = 1
    
    for(stage in 1:n.stage){
      # cat("stage =", stage, "\n")
      for(i in 1:stage.length[stage]){
        # cat("i =", i, "\n")
        alpha.tmp = get.simplex(J)
        alpha.new = alpha.old * (1-e) + alpha.tmp * e
        f.new = obj.function(type, alpha.new, outcome, cvpred.files, num_reps, prop_scores)
        
        if(runif(1) < exp((f.old - f.new)/temp[stage])){
          # cat(f.new, f.old, runif(1) < exp((f.old - f.new)/temp[stage]), "\n")
          f.old = f.new
          alpha.old = alpha.new
        }#acceptance condition
        
        f.record[step] = f.old
        step = step + 1
      }#i, within each stage
    }#stage, go through all stages
    
    alpha.final[start,] = alpha.old
    f.final[start] = f.old
  }#start, go through n.picks
  
  # Find the best f and alpha
  f.best = min(f.final)
  
  if (length(which(f.final == f.best)) == 1){ # if no ties in the smallest f.final
    alpha.best = alpha.final[which(f.final == f.best),]
  } else{ # if ties in the smallest f.final, select the very first alpha.final
    alpha.best = alpha.final[which(f.final == f.best),][1,]
  }  
  
  ## Reevaluate with alpha.best ## 
  print("Reevaluate with alpha.best")
  cat("f.best is: ", f.best, "\n")
  cat("alpha.best is: ", alpha.best, "\n############################################################\n")
  
  ## Loop through each treatment to get SL predictions ##
  sl.preds <- matrix(NA, nrow = length(ids), ncol = num_trt) %>% as.data.frame()
  colnames(sl.preds) <- 1:num_trt
  
  # Calculate the weighted average of for all methods
  # Loop through each model to get full predictions (no second layer of lists because no repetitions)
  pred.files <- list()
  
  for (j in 1:length(base.learners)){ 
    # Calculate prediction data (full data, no CV)
    base.learner <- base.learners[j]
    pred_result <- NULL
    
    dat <- dat %>% filter(!(pid %in% remove_ids))
    dat.train <- dat
    dat.test <- dat
    x.train <- dat.train %>% select(contains("x"), a)
    x.test <- dat.test %>% select(contains("x"), a)

    ## Start fitting training and testing models
    if (base.learner == "rf"){
      pred_result = RF(outcome, x.train, x.test, dat.train, dat.test)
    } else if (base.learner == "rlt"){
      pred_result = rlt(outcome, x.train, x.test, dat.train, dat.test)
    } else if (base.learner == "krr"){
      pred_result = KRR(outcome, x.train, x.test, dat.train, dat.test)
    } else if (base.learner == "lasso"){
      pred_result = PENAL(outcome, x.train, x.test, dat.train, dat.test, 1)
    } else if (base.learner == "ridge"){
      pred_result = PENAL(outcome, x.train, x.test, dat.train, dat.test, 0)
    } else if (base.learner == "elastic_net"){
      pred_result = PENAL(outcome, x.train, x.test, dat.train, dat.test, 0.5)
    } else if (base.learner == "bart"){ 
      pred_result = BART(outcome, x.train, x.test, dat.train, dat.test)
    } else{
      print("ERROR: Invalid model name!")
    }
    
    pred.files[[j]] <- pred_result %>% arrange(id)
  }#j, end of loop for models
  
  for (t in 1:num_trt){
    O = pred.files[[1]][[paste0("pred", t)]] * alpha.best[1]     
    for (j in 2:J){
      O = O + pred.files[[j]][[paste0("pred", t)]] * alpha.best[j]
    }#j, base learner
    
    sl.preds[,t] <- O
  }#t, num_trts
  
  # write_csv(sl.preds, paste0("./superlearner/C1_listsl_pred_outcomes_", datacode, "_", outcome, ".csv"))
  return(sl.preds)
}

