#####
# PROJECT. Precision Medicine for Knee OA
# AUTHOR. Xiaotong Jiang
#####
# USAGE. Compare the optimal PMM and optimal ZOM. Run main.R for PMOA.
#####

F4_PMM_vs_ZOM <- function(estimator, mean_var.pmm, mean_var.zom, UW.pmm, UW.zom, data.type){
  
  #' Find and compare the optimal PMM and optimal ZOM with a Z-test
  #' 
  #' @param estimator string, which estimator to use for the value function, originally {"empirical", "jackknife", "empirical+test", "jackknife+test"}; for the purpose of this example code, estimator is set to be "jackknife"
  #' @param mean_var.pmm data.frame, mean and variance of PMMs 
  #' @param mean_var.zom data.frame, mean and variance of ZOMs
  #' @param UW.pmm data.frame, long version of PMM results UW
  #' @param UW.zom data.frame, long version of ZOM results UW
  #' @param data.type string, type of data, {"circle", "steps", "line", "quadratic", "null"}; "null" means one treatment is better than the other for all X values 
  #' 
  #' @return 
  
  if (estimator == "jackknife"){
    n = unique(mean_var.pmm$n) # sample size
    num_folds = n # jackknife is n-fold CV
    num_reps = 1 # rep if CV
  }
 
  optPMM <- mean_var.pmm %>% filter(V == max(V))
  optZOM <- mean_var.zom %>% filter(V == max(V))
  
  UW.PMM <- UW.pmm %>% filter(model == optPMM$model) %>% arrange(sim, test.n, model, which.Z)
  UW.ZOM <- UW.zom %>% filter(model == optZOM$model) %>% arrange(sim, test.n, model, which.Z)
  
  UW.PMM.list <- list(Umk = UW.PMM %>% filter(which.Z == "Ui") %>% select(value) %>% unlist() %>% unname(),
                      Wmk = UW.PMM %>% filter(which.Z == "Wi") %>% select(value) %>% unlist() %>% unname())
  UW.ZOM.list <- list(Umk = UW.ZOM %>% filter(which.Z == "Ui") %>% select(value) %>% unlist() %>% unname(),
                      Wmk = UW.ZOM %>% filter(which.Z == "Wi") %>% select(value) %>% unlist() %>% unname())
  
  Vhat.PMM <- Get_Value_CV(UW.PMM.list, num_folds)[["mean_value"]]
  Vhat.ZOM <- Get_Value_CV(UW.ZOM.list, num_folds)[["mean_value"]]

  Ri.PMM <- Get_Rmk_CV(UW.PMM.list)
  Ri.ZOM <- Get_Rmk_CV(UW.ZOM.list)
  
  T.statistic <- ifelse(Vhat.PMM == Vhat.ZOM, 0,
                        (Vhat.PMM - Vhat.ZOM) / sqrt( sum((Ri.PMM - Ri.ZOM)^2) / (num_folds * (num_reps*num_folds-1))))
  
  p.value <- 2 *(1-pnorm(abs(T.statistic))) # one-side, T statistic cannot be negative, but multiply by 2 because the distribution under the null is the absolute value of a standard normal.
  
  test.result <- data.frame(n = n, 
                            data.type = data.type,
                            optimalPMM = optPMM$model, 
                            optimalZOM = optZOM$model, 
                            vhatPMM = Vhat.PMM, 
                            vhatZOM = Vhat.ZOM, 
                            test.stat = T.statistic, 
                            p.value = p.value) 
  
  return(test.result)
}

