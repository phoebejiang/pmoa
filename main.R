#####
# PROJECT. Precision Medicine for Knee OA
# AUTHOR. Xiaotong Jiang
#####
# USAGE. Execute all precision medicine analysis (using simulation data). Run this script main.R for PMOA; the rest of the scripts are supporting functions.
#####
# It is highly recommended to run this script on cluster. Step 2 will take the longest time.
#####


wd <- "/Users/Jiang/Documents/GitHub/PMOA/scripts/simulation_public/"
setwd(wd)
source("./F_utils.R")
source("./F_superlearner.R")
source("./F1_gensimdata.R")
source("./F2_PMM.R")
source("./F3_ZOM.R")
source("./F4_PMM_vs_ZOM.R")

## User-specified constants
n = 200
data.type = "circle"
num_rep = 10
num_folds = 5
n.try = 1000
n.pick = 10
estimator = "jackknife"
size.rule <- c(2L, 3L, 5L, 10L)       # number of list nodes 
models <- c("krr", 
            paste0("list.",size.rule), 
            paste0("list", size.rule, "+rf"), 
            "list10+en", 
            "lasso", "ridge", "elastic_net", 
            "rf", "rlt", 
            "rwl-linear", "rwl-radial",
            paste0("list", size.rule, "+sl"), "bart")  # models listed in Table 2
models <- c("krr", paste0("list.",size.rule), 
            paste0("list", size.rule, "+rf"), "list10+en",  
            "lasso", "ridge", "elastic_net", "rf", "rwl-linear")

## Step 1. Generate simulation data
f1.list <- F1_genSimData(n.sim = 1, n = n, data.type = data.type)
dat <- f1.list[[1]]
test.dat <- f1.list[[2]]

## Step 2. Fit each precision medicine model (PMM) with jackknife estimator
f2.list <- F2_PMM(n.sim = 1, n = n, data.type = data.type, models = models,
                       dat = dat, test.dat = test.dat, estimator = estimator, 
                       num_rep = num_rep, num_folds = num_folds, n.try = n.try, n.pick = n.pick)
UW.pmm <- f2.list[[1]]
mean_var.pmm <- f2.list[[2]]

## Step 3. Fit each ZOM 
f3.list <- F3_ZOM(n.sim = 1, n = n, data.type = data.type, dat = dat, test.dat = test.dat, estimator = estimator)
UW.zom <- f3.list[[1]]
mean_var.zom <- f3.list[[2]]

## Step 4. Compare the optimal PMM and ZOM with Z-test
test.result <- F4_PMM_vs_ZOM(estimator, mean_var.pmm, mean_var.zom, UW.pmm, UW.zom, data.type)

## Step 5. Determine optimal decision rule 
optPMM <- test.result$optimalPMM

model <- ifelse(grepl("^list[0-9]+\\+.*", optPMM), 
                sub("[0-9]+\\+.*", "", optPMM), # if list+rf
                sub("[\\.|-].*$", "", optPMM)) # if not list+rf

node <- ifelse(grepl("^list[0-9]+\\+.*$", optPMM), 
               optPMM %>% str_match_all("[0-9]+") %>% unlist(),  # if list+rf
               str_split(optPMM, "\\.")[[1]][2]) # if not list+rf

Q.model <- ifelse(grepl("^list[0-9]+\\+", optPMM), 
                  sub("^list[0-9]+\\+", "", optPMM),  NA)

sl.bundle <- NULL # only used for super learning model
if ((Q.model == "sl") & !is.na(Q.model)){sl.bundle = list(dat, num_rep, num_folds, n.try, n.pick, NA)} # remove_ids initialized to be NA for right now

which_wl <- str_split(optPMM, "-")[[1]][2]

train.mod <- Training(data = dat, 
         model = optPMM, 
         size = ifelse(!is.na(node), as.integer(node), NA), 
         kernel = ifelse(!is.na(which_wl), which_wl, NA), 
         Q.model = Q.model, sl.bundle = sl.bundle)
  