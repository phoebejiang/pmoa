# pmoa

## Instructions

This repository provides example code for readers to run and test the pipeline developed in Jiang et al. 2020. 

Run `main.R` and pass the argument `models` with the following values such as `elastic_net, lasso, ridge, krr, list.5, list5+rf, list5+en, list5+sl, rf, rlt, bart, rwl-poly2, rwl-poly3, rwl-radial`. The number of node in list-based models, 5, could be replaced by any reasonable integer (usually between 2 and 10). 

High-performance computing is strongly recommended if candidate models include Bayesian additive regression trees (BART), reinforcement learning tree (RLT), or super learning (SL). 

## Related Publications

1. The original paper.

Jiang et al. "A precision medicine approach to develop and validate optimal exericse and weight loss treatments for overweight and obese adults with knee osteoarthritis (Data from the Intensive Diet and Exercise for Arthritis (IDEA) trial)." Arthritis Care and Research (2020).

2. Statistical report.

Jiang et al. "Technical Background for 'A Precision Medicine Approach to Develop and Internally Validate Optimal Exercise and Weight Loss Treatments for Overweight and Obese Adults with Knee Osteoarthritis'". 
[arXiv:2001.09930](https://arxiv.org/abs/2001.09930) (2020). 


