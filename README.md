# TL-SCP
A transfer-learning based statistical framework for survival-based clustering of predictors, with its application to germline *TP53* mutation annotation.

*Jointly developed by [Wenyi Wang](https://odin.mdacc.tmc.edu/~wwang7/index.html)'s lab and [Xiaoqian Liu](https://xiaoqian-liu.github.io/), 
actively maintained by [Xiaoqian Liu](https://xiaoqian-liu.github.io/).*

## Introduction
This repository includes all codes and simulation scripts for reproducing the results in the manuscript titled "Transfer Learning for Survival-based Clustering of Predictors with an Application to *TP53* Mutation Annotation". The `Code` folder includes all functions needed to implement the basic SCP method and its transfer learning extension, TL-SCP. The `Simulations` folder includes all scripts for reproducing the simulation results in the paper.  The `TP53 Annotation` folder includes all scripts for reproducing the results for our study of *TP53* mutation annotation. 

## Prerequisites

To run SCP/TL-SCP with the Lasso, SCAD, or MCP penalty, you need to use the [Gurobi R interface](https://docs.gurobi.com/projects/optimizer/en/current/reference/r.html). 
To do that, you need to first install the Gurobi optimizer in your personal machine. Academics can obtain a free academic named-user license by following this [instruction page](https://www.gurobi.com/features/academic-named-user-license/) . 
After that, you can install the gurobi R package by referring to this [installation guidance](https://docs.gurobi.com/projects/optimizer/en/current/reference/r/setup.html).

## Implementation

To run SCP, use the `SCP` function with the required input.

To run TL-SCP, use the `trans_SCP` function withe the required input.

We direct users to the function files under the `Code` folder for detailed instruction on using the two main functions.




## Contact

Please reach out to xiaoqian.liu@ucr.edu for any questions. 

## Citation

Xiaoqian Liu, Hao Yan, Haoming Shi, Emilie Montellier, Eric C. Chi, Pierre Hainaut, Wenyi Wang (2025). Transfer Learning for Survival-based Clustering of Predictors with an Application to TP53 Mutation Annotation. 
bioRxiv 2025.10.06.680732; doi: https://doi.org/10.1101/2025.10.06.680732

