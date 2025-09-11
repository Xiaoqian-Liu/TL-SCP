# TL-SCP
A transfer-learning based statistical framework for survival-based clustering of predictors, with its application to germline TP53 mutation annotation.

## Introduction
This repository includes all codes and simulation scripts for reproducing the results in the manuscript titled "Transfer Learning for Survival-based Clustering of Predictors with an Application to TP53 Mutation Annotation". The `Code` folder includes all functions needed to implement the basic SCP method and its transfer learning extension, TL-SCP. The `Simulations` folder includes all scripts for reproducing the simulation results in the paper.  The `TP53 Annotation` folder includes all scripts for reproducing the results for our study of TP53 mutation annotation. 

## Prerequisites

To run SCP/TL-SCP with the penalty Lasso, SCAD, or MCP, please install the Gurobi R interface in your RStudio. Details of the installation are referred to [here](https://support.gurobi.com/hc/en-us/articles/14462206790033-How-do-I-install-Gurobi-for-R).

## Implementation

To run SCP, use the `SCP` function with the required input.

To run TL-SCP, use the `trans_SCP` function withe the required input.

We direct users to the function files under the `Code` folder for detailed instruction on using the two main functions.



## Contact

Please reach out to xiaoqian.liu@ucr.edu for any questions. 



