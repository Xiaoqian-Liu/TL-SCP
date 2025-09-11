source('~/Code/cluster_label.R')
source('~/Code/fused_Dist.R')
source('~/Code/fused_Lasso.R')
source('~/Code/fused_MCP.R')
source('~/Code/fused_SCAD.R')
source('~/Code/myGetDkn.R')
source('~/Code/pcox_path.R')
source('~/Code/pcox.R')
source('~/Code/pll_cox.R')
source('~/Code/wz.R')
source('~/Code/wRankAveg.R')
source('~/Code/trans_SCP.R')
source('~/Code/SCP.R')


library(glmnet)
library(Matrix)
library(survival)
library(gurobi)
library(gtools)
library(dplyr)
library(ggplot2)
library(factoextra)

# set dimensions
n <- 300
p <- 30
sig <- 0.5

cs_rate <- c(0.1, 0.3, 0.5, 0.7)
nRate <- length(cs_rate)

nReps <- 50

# parameters
lambdaSeq <- 10^seq(0, -4, length.out=30)
lambdaSeq2 <- 10^seq(0, -4, length.out=30)
kSeq <- 1:5
rho <- 1
maxiters <- 1e3
tol <- 1e-5

# define outputs
Err_Lasso <- matrix(NA, nrow = nRate, ncol = nReps)
Err_SCAD <- matrix(NA, nrow = nRate, ncol = nReps)
Err_MCP <- matrix(NA, nrow = nRate, ncol = nReps)
Err_Dist <- matrix(NA, nrow = nRate, ncol = nReps)
Err_cox <- matrix(NA, nrow = nRate, ncol = nReps)


# AMI
AMI_Lasso <- matrix(NA, nrow = nRate, ncol = nReps)
AMI_SCAD <- matrix(NA, nrow = nRate, ncol = nReps)
AMI_MCP <- matrix(NA, nrow = nRate, ncol = nReps)
AMI_Dist <- matrix(NA, nrow = nRate, ncol = nReps)
AMI_cox <- matrix(NA, nrow = nRate, ncol = nReps)


# max variant of NMI
NMI_M_Lasso <- matrix(NA, nrow = nRate, ncol = nReps)
NMI_M_SCAD <- matrix(NA, nrow = nRate, ncol = nReps)
NMI_M_MCP <- matrix(NA, nrow = nRate, ncol = nReps)
NMI_M_Dist <- matrix(NA, nrow = nRate, ncol = nReps)
NMI_M_cox <- matrix(NA, nrow = nRate, ncol = nReps)

# record beta paramter estimates
Beta_Lasso <- list()
Beta_SCAD <- list()
Beta_MCP <- list()
Beta_cox <- list()
Beta_Dist <- list()

k_Dist <- matrix(NA, nrow = nRate, ncol = nReps)

# record running time 
time_Lasso <- matrix(NA, nrow = nRate, ncol = nReps)
time_SCAD <- matrix(NA, nrow = nRate, ncol = nReps)
time_MCP <- matrix(NA, nrow = nRate, ncol = nReps)
time_Dist <- matrix(NA, nrow = nRate, ncol = nReps)
time_cox <- matrix(NA, nrow = nRate, ncol = nReps)



for (i in 1:nRate) {
  
  set.seed(2025)
  # fix p and beta
  beta <- c(rep(2, p/3), rep(0.5, p/3), rep(-1, p/3))
  beta <- permute(beta)
  cluster.true <- cluster_label(beta)
  # fix X
  X <- matrix(runif(n*p, -1/2, 1/2), nrow = n, ncol = p)
  
  # fix censor_rate
  censor_rate <- cs_rate[i]
  
  beta_Lasso <- matrix(NA, nrow = p, ncol = nReps)
  beta_SCAD <- matrix(NA, nrow = p, ncol = nReps)
  beta_MCP <- matrix(NA, nrow = p, ncol = nReps)
  beta_cox <- matrix(NA, nrow = p, ncol = nReps)  
  beta_Dist <- matrix(NA, nrow = p, ncol = nReps)
  

  for (j in 1:nReps) {
    
    cat(paste('Current stage:', i, 'case', j, 'reps '))
    
    set.seed(2025+j)
    eta <- as.vector(X%*%beta+sig*rnorm(n)) 
    y <- rep(NA, n)
    for (k in 1:n) {
      y[k] <- rexp(1, rate=exp(eta[k]))
    }
    # get the status 
    status <- rep(1, n)
    status[sample(1:n, size=censor_rate*n)] <- 0
    
    Y <- cbind(time=y, status=status)
    
    ##########################################
    ##########################################
    
    ### Lasso
    t0 <- proc.time()
    fit_bcards_lasso <- SCP(y, status, X, penalty="Lasso", lambdaSeq = lambdaSeq,
                                   ebic.gm=0, maxiters=maxiters, tol=tol)
    t1 <-proc.time()-t0
    time_Lasso[i, j] <- t1[[3]]
    idx1 <- which.min(fit_bcards_lasso$IC_path)
    beta_Lasso[, j] <- fit_bcards_lasso$Beta_path[, idx1]
    Err_Lasso[i, j] <- norm(beta_Lasso[, j] - beta, "2")^2/norm(beta, "2")^2
    k_select <- fviz_nbclust(as.matrix(beta_Lasso[, j]), kmeans, method = "silhouette")
    cluster.Lasso <- cluster_label_kmeans(beta_Lasso[, j], which.max(k_select$data$y))
    AMI_Lasso[i, j] <- aricode::AMI(cluster.true, cluster.Lasso)
    NMI_M_Lasso[i, j] <- aricode::NMI(cluster.true, cluster.Lasso)
    cat('Lasso Done ')
    
    
    

    
    ### SCAD
    t0 <- proc.time()
    fit_bcards_scad <- SCP(y, status, X, penalty="SCAD", lambdaSeq = lambdaSeq,
                                  ebic.gm=0, maxiters=maxiters, tol=tol)
    t3 <-proc.time()-t0
    time_SCAD[i, j] <- t3[[3]]
    idx3 <- which.min(fit_bcards_scad$IC_path)
    beta_SCAD[, j] <- fit_bcards_scad$Beta_path[, idx3]
    Err_SCAD[i, j] <- norm(beta_SCAD[, j] - beta, "2")^2/norm(beta, "2")^2
    k_select <- fviz_nbclust(as.matrix(beta_SCAD[, j]), kmeans, method = "silhouette")
    cluster.SCAD <- cluster_label_kmeans(beta_SCAD[, j], which.max(k_select$data$y))
    AMI_SCAD[i, j] <- aricode::AMI(cluster.true, cluster.SCAD)
    NMI_M_SCAD[i, j] <- aricode::NMI(cluster.true, cluster.SCAD)
    
    cat('SCAD Done ')
    
    
    ### MCP
    t0 <- proc.time()
    fit_bcards_mcp <- SCP(y, status, X, penalty="MCP", lambdaSeq = lambdaSeq,
                                 ebic.gm=0, maxiters=maxiters, tol=tol)
    t4 <-proc.time()-t0
    time_MCP[i, j] <- t4[[3]]
    idx4 <- which.min(fit_bcards_mcp$IC_path)
    beta_MCP[, j] <- fit_bcards_mcp$Beta_path[, idx4]
    Err_MCP[i, j] <- norm(beta_MCP[, j] - beta, "2")^2/norm(beta, "2")^2
    k_select <- fviz_nbclust(as.matrix(beta_MCP[, j]), kmeans, method = "silhouette")
    cluster.MCP <- cluster_label_kmeans(beta_MCP[, j], which.max(k_select$data$y))
    AMI_MCP[i, j] <- aricode::AMI(cluster.true, cluster.MCP)
    NMI_M_MCP[i, j] <- aricode::NMI(cluster.true, cluster.MCP)
    
    cat('MCP Done ')
    
    
    
    ### Cox together
    t0 <- proc.time()
    fit_cox <- glmnet(X, Y, family = "cox", lambda = 0)
    beta0 <- as.vector(fit_cox$beta) 
    
    # check for zero beta estimates
    if(norm(beta0, "2")==0){
      fit_cox <- glmnet(X, Y, family = "cox", lambda = 1e-2)
      beta0 <- as.vector(fit_cox$beta)
    }
    t5 <-proc.time()-t0
    time_cox[i, j] <- t5[[3]]
    beta_cox[, j] <- beta0
    Err_cox[i, j] <- norm(beta_cox[, j] - beta, "2")^2/norm(beta, "2")^2
    k_select <- fviz_nbclust(as.matrix(beta_cox[, j]), kmeans, method = "silhouette")
    cluster.cox <- cluster_label_kmeans(beta0, which.max(k_select$data$y))
    AMI_cox[i, j] <- aricode::AMI(cluster.true, cluster.cox)
    NMI_M_cox[i, j] <- aricode::NMI(cluster.true, cluster.cox)
    
    cat('Cox Done ')
    
    
    ### Distance
    t0 <- proc.time()
    fit_bcards_dist <- SCP(y, status, X, penalty="distance", nclusters=kSeq,
                                  ebic.gm=0, maxiters=maxiters, tol=tol)
    
    t6 <-proc.time()-t0
    time_Dist[i, j] <- t6[[3]]
    idx6 <- which.min(fit_bcards_dist$IC_path)
    beta_Dist[, j] <- fit_bcards_dist$Beta_path[, idx6]
    Err_Dist[i, j] <- norm(beta_Dist[, j] - beta, "2")^2/norm(beta, "2")^2 
    k_Dist[i, j] <- kSeq[idx6]
    cluster.Dist <- cluster_label_dist(beta_Dist[, j], k_Dist[i, j])
    AMI_Dist[i, j] <- aricode::AMI(cluster.true, cluster.Dist)
    NMI_M_Dist[i, j] <- aricode::NMI(cluster.true, cluster.Dist)
    
    cat('Dist Done \n')
    
    save.image(file = "SCP_censor.RData")
  }
  
  cat('case', i, 'finish \n')
  
  Beta_Lasso[[i]] <- beta_Lasso
  Beta_SCAD[[i]] <- beta_SCAD
  Beta_MCP[[i]] <- beta_MCP
  Beta_cox[[i]] <- beta_cox
  Beta_Dist[[i]] <- beta_Dist
  
  save.image(file = "SCP_censor.RData")
}
 