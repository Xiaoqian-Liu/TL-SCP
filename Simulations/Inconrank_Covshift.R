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
library(latex2exp)
library(factoextra)


#set dimensions
n_T <- 200
n_S <- 2000
n <- n_T+n_S
p <- 90
sSeq <- c(5, 10, 15, 20) # We will permute the first s elements of beta to have a different order
nS <- length(sSeq)

sig <- 0.5
censor_rate <- 0.3

nReps <- 50

# parameters
kSeq <- 1:5
maxiters <- 1e3
tol <- 1e-4



# define outputs
Err_SCP <- matrix(NA, nrow = nS, ncol = nReps)
Err_tSCP <- matrix(NA, nrow = nS, ncol = nReps)
Err_cox <- matrix(NA, nrow = nS, ncol = nReps)

# AMI
AMI_SCP <- matrix(NA, nrow = nS, ncol = nReps)
AMI_tSCP <- matrix(NA, nrow = nS, ncol = nReps)
AMI_cox <- matrix(NA, nrow = nS, ncol = nReps)

# max variant of NMI
NMI_M_SCP <- matrix(NA, nrow = nS, ncol = nReps)
NMI_M_tSCP <- matrix(NA, nrow = nS, ncol = nReps)
NMI_M_cox <- matrix(NA, nrow = nS, ncol = nReps)


# sum variant of NMI
NMI_S_SCP <- matrix(NA, nrow = nS, ncol = nReps)
NMI_S_tSCP <- matrix(NA, nrow = nS, ncol = nReps)
NMI_S_cox <- matrix(NA, nrow = nS, ncol = nReps)

# ARI
ARI_SCP <- matrix(NA, nrow = nS, ncol = nReps)
ARI_tSCP <- matrix(NA, nrow = nS, ncol = nReps)
ARI_cox <- matrix(NA, nrow = nS, ncol = nReps)



# record beta paramter estimates
Beta_SCP <- list()
Beta_tSCP <- list()
Beta_cox <- list()

# record the selected k
k_SCP <- matrix(NA, nrow = nS, ncol = nReps)
k_tSCP <- matrix(NA, nrow = nS, ncol = nReps)
k_cox <- matrix(NA, nrow = nS, ncol = nReps)

# record running time 
time_SCP <- matrix(NA, nrow = nS, ncol = nReps)
time_tSCP <- matrix(NA, nrow = nS, ncol = nReps)
time_cox <- matrix(NA, nrow = nS, ncol = nReps)


for (i in 1:nS) {
  
  set.seed(2024)

  #########################
  ##### Target data #######
  #########################
  beta_T <- c(rep(2, p/3), rep(0.5, p/3), rep(-1, p/3))
  beta_T <- permute(beta_T)
  cluster.true <- cluster_label(beta_T)
  X_T <- matrix(runif(n_T*p, -1/2, 1/2), nrow = n_T, ncol = p)
  
  #########################
  ##### Source data #######
  #########################
  # We assume beta is different, but has partially different ranking
  rk <- rank(beta_T, ties.method = "first")
  #change the ranking of the first s elements
  s <- sSeq[i]
  rk[1:s] <- permute(rk[1:s])
  beta_S <- log(rk)
  X_S <- matrix(runif(n_S*p, -1/2, 1/2), nrow = n_S, ncol = p)+
              matrix(rnorm(n_S*p, 0, 1/4), nrow = n_S, ncol = p)
  
  beta_SCP <- matrix(NA, nrow = p, ncol = nReps)  
  beta_tSCP <- matrix(NA, nrow = p, ncol = nReps)
  beta_cox <- matrix(NA, nrow = p, ncol = nReps) 
  
  for (j in 1:nReps) {
    
    cat(paste('Current stage:', i, 'case', j, 'reps '))
    
    set.seed(2024+j)
    #########################
    ##### Target data #######
    #########################
    eta_T <- as.vector(X_T%*%beta_T+sig*rnorm(n_T)) 
    y_T <- rep(NA, n_T)
    for (k in 1:n_T) {
      y_T[k] <- rexp(1, rate=exp(eta_T[k]))
    }
    # get the status 
    status_T <- rep(1, n_T)
    status_T[sample(1:n_T, size=censor_rate*n_T)] <- 0
    Y_T <- cbind(time=y_T, status=status_T)
    
    #########################
    ##### Source data #######
    #########################   
    eta_S <- as.vector(X_S%*%beta_S+sig*rnorm(n_S)) 
    y_S <- rep(NA, n_S)
    for (k in 1:n_S) {
      y_S[k] <- rexp(1, rate=exp(eta_S[k]))
    }
    # get the status 
    status_S <- rep(1, n_S)
    status_S[sample(1:n_S, size=censor_rate*n_S)] <- 0
    Y_S <- cbind(time=y_S, status=status_S)
    
    
    #########################
    ### SCP + distance
    #########################
    t0 <- proc.time()
    fit_SCP <- SCP(y=Y_T[, 1], status=Y_T[, 2], X_T, penalty="distance", nclusters=kSeq,
                   ebic.gm=0, maxiters=maxiters, tol=tol)
    
    t1 <-proc.time()-t0
    time_SCP[i, j] <- t1[[3]]
    idx1 <- which.min(fit_SCP$IC_path)
    beta_SCP[, j] <- fit_SCP$Beta_path[, idx1]
    Err_SCP[i, j] <- norm(beta_SCP[, j] - beta_T, "2")^2/norm(beta_T, "2")^2 
    k_SCP[i, j] <- kSeq[idx1]
    cluster.SCP <- cluster_label_dist(beta_SCP[, j], k_SCP[i, j])
    AMI_SCP[i, j] <- aricode::AMI(cluster.true, cluster.SCP)
    NMI_M_SCP[i, j] <- aricode::NMI(cluster.true, cluster.SCP)
    NMI_S_SCP[i, j] <- aricode::NMI(cluster.true, cluster.SCP, variant = 'sum')
    ARI_SCP[i, j] <- aricode::ARI(cluster.true, cluster.SCP)
    cat('SCP Done \n')
    
    
    #########################
    ### trans-SCP + distance
    #########################
    t0 <- proc.time()
    fit_tSCP <- trans_SCP(Y_T, X_T, Y_S, X_S, penalty="distance", nclusters=kSeq,
                          ebic.gm=0, maxiters=maxiters, tol=tol)
    
    t2 <-proc.time()-t0
    time_tSCP[i, j] <- t2[[3]]
    idx2 <- which.min(fit_tSCP$IC_path)
    beta_tSCP[, j] <- fit_tSCP$Beta_path[, idx2]
    Err_tSCP[i, j] <- norm(beta_tSCP[, j] - beta_T, "2")^2/norm(beta_T, "2")^2 
    k_tSCP[i, j] <- kSeq[idx2]
    cluster.tSCP <- cluster_label_dist(beta_tSCP[, j], k_tSCP[i, j])
    AMI_tSCP[i, j] <- aricode::AMI(cluster.true, cluster.tSCP)
    NMI_M_tSCP[i, j] <- aricode::NMI(cluster.true, cluster.tSCP)
    NMI_S_tSCP[i, j] <- aricode::NMI(cluster.true, cluster.tSCP, variant = 'sum')
    ARI_tSCP[i, j] <- aricode::ARI(cluster.true, cluster.tSCP)
    cat('tSCP Done \n')
    
    
    
    #########################
    ### cox-kmeans on the target only
    #########################
    t0 <- proc.time()
    fit_cox <- glmnet(X_T, Y_T, family = "cox", lambda = 0)
    beta0 <- as.vector(fit_cox$beta) 
    # check for zero beta estimates
    if(norm(beta0, "2")==0){
      fit_cox <- glmnet(X_T, Y_T, family = "cox", lambda = 1e-2)
      beta0 <- as.vector(fit_cox$beta)
    }
    t3 <-proc.time()-t0
    
    time_cox[i, j] <- t3[[3]]
    beta_cox[, j] <- beta0
    Err_cox[i, j] <- norm(beta_cox[, j] - beta_T, "2")^2/norm(beta_T, "2")^2
    k_select <- fviz_nbclust(as.matrix(beta_cox[, j]), kmeans, method = "silhouette")
    cluster.cox <- cluster_label_kmeans(beta0, which.max(k_select$data$y))
    k_cox[i, j] <- which.max(k_select$data$y)
    AMI_cox[i, j] <- aricode::AMI(cluster.true, cluster.cox)
    NMI_M_cox[i, j] <- aricode::NMI(cluster.true, cluster.cox)
    NMI_S_cox[i, j] <- aricode::NMI(cluster.true, cluster.cox, variant = 'sum')
    ARI_cox[i, j] <- aricode::ARI(cluster.true, cluster.cox)
    
    cat('cox-kmeans Done \n')
    
    save.image(file = "Inconrank_Covshift.RData")
  }
  
  cat('case', i, 'finish \n')
  
  Beta_SCP[[i]] <- beta_SCP
  Beta_tSCP[[i]] <- beta_tSCP
  Beta_cox[[i]] <- beta_cox
  
  save.image(file = "Inconrank_Covshift.RData")
}
