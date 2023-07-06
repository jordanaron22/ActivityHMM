
#### Functions ####
SimulateHMM <- function(day_length,num_ind,covar_num,init,tran,emit_act,emit_light, wake_sleep_corr){
  
  covar_mat <- matrix(rbinom(num_ind*covar_num,1,0.5),num_ind,covar_num)
  covar_mat <- cbind((numeric(num_ind) + 1),covar_mat)
  activity_array <- array(0,dim  = c(day_length,covar_num+1,num_ind))
  
  for (ind in 1:num_ind){
    hidden_states <- numeric(day_length)
    activity_mat_ind <- matrix(NA,day_length,covar_num+1)
    light <- numeric(day_length)
    
    for (i in 1:day_length){
      if (i == 1) {
        hidden_states[1] <- rbinom(1,1,init[2])
      } else {
        hidden_states[i] <- rbinom(1,1,tran[hidden_states[i-1] + 1,2])
      }
      
      
      mu_act <- (covar_mat[ind,] * emit_act[,1,hidden_states[i]+1])
      sig_act <- (covar_mat[ind,] * emit_act[,2,hidden_states[i]+1]) 
      
      mu_light <- emit_light[hidden_states[i] + 1,1]
      sig_light <- emit_light[hidden_states[i] + 1,2]
      
      sigma_mat <- matrix(c(sig_act[1]^2,wake_sleep_corr[hidden_states[i] + 1]*sig_act[1]*sig_light,
                            wake_sleep_corr[hidden_states[i] + 1]*sig_act[1]*sig_light, sig_light^2), nrow = 2, byrow=T)
      act_light <- mvrnorm(1,c(mu_act[1],mu_light),sigma_mat)
      
      activity_mat_ind[i,1] <- act_light[1]
      activity_mat_ind[i,2:3] <- rnorm(covar_num,mu_act[2:3],sig_act[2:3])
      light[i] <- act_light[2]
    }
    
    activity_array[,,ind] <- activity_mat_ind
    activity_vec <- rowSums(activity_mat_ind)
    
    if (ind == 1){
      hidden_states_matrix <- hidden_states
      activity_matrix <- activity_vec
      light_matrix <- light
    } else {
      hidden_states_matrix <- cbind(hidden_states_matrix,hidden_states)
      activity_matrix <- cbind(activity_matrix,activity_vec)
      light_matrix <- cbind(light_matrix,light)
    }
    
    
  }
  
  
  return(list(hidden_states_matrix,activity_matrix,light_matrix,covar_mat,activity_array))
}

logClassification <- function(time,current_state,act,light,emit_act,emit_light,covar_ind,wake_sleep_corr, BiVar = T){
  
  mu_act <- emit_act[1,1,current_state+1]
  covar_act <- (covar_ind[2:3] %*% emit_act[2:3,1,current_state+1])
  sig_act <- sqrt(sum((covar_ind * emit_act[,2,current_state+1])^2))
  
  mu_light <- emit_light[current_state+1,1]
  sig_light <- emit_light[current_state+1,2]
  
  sigma_mat <- matrix(c(sig_act^2,wake_sleep_corr[current_state + 1]*sig_act*sig_light,
                        wake_sleep_corr[current_state + 1]*sig_act*sig_light, sig_light^2), nrow = 2, byrow=T)
  
  ####HERE
  #
  #
  if (!is.na(act[time]) & !is.na(light[time])) {
    if (BiVar){
      lognorm_dens <- log(mvtnorm::dmvnorm(x = c(act[time] - covar_act,light[time]), mean = c(mu_act,mu_light), sigma = sigma_mat, checkSymmetry = F))
    } else {
      lognorm_dens <- log(dnorm(act[time] - covar_act,mu_act,sig_act)) + log(dnorm(light[time],mu_light,sig_light))
    }
    
  } else if (!is.na(act[time]) & is.na(light[time])){
    lognorm_dens <- log(dnorm(act[time] - covar_act,mu_act,sig_act)) 
    
  } else if (is.na(act[time]) & !is.na(light[time])){
    lognorm_dens <- log(dnorm(light[time],mu_light,sig_light))
    
  } else {lognorm_dens <- 0}
    

  if (lognorm_dens == -Inf){
    lognorm_dens <- log(.Machine$double.xmin)
  }
  
  return(lognorm_dens)
}

Backward <- function(act,light,tran,emit_act,emit_light,covar_mat,wake_sleep_corr){
  # beta_list <- list()
  
  beta_list <- foreach(ind = 1:dim(act)[2]) %dopar% {
    act_ind <- act[,ind]
    light_ind <- light[,ind]
    BackwardInd(act_ind,light_ind,tran,emit_act,emit_light,covar_mat[ind,],wake_sleep_corr)
  }
  return(beta_list)
}

BackwardInd <- function(act, light, tran, emit_act, emit_light,covar_ind,wake_sleep_corr) {
  n <- length(act)
  beta <- matrix(0, ncol = 2, nrow = n)
  
  beta[n,1] <- log(1)
  beta[n,2] <- log(1)
  
  
  for (i in (n-1):1){
    #State 0 from 0
    bp_00 <- log(tran[1,1]) + logClassification(i+1,0,act,light,emit_act,emit_light,covar_ind,wake_sleep_corr) + beta[i+1,1]
    
    #State 1 from 0
    bp_01 <- log(tran[1,2]) + logClassification(i+1,1,act,light,emit_act,emit_light,covar_ind,wake_sleep_corr) + beta[i+1,2] 
    
    #State 0 from 1
    bp_10 <- log(tran[2,1]) + logClassification(i+1,0,act,light,emit_act,emit_light,covar_ind,wake_sleep_corr) + beta[i+1,1] 
    
    #State 1 from 1
    bp_11 <- log(tran[2,2]) + logClassification(i+1,1,act,light,emit_act,emit_light,covar_ind,wake_sleep_corr) + beta[i+1,2] 
    
    
    beta[i,1] <- logSumExp(c(bp_00,bp_01))
    beta[i,2] <- logSumExp(c(bp_10,bp_11))
  }
  
  return(beta)
  
}

Forward <- function(act,light,init,tran,emit_act,emit_light,covar_mat,wake_sleep_corr){
  # alpha_list <- list()
  
  alpha_list <- foreach(ind = 1:dim(act)[2]) %dopar% {
    act_ind <- act[,ind]
    light_ind <- light[,ind]
  
    ForwardInd(act_ind,light_ind,init,tran,emit_act,emit_light,covar_mat[ind,],wake_sleep_corr)
  }
  return(alpha_list)
}

ForwardInd <- function(act, light, init, tran, emit_act, emit_light,covar_ind,wake_sleep_corr) {
  alpha <- matrix(0, ncol = 2, nrow=length(act))
  
  alpha[1,1] <- log(init[1]) + logClassification(1,0,act,light,emit_act,emit_light,covar_ind,wake_sleep_corr)
  
  alpha[1,2] <- log(init[2]) + logClassification(1,1,act,light,emit_act,emit_light,covar_ind,wake_sleep_corr)
  
  for (i in 2:length(act)){
    
    #From state 0 to 0
    fp_00 <- alpha[i-1,1] + log(tran[1,1]) + logClassification(i,0,act,light,emit_act,emit_light,covar_ind,wake_sleep_corr)
    
    #From state 1 to 0
    fp_10 <- alpha[i-1,2] + log(tran[2,1]) + logClassification(i,0,act,light,emit_act,emit_light,covar_ind,wake_sleep_corr)
    
    #From state 0 to 1
    fp_01 <- alpha[i-1,1] + log(tran[1,2]) + logClassification(i,1,act,light,emit_act,emit_light,covar_ind,wake_sleep_corr)
    
    #From state 1 to 1
    fp_11 <- alpha[i-1,2] + log(tran[2,2]) + logClassification(i,1,act,light,emit_act,emit_light,covar_ind,wake_sleep_corr)
    
    alpha[i,1] <- logSumExp(c(fp_00,fp_10))
    alpha[i,2] <- logSumExp(c(fp_01,fp_11))
  }
  
  return(alpha)
}

CalcInitInd <- function(alpha,beta){
  return(ProbMC(alpha,beta,1))
}

CalcInit <- function(alpha,beta){
  
  init_0_vec <- c()
  init_1_vec <- c()
  
  for(ind in 1:length(alpha)){
    
    alpha_ind <- alpha[[ind]]
    beta_ind <- beta[[ind]]
    likelihood <- logSumExp(c(alpha_ind[dim(alpha_ind)[1],]))
    
    init_0_vec <- c(init_0_vec,alpha_ind[1,1] + beta_ind[1,1] - likelihood)
    init_1_vec <- c(init_1_vec,alpha_ind[1,2] + beta_ind[1,2] - likelihood)
    
  }
  
  init <- c(logSumExp(init_0_vec),logSumExp(init_1_vec))
  
  return(exp(init - logSumExp(init)))
}

ProbMC <- function(alpha,beta,time){
  denom <- logSumExp(c(alpha[time,] + beta[time,]))
  num <- alpha[time,] + beta[time,]
  return(exp(num-denom))
}

ProbWeights <- function(alpha,beta){
  weights_mat <- matrix(0,dim(alpha[[1]])[1],length(alpha))
  for (i in 1:length(alpha)){
    weights_mat[,i] <- ProbWeightsInd(alpha[[i]],beta[[i]])
  }
  return(weights_mat)
}

ProbWeightsInd <- function(alpha,beta){
  weights <- numeric(dim(alpha)[1])
  for (i in 1:dim(alpha)[1]){
    weights[i] <- ProbMC(alpha,beta,i)[2]
  }
  return(weights)
}

CalcTran <- function(alpha,beta,act,light,tran,emit_act,emit_light,covar_mat,wake_sleep_corr){
  
  tran_matrix <- matrix(0,2,2)
  len <- dim(act)[1]
  
  for (init_state in 1:2){
    for (new_state in 1:2){
        
      tran_vals <- foreach(ind = 1:length(alpha), .combine = 'c') %dopar% {
        
        alpha_ind <- alpha[[ind]]
        beta_ind <- beta[[ind]]
        likelihood <- logSumExp(c(alpha_ind[len,]))
        
        act_ind <- act[,ind]
        light_ind <- light[,ind]
        
        alpha_ind[1:(len-1),init_state] + 
          beta_ind[2:len,new_state] + 
          log(tran[init_state,new_state]) + 
          unlist(lapply(c(2:len),
                        logClassification,current_state = new_state-1,act = act_ind,light=light_ind,
                        emit_act=emit_act,emit_light=emit_light,covar_ind=covar_mat[ind,], wake_sleep_corr = wake_sleep_corr)) - 
          likelihood
        

        
      }
      tran_matrix[init_state,new_state] <- logSumExp(c(tran_vals))
    }
  }
  # return(tran_matrix)
  
  rowsum_1 <- logSumExp(c(tran_matrix[1,]))
  rowsum_2 <- logSumExp(c(tran_matrix[2,]))
  
  est_tran_matrix <- matrix(0,2,2)
  est_tran_matrix[1,1] <- exp(tran_matrix[1,1] - rowsum_1)
  est_tran_matrix[1,2] <- exp(tran_matrix[1,2] - rowsum_1)
  est_tran_matrix[2,1] <- exp(tran_matrix[2,1] - rowsum_2)
  est_tran_matrix[2,2] <- exp(tran_matrix[2,2] - rowsum_2)
  
  return(est_tran_matrix)
}

#Used for optim, so not currently in use
LogLikeNorm <- function(params,data,weights,activity_bool, covar_mat = NA){
  loglike <- 0
  for (ind in 1:dim(data)[2]){
    
    alpha_ind <- alpha[[ind]]
    ind_likelihood <- logSumExp(c(alpha_ind[dim(alpha_ind)[1],]))
    
    for (i in 1:dim(data)[1]){
      #Take absolute value of sd b/c optim was estimating negative variance
      #Check with mark and paul on this
      
      if (activity_bool){
        mu1 <- sum(covar_mat[ind,] * params[c(1,3,5)])
        sig1 <- sqrt(sum((covar_mat[ind,] * params[c(2,4,6)])^2))
        
        mu2 <- sum(covar_mat[ind,] * params[c(7,9,11)])
        sig2 <- sqrt(sum((covar_mat[ind,] * params[c(8,10,12)])^2))
      } else {
        mu1 <- params[1]
        sig1 <- params[2]
        
        mu2 <- params[3]
        sig2 <- params[4]
      }
      
      if (!is.na(data[i,ind])){
        loglike_norm1 <- log(dnorm(data[i,ind],mean = mu1,sd = abs(sig1)))
        loglike_norm2 <- log(dnorm(data[i,ind],mean = mu2,sd = abs(sig2)))
        
        if (loglike_norm1 < log(.Machine$double.xmin)){loglike_norm1 <- log(.Machine$double.xmin)}
        if (loglike_norm2 < log(.Machine$double.xmin)){loglike_norm2 <- log(.Machine$double.xmin)}
        
        loglike <- loglike + (1-weights[i,ind])*loglike_norm1 + weights[i,ind]*loglike_norm2
        
      }
    }
  }
  return(-loglike)
}

LogLikeBiNorm <- function(params,act,light,weights, covar_mat){
  loglike <- 0
  
  mu_act_wake <- params[1]
  sig_act_wake <- params[4]
  
  mu_act_sleep <- params[7]
  sig_act_sleep <- params[10]
  
  
  mu_light_wake <- params[13]
  mu_light_sleep <- params[14]
  
  sig_light_wake <- params[15]
  sig_light_sleep <- params[16]
  
  wake_corr <- params[17]
  sleep_corr <- params[18]
  
  
  sigma_mat_wake <- diag(c(sig_act_wake,sig_light_wake)) %*% 
    nearPD(matrix(c(1,wake_corr,wake_corr,1),nrow = 2, byrow = T),corr = T, base.matrix = T)$mat %*% 
    diag(c(sig_act_wake,sig_light_wake))
  
  sigma_mat_sleep <- diag(c(sig_act_sleep,sig_light_sleep)) %*% 
    nearPD(matrix(c(1,sleep_corr,sleep_corr,1),nrow = 2, byrow = T),corr = T, base.matrix = T)$mat %*% 
    diag(c(sig_act_sleep,sig_light_sleep))
  
  for (ind in 1:dim(act)[2]){
    covar1_act_wake <- params[2] * covar_mat[ind,2]
    covar2_act_wake <- params[3] * covar_mat[ind,3]
    covar1_act_sleep <- params[8] * covar_mat[ind,2]
    covar2_act_sleep <- params[9] * covar_mat[ind,3]
    
    covar_act_wake <- covar1_act_wake + covar2_act_wake
    covar_act_sleep <- covar1_act_sleep + covar2_act_sleep
    
    # print(sigma_mat_wake)
    
    loglike_sleep <- dmvnorm(cbind(act[,ind] - covar_act_sleep,light[,ind]),
                             mu = c(mu_act_sleep,mu_light_sleep),
                             Sigma = sigma_mat_sleep,
                             log = T)
    
    loglike_wake <- dmvnorm(cbind(act[,ind] - covar_act_wake,light[,ind]),
                            mu = c(mu_act_wake,mu_light_wake),
                            Sigma = sigma_mat_wake,
                            log = T)
    
    
    
    loglike <- loglike + (1-weights[,ind])%*%loglike_wake + weights[,ind]%*%loglike_sleep
    
  }
  return(-loglike)
}


CalcLikelihood <- function(alpha){
  num_obs <- dim(alpha[[1]])[1]
  like_vec <- c()
  for (i in 1:length(alpha)){
    ind_like <- logSumExp(c(alpha[[i]][num_obs,]))
    like_vec <- c(like_vec,ind_like)
  }
  return(sum(like_vec))
  
}

CalcEmpiricTran <- function(mc){
  tran <- matrix(0,2,2)
  for (i in 2:dim(mc)[1]){
    for (j in 1:dim(mc)[2]){
      tran[mc[i-1,j]+1,mc[i,j] + 1] <- tran[mc[i-1,j]+1,mc[i,j] + 1] + 1
    }
  }
  return(tran/rowSums(tran))
}

CalcEmpiricLight <- function(mc,data){
  emit <- matrix(0,2,2)
  mc_vec <- as.vector(mc)
  data_vec <- as.vector(data)
  
  emit_param <- c(mean(data_vec[mc_vec == 0]),
                  sqrt(var(data_vec[mc_vec == 0])),
                  mean(data_vec[mc_vec == 1]),
                  sqrt(var(data_vec[mc_vec == 1])))
  
  emit_mat <- matrix(emit_param, byrow = T,2,2)
  return(emit_mat)
  
}

CalcEmpiricAct <- function(mc,act_array,covar_mat){
  
  state0_norm <- c()
  state0_covar1 <- c()
  state0_covar2 <- c()
  
  state1_norm <- c()
  state1_covar1 <- c()
  state1_covar2 <- c()
  
  
  for (ind in 1:dim(mc)[2]){
    current_act <- act_array[,,ind]
    current_mc <- mc[,ind]
    
    state0_norm <- c(state0_norm,current_act[current_mc == 0,1])
    state1_norm <- c(state1_norm,current_act[current_mc == 1,1])
    
    if (covar_mat[ind,2] == 1){
      state0_covar1 <- c(state0_covar1,current_act[current_mc == 0,2])
      state1_covar1 <- c(state1_covar1,current_act[current_mc == 1,2])
    }
    
    if (covar_mat[ind,3] == 1){
      state0_covar2 <- c(state0_covar2,current_act[current_mc == 0,3])
      state1_covar2 <- c(state1_covar2,current_act[current_mc == 1,3])
    }
  }
  
  emit_act_wake_emptrue <- matrix(c(mean(state0_norm),sqrt(var(state0_norm)),
                                    mean(state0_covar1),sqrt(var(state0_covar1)),
                                    mean(state0_covar2),sqrt(var(state0_covar2))),byrow=T,3)
  
  emit_act_sleep_emptrue  <- matrix(c(mean(state1_norm),sqrt(var(state1_norm)),
                                      mean(state1_covar1),sqrt(var(state1_covar1)),
                                      mean(state1_covar2),sqrt(var(state1_covar2))),byrow=T,3)
  
  
  emit_act_emptrue <- array(0,dim  = c(3,2,2))
  emit_act_emptrue[,,1] <- emit_act_wake_emptrue
  emit_act_emptrue[,,2] <- emit_act_sleep_emptrue
  return(emit_act_emptrue)
  
}

ParamsAct2EmitAct <- function(param_act){
  
  emit_act_placeholder <- matrix(param_act,byrow = T,6)
  emit_act <- abind(emit_act_placeholder[1:3,],emit_act_placeholder[4:6,], along = 3)
  
  return(emit_act)
}

ParamsAct2EmitActBiv <- function(param_act){
  
  
  emit_act1 <- matrix(param_act[1:6],byrow = F,3)
  emit_act2 <- matrix(param_act[7:12],byrow = F,3)
  emit_act <- abind(emit_act1,emit_act2, along = 3)
  
  return(emit_act)
}

EmitAct2ParamsAct <- function(emit_act){
  return(c(as.vector(t(emit_act[,,1])),as.vector(t(emit_act[,,2]))))
}

InduceMissingVec <- function(vec,prob){
  for (i in 1:length(vec)){
    if (runif(1) < prob){vec[i] <- NA}
  }
  return(vec)
}

WeightedSE <- function(model, weights_vec){
  k=length(model$coefficients)
  SSE <- model$residuals^2 %*% weights_vec
  n <- sum(weights_vec)
  return(sqrt(SSE/(n-k+1)))
}

unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

#Gets predicted act without covariates buy subtracting expected value of covariate 
#Can then calculate correlation since we simulate correlation not on covariates
RecoverAct <- function(act,emit_act,covar_mat,weights_mat){
  recovered_act <- act
  for (ind in 1:dim(act)[2]){
    recovered_act[,ind] <- act[,ind] - 
      sum((covar_mat[ind,] * emit_act[,1,1])[2:3]) * (1 - weights_mat[,ind]) -
      sum((covar_mat[ind,] * emit_act[,1,2])[2:3]) * (weights_mat[,ind])
  }
  return(recovered_act)
}

CalcWeightCorr <- function(recovered_act,light,weights_vec){
  recovered_act_vec <- as.vector(recovered_act)
  light_vec <- as.vector(light)
  
  inds_to_keep <- !is.na(light_vec) & !is.na(recovered_act_vec)
  weights_vec <- weights_vec[inds_to_keep]
  recovered_act_vec <- recovered_act_vec[inds_to_keep]
  light_vec <- light_vec[inds_to_keep]
  

  return(c(weightedCorr(recovered_act_vec,light_vec,"Pearson",weights = (1-weights_vec), fast = F),
           weightedCorr(recovered_act_vec,light_vec,"Pearson",weights = (weights_vec), fast = F)))
}

CalcWeightCorrFish <- function(recovered_act,light,weights_mat){
#   fish_z_vec <- c()
#   for (i in 1:dim(recovered_act)[2]){
#     fish_z <- FisherZ(weightedCorr(recovered_act[,i],light[,i],"Pearson",weights = weights_mat[,i]))
#     fish_z_vec <- c(fish_z_vec, fish_z)
#   }
#   sleep_corr <- FisherZInv(mean(fish_z_vec))
#   
#   fish_z_vec <- c()
#   for (i in 1:dim(recovered_act)[2]){
#     fish_z <- FisherZ(weightedCorr(recovered_act[,i],light[,i],"Pearson",weights = (1-weights_mat[,i])))
#     fish_z_vec <- c(fish_z_vec, fish_z)
#   }
#   wake_corr <- FisherZInv(mean(fish_z_vec))
#   
#   return(c(wake_corr,sleep_corr))
}
  
#### User Settings Start Here ####

library(matrixStats)
library(abind)
library(foreach)
library(doParallel)
library(MASS)
library(wCorr)
library(mvtnorm)
library(Matrix)
library(emdbook)


sim_num <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
if (is.na(sim_num)){sim_num <- 14}
# set.seed(sim_num)

#### Set True Parameters ####
init_true <- c(.292,.708)
tran_true <- matrix(c(.846,.154,
                      .074,.926),byrow = T, ncol = 2)

# emit_act_wake_true <- matrix(c(15.923,11.797,
#                                .391,0,
#                                1.549,0),byrow=T,3)
# 
# emit_act_sleep_true <- matrix(c(3.564,7.39,
#                                 .031,0,
#                                 .574,0),byrow=T,3)


emit_act_wake_true <- matrix(c(15.923,11.797,
                               0,0,
                               0,0),byrow=T,3)

emit_act_sleep_true <- matrix(c(3.564,7.39,
                                0,0,
                                0,0),byrow=T,3)


emit_act_true <- abind(emit_act_wake_true,emit_act_sleep_true,along=3)

dimnames(emit_act_true)[[1]] <- c("Norm Act","Covar 1", "Covar 2")
dimnames(emit_act_true)[[2]] <- c("Mean","Std Dev")
dimnames(emit_act_true)[[3]] <- c("Wake","Sleep")
# 
# emit_light_true <- matrix(c(272.956,562.920,
#                             0,3e-50),byrow = T,ncol = 2)

emit_light_true <- matrix(c(272.956,562.920,
                            0,1),byrow = T,ncol = 2)

colnames(emit_light_true) <- c("Mean","Std Dev")
rownames(emit_light_true) <- c("Wake","Sleep")

# wake_sleep_corr_true <- c(.15,0)
# wake_sleep_corr_true <- c(.3,.1)
wake_sleep_corr_true <- c(.8,.8)

#### Simulate True Data ####

day_length <- 96 * 1/2
num_of_people <- 1000

simulated_hmm <- SimulateHMM(day_length,num_of_people,2,init_true,tran_true,emit_act_true,emit_light_true,wake_sleep_corr_true)

mc <- simulated_hmm[[1]]
act <- simulated_hmm[[2]]
light <- simulated_hmm[[3]]
covar_mat <- simulated_hmm[[4]]
act_array <- simulated_hmm[[5]]


tran_true_emp <- CalcEmpiricTran(mc)
init_true_emp <- c(sum(mc[1,] == 0),sum(mc[1,] == 1)) / dim(mc)[2]
emit_light_true_emp <- CalcEmpiricLight(mc,light)
params_light_true_emp <- as.vector(t(emit_light_true_emp))
emit_act_true_emp <- CalcEmpiricAct(mc,act_array,covar_mat)
params_act_true_emp <- EmitAct2ParamsAct(emit_act_true_emp)
wake_sleep_corr_true_emp <- c(cor(as.vector(light)[as.vector(mc) == 0],as.vector(act_array[,1,])[as.vector(mc) == 0]),
                              cor(as.vector(light)[as.vector(mc) == 1],as.vector(act_array[,1,])[as.vector(mc) == 1]))

# act <- apply(act,2,InduceMissingVec, prob = .1)
# light <- apply(light,2,InduceMissingVec, prob = .1)

#### Initialize starting parameters ####
# init  <- c(runif(1,0,.5),0)
# init[2] <- 1 - init[1]
# 
# tran <- matrix(c(runif(1,.85,1),0,0,runif(1,.85,1)), byrow = T, 2)
# tran[1,2] <- 1 - tran[1,1]
# tran[2,1] <- 1 - tran[2,2]
# 
# emit_act_wake <- matrix(c(runif(1,15,19),runif(1,9,13),
#                           runif(1,0,1),runif(1,0,0),
#                           runif(1,1,2),runif(1,0,0)),byrow=T,3)
# 
# emit_act_sleep <- matrix(c(runif(1,0,3),runif(1,2,4),
#                            runif(1,-1,0),runif(1,0,0),
#                            runif(1,0,1/5),runif(1,0,0)),byrow=T,3)
# 
# emit_act <- abind(emit_act_wake,emit_act_sleep,along=3)
# 
# emit_light <- matrix(c(runif(1,225,275),runif(1,300,350),
#                        runif(1,3,7),runif(1,3,7)), byrow = T, 2)
# 
# params_act <- EmitAct2ParamsAct(emit_act)
# params_light <- as.vector(t(emit_light))
# 
# wake_sleep_corr <- c(runif(1,0,.2),runif(1,.1,.3))

#### Load in Real Data ####

# load("WaveHdata.rda")
# 
# act <- t(wave_data[[1]])[2:865,]
# light <- t(wave_data[[2]])[2:865,]
# covar_mat <- as.matrix(wave_data[[3]])
# 
# 
# act <- act[1:96,1:200]
# light <- light[1:96,1:200]
# covar_mat <- covar_mat[1:200,]
# 
# day_length <- dim(act)[1]
# 
# print("Loaded Data")
# 

#### EM #### 

init <- init_true
tran <- tran_true
emit_act <- emit_act_true
emit_light <- emit_light_true
wake_sleep_corr <- wake_sleep_corr_true

# apply(alpha[[2]]+beta[[2]],1,logSumExp)

RepCovarInd <- function(covar_ind){return(t(replicate(day_length, covar_ind)))}
covar_mat_rep <- apply(covar_mat,2,RepCovarInd)

act_cv.df <- data.frame(activity = as.vector(act),
                        cv1 = covar_mat_rep[,2],
                        cv2 = covar_mat_rep[,3])
likelihood_vec <- c()

print("PRE PAR")

print(parallel::detectCores() - 1)
n.cores <- 7
cl <- makeCluster(n.cores)

print("PRE REG")

registerDoParallel(cl)

print("POST REG")

clusterExport(cl,c('ForwardInd','BackwardInd','logClassification','logSumExp','dmvnorm'))

foreach::getDoParRegistered()
foreach::getDoParWorkers()

print("PRE ALPHA")
alpha <- Forward(act,light,init,tran,emit_act,emit_light,covar_mat,wake_sleep_corr)
print("POST ALPHA")
beta <- Backward(act,light,tran,emit_act,emit_light,covar_mat,wake_sleep_corr)

new_likelihood <- CalcLikelihood(alpha)
likelihood_vec <- c(new_likelihood)
likelihood <- -Inf
like_diff <- new_likelihood - likelihood

# break

while(abs(like_diff) > .00000000000001){
  # tic()
  likelihood <- new_likelihood
  
  
  # init <- CalcInit(alpha,beta)
  # tran <- CalcTran(alpha,beta,act,light,tran,emit_act,emit_light,covar_mat,wake_sleep_corr)
  
  weights_mat <- ProbWeights(alpha,beta)
  weights_vec <- as.vector(weights_mat)
  
  
  ########
  #Est normal distributions
  # wake_act_lm <- lm(activity ~cv1+cv2,data = act_cv.df, weights = (1-weights_vec))
  # sleep_act_lm <- lm(activity ~cv1+cv2,data = act_cv.df, weights = weights_vec)
  wake_act_lm <- lm(activity ~1,data = act_cv.df, weights = (1-weights_vec))
  sleep_act_lm <- lm(activity ~1,data = act_cv.df, weights = weights_vec)
  wake_act_sigma <- WeightedSE(wake_act_lm,1-weights_vec[!is.na(act_cv.df$activity)])
  sleep_act_sigma <- WeightedSE(sleep_act_lm,weights_vec[!is.na(act_cv.df$activity)])


  wake_light_lm <- lm(as.vector(light) ~ 1, weights = (1-weights_vec))
  sleep_light_lm <- lm(as.vector(light) ~ 1, weights = (weights_vec))
  wake_light_sigma <- WeightedSE(wake_light_lm,1-weights_vec[!is.na(light)])
  sleep_light_sigma <- WeightedSE(sleep_light_lm,weights_vec[!is.na(light)])
  

  # Y <- as.vector(light)[!is.na(as.vector(light))]
  # W <- diag(c(1-weights_vec))
  # X <- numeric(length(weights_vec)) + 1
  # B <- solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% Y
  # sig <- sqrt(t(Y - (X*B)) %*% W %*% (Y - (X*B)) / sum(W))
  # 
  
  # wake_light_sigma <- summary(wake_light_lm)$sigma
  # sleep_light_sigma <- summary(sleep_light_lm)$sigma

  # wake_act_sigma <- emit_act[1,2,1]
  # sleep_act_sigma <- emit_act[1,2,2]
  # wake_light_sigma <- emit_light[1,2]
  # sleep_light_sigma <- emit_light[2,2]
  # sleep_light_sigma <- 3e-150

  # wake_light_lm[[1]] <- emit_light[1,1]
  # sleep_light_lm[[1]] <- emit_light[2,1]

  params_act <- c(c(rbind(summary(wake_act_lm)$coefficients[,1],c(wake_act_sigma,0,0))),
                  c(rbind(summary(sleep_act_lm)$coefficients[,1],c(sleep_act_sigma,0,0))))
  params_light <- c(wake_light_lm[[1]],wake_light_sigma,sleep_light_lm[[1]],sleep_light_sigma)
  
  
  
  params_act <- c(c(rbind(c(summary(wake_act_lm)$coefficients[,1],0,0),c(wake_act_sigma,0,0))),
                  c(rbind(c(summary(sleep_act_lm)$coefficients[,1],0,0),c(sleep_act_sigma,0,0))))


  # params_act <- optim(params_act, LogLikeNorm, data = act, weights = weights_mat,activity_bool = T, covar_mat = covar_mat)$par
  # params_light <- optim(params_light, LogLikeNorm, data = light, weights = weights_mat,activity_bool = F,method="BFGS")$par


  emit_act <- ParamsAct2EmitAct(params_act)
  emit_light <- matrix(params_light,byrow = T,2)

  emit_act[,2,] <- abs(emit_act[,2,])
  emit_light[,2] <- abs(emit_light[,2])


  recovered_act <- RecoverAct(act,emit_act,covar_mat,weights_mat)
  wake_sleep_corr <- CalcWeightCorr(recovered_act,light,weights_vec)
  # # wake_sleep_corr <- CalcWeightCorrFish(act,light,weights_mat)
  
  ######## Bivariate Optim
  ########
  # #Doesnt work
  # #Need to use nearPD for covar matrix which I think breaks it
  # params <- c(as.vector(emit_act),as.vector(emit_light),wake_sleep_corr)
  # params <- optim(params,LogLikeBiNorm,act = act, light = light, weights = weights_mat, covar_mat = covar_mat,method="BFGS")$par
  # emit_act <- ParamsAct2EmitActBiv(params[1:12])
  # emit_light <-  matrix(params[13:16],byrow = F,2)
  # wake_sleep_corr <- params[17:18]
  ########
  
  
  alpha <- Forward(act,light,init,tran,emit_act,emit_light,covar_mat,wake_sleep_corr)
  beta <- Backward(act,light,tran,emit_act,emit_light,covar_mat,wake_sleep_corr)
  
  new_likelihood <- CalcLikelihood(alpha)
  like_diff <- new_likelihood - likelihood
  print(like_diff)
  likelihood_vec <- c(likelihood_vec,new_likelihood)
  # toc()
  
  # break
}

  

break 

true_params <- list(init_true_emp,tran_true_emp,emit_light_true_emp,emit_act_true_emp,wake_sleep_corr_true_emp)

est_params <- list(init,tran,emit_light,emit_act,wake_sleep_corr)

params_to_save <- list(true_params,est_params,likelihood_vec)
# save(params_to_save,file = paste0("SimAHMM",sim_num,".rda"))


alpha <- Forward(act,light,init,tran,emit_act,emit_light_true,covar_mat,wake_sleep_corr)
beta <- Backward(act,light,tran,emit_act,emit_light_true,covar_mat,wake_sleep_corr)
weights_mat <- ProbWeights(alpha,beta)
init_like_biv <- CalcLikelihood(alpha)
init_normlike_biv <- -LogLikeBiNorm(c(as.vector(emit_act),as.vector(emit_light_true),wake_sleep_corr),act,light,weights_mat,covar_mat)
init_normlike_uni <- -LogLikeNorm(emit_light_true,light,weights_mat,F)

alpha <- Forward(act,light,init,tran,emit_act,emit_light,covar_mat,wake_sleep_corr)
beta <- Backward(act,light,tran,emit_act,emit_light,covar_mat,wake_sleep_corr)
weights_mat <- ProbWeights(alpha,beta)
em_like_biv <- CalcLikelihood(alpha)
em_normlike_biv <- -LogLikeBiNorm(c(as.vector(emit_act),as.vector(emit_light),wake_sleep_corr),act,light,weights_mat,covar_mat)
em_normlike_uni <- -LogLikeNorm(emit_light,light,weights_mat,F)



parallel::stopCluster(cl)
unregister_dopar()


library(MGMM)

x <- cbind(as.vector(act)[as.vector(mc) == 0],as.vector(light)[as.vector(mc) == 0])
y1 <- FitGMM(x, k = 1)

x <- cbind(as.vector(act)[as.vector(mc) == 1],as.vector(light)[as.vector(mc) == 1])
y2 <- FitGMM(x, k = 1)


emit_act_est <- emit_act
emit_act_est[1,1,1] <- y1@Mean[1]
emit_act_est[1,1,2] <- y2@Mean[1]
emit_act_est[1,2,1] <- sqrt(y1@Covariance[1,1])
emit_act_est[1,2,2] <- sqrt(y2@Covariance[1,1])

emit_light_est <- emit_light
emit_light_est[1,1] <- y1@Mean[2]
emit_light_est[2,1] <- y2@Mean[2]
emit_light_est[1,2] <- sqrt(y1@Covariance[2,2])
emit_light_est[2,2] <- sqrt(y2@Covariance[2,2])

wake_sleep_corr_est <- c(y1@Covariance[1,2] / sqrt(y1@Covariance[1,1]) / sqrt(y1@Covariance[2,2]),
                         y2@Covariance[1,2] / sqrt(y2@Covariance[1,1]) / sqrt(y2@Covariance[2,2]))

alpha <- Forward(act,light,init,tran,emit_act,emit_light_est,covar_mat,wake_sleep_corr)
beta <- Backward(act,light,tran,emit_act,emit_light_est,covar_mat,wake_sleep_corr)
est_like <- CalcLikelihood(alpha)
