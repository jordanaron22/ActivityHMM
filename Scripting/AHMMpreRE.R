use_bivar <- T
use_covar <- T

large_sim <- F
set_seed <- T

real_data <- F
epsilon <- .00001


#### Functions ####
SimulateHMM <- function(day_length,num_ind,covar_num,init,params_tran,emit_act,emit_light, wake_sleep_corr, act_light_binom,use_covar){
  
  if (use_covar == F){
    # wake_sleep_corr <- c(0,0)
  }
  
  covar_mat_emis <- matrix(rbinom(num_ind*covar_num,1,0.5),num_ind,covar_num)
  covar_mat_emis <- cbind((numeric(num_ind) + 1),covar_mat_emis)
  activity_array <- array(0,dim  = c(day_length,covar_num+1,num_ind))
  
  covar_mat_tran <- t(rmultinom(num_ind,1,rep(1/6,6)))
  covar_mat_tran <- cbind((numeric(num_ind) + 1),covar_mat_tran[,2:6])
  
  for (ind in 1:num_ind){
    hidden_states <- numeric(day_length)
    activity_mat_ind <- matrix(NA,day_length,covar_num+1)
    light <- numeric(day_length)
    
    act_lod_vec <- numeric(day_length)
    light_lod_vec <- numeric(day_length)
    
    for (i in 1:day_length){
      tran <- Params2Tran(params_tran,i,ChooseTran(covar_mat_tran[ind,]))
      
      if (i == 1) {
        hidden_states[1] <- rbinom(1,1,init[2])
      } else {
        hidden_states[i] <- rbinom(1,1,tran[hidden_states[i-1] + 1,2])
      }
      
      mu_act <- (covar_mat_emis[ind,] * emit_act[,1,hidden_states[i]+1])
      sig_act <- (covar_mat_emis[ind,] * emit_act[,2,hidden_states[i]+1]) 
      
      mu_light <- emit_light[hidden_states[i] + 1,1]
      sig_light <- emit_light[hidden_states[i] + 1,2]
      
      sigma_mat <- matrix(c(sig_act[1]^2,wake_sleep_corr[hidden_states[i] + 1]*sig_act[1]*sig_light,
                            wake_sleep_corr[hidden_states[i] + 1]*sig_act[1]*sig_light, sig_light^2), nrow = 2, byrow=T)
      act_light <- mvrnorm(1,c(mu_act[1],mu_light),sigma_mat)
      
      activity_mat_ind[i,1] <- act_light[1]
      activity_mat_ind[i,2:3] <- mu_act[2:3]
      light[i] <- act_light[2]
      
      if (hidden_states[i] == 1){
        if (rbinom(1,1,act_light_binom[1]) == 1){
          activity_mat_ind[i,1] <- log(epsilon)
          activity_mat_ind[i,2:3] <- 0
          act_lod_vec[i] <- 1
        }
        
        if (rbinom(1,1,act_light_binom[2]) == 1){
          light[i] <- log(epsilon)
          light_lod_vec[i] <- 1
        } 
      }
    }
    
    activity_array[,,ind] <- activity_mat_ind
    activity_vec <- rowSums(activity_mat_ind)
    
    if (ind == 1){
      hidden_states_matrix <- hidden_states
      activity_matrix <- activity_vec
      light_matrix <- light
      act_lod <- act_lod_vec
      light_lod <- light_lod_vec
    } else {
      hidden_states_matrix <- cbind(hidden_states_matrix,hidden_states)
      activity_matrix <- cbind(activity_matrix,activity_vec)
      light_matrix <- cbind(light_matrix,light)
      act_lod <- cbind(act_lod,act_lod_vec)
      light_lod <- cbind(light_lod,light_lod_vec)
    }
    
    
  }
  
  
  return(list(hidden_states_matrix,activity_matrix,light_matrix,covar_mat_emis,covar_mat_tran,activity_array,act_lod,light_lod))
}

logClassification <- function(time,current_state,act,light,emit_act,emit_light,covar_ind,wake_sleep_corr, act_light_binom,BiVar = T){
  mu_act <- emit_act[,1,current_state+1] %*% covar_ind
  sig_act <- sqrt(sum((covar_ind * emit_act[,2,current_state+1])^2))
  
  mu_light <- emit_light[current_state+1,1]
  sig_light <- emit_light[current_state+1,2]
  
  sigma_mat <- matrix(c(sig_act^2,wake_sleep_corr[current_state + 1]*sig_act*sig_light,
                        wake_sleep_corr[current_state + 1]*sig_act*sig_light, sig_light^2), nrow = 2, byrow=T)
  
  if (!is.na(act[time]) & !is.na(light[time])) {
    if (BiVar){
      if (current_state == 0){
        lognorm_dens <- log(mvtnorm::dmvnorm(x = c(act[time],light[time]), mean = c(mu_act,mu_light), sigma = sigma_mat, checkSymmetry = F))
      } else {
        lognorm_dens <- log(
          ((1-act_light_binom[2]) * (light[time]!=log(epsilon)) * (1-act_light_binom[1]) * (act[time]!=log(epsilon)) *
             mvtnorm::dmvnorm(x = c(act[time],light[time]), mean = c(mu_act,mu_light), sigma = sigma_mat, checkSymmetry = F)) +
            
            ((act_light_binom[2] * (light[time]==log(epsilon))) * (1-act_light_binom[1]) * (act[time]!=log(epsilon)) * 
            dnorm(act[time],mu_act,sig_act)) + 
            
            ((1-act_light_binom[2]) * (light[time]!=log(epsilon)) * (act_light_binom[1]) * (act[time]==log(epsilon)) * 
            dnorm(light[time],mu_light,sig_light)) +
            
            ((act_light_binom[2] * (light[time]==log(epsilon))) * (act_light_binom[1]) * (act[time]==log(epsilon)))) 
      }
        
    } else {
      if (current_state == 0){
        lognorm_dens <- log(dnorm(act[time],mu_act,sig_act)) + log(dnorm(light[time],mu_light,sig_light))
      } else {
        lognorm_dens <- log(dnorm(act[time],mu_act,sig_act)) + 
          #need to add binom for act here
          log(
            ((1-act_light_binom[2]) * (light[time]!=log(epsilon)) * dnorm(light[time],mu_light,sig_light))+
              (act_light_binom[2] * (light[time]==log(epsilon))))
      }
        
    }
    
  
  } else if (!is.na(act[time]) & is.na(light[time])){
    if (current_state == 0){
      lognorm_dens <- log(dnorm(act[time],mu_act,sig_act))
    } else {
      lognorm_dens <-log(
        ((1-act_light_binom[1]) * (act[time]!=log(epsilon)) * dnorm(act[time],mu_act,sig_act)) +
          (act_light_binom[1] * (act[time]==log(epsilon))))
    }

  } else if (is.na(act[time]) & !is.na(light[time])){
    if (current_state == 0){
      lognorm_dens <- log(dnorm(light[time],mu_light,sig_light))
    } else {
      lognorm_dens <-log(
        ((1-act_light_binom[2]) * (light[time]!=log(epsilon)) * dnorm(light[time],mu_light,sig_light)) +
          (act_light_binom[2] * (light[time]==log(epsilon))))
    }


  } else {lognorm_dens <- 0}
    

  if (lognorm_dens == -Inf){
    lognorm_dens <- log(.Machine$double.xmin)
  }
  
  return(lognorm_dens)
}

Param2TranHelper <- function(p12,p21){
  tran <- matrix(0,2,2)
  tran[1,2] <- expit(p12)
  tran[1,1] <- 1- tran[1,2]
  tran[2,1] <- expit(p21)
  tran[2,2] <- 1 - tran[2,1]
  return(tran)
}

Params2Tran <- function(params_tran,time,index){
  param_matrix <- matrix(params_tran,ncol=8,nrow=2, byrow = T)
  if (index == 1){
    tran <- Param2TranHelper(param_matrix[1,1]+param_matrix[1,7]*cos(2*pi*time/96)+param_matrix[1,8]*sin(2*pi*time/96),
                             param_matrix[2,1]+param_matrix[2,7]*cos(2*pi*time/96)+param_matrix[2,8]*sin(2*pi*time/96))
  } else {
    tran <- Param2TranHelper(param_matrix[1,1]+param_matrix[1,index]+param_matrix[1,7]*cos(2*pi*time/96)+param_matrix[1,8]*sin(2*pi*time/96),
                             param_matrix[2,1]+param_matrix[2,index]+param_matrix[2,7]*cos(2*pi*time/96)+param_matrix[2,8]*sin(2*pi*time/96))
  }
  return(tran)
}

ChooseTran <- function(covar_tran_bool){
  covar_ind <- which(covar_tran_bool== 1)
  if (length(covar_ind) == 2){
    return(covar_ind[2])
  } else {
    return(1)
  }
}

Backward <- function(act,light,params_tran,emit_act,emit_light,covar_mat_emis,covar_mat_tran,wake_sleep_corr,act_light_binom){
  beta_list <- foreach(ind = 1:dim(act)[2]) %dopar% {
    act_ind <- act[,ind]
    light_ind <- light[,ind]
    BackwardInd(act_ind,light_ind,params_tran,emit_act,emit_light,covar_mat_emis[ind,],covar_mat_tran[ind,],wake_sleep_corr,act_light_binom)
  }
  return(beta_list)
}

BackwardInd <- function(act, light, params_tran, emit_act, emit_light,covar_ind_emis,covar_ind_tran,wake_sleep_corr,act_light_binom) {
  
  n <- length(act)
  beta <- matrix(0, ncol = 2, nrow = n)
  tran_ind <- ChooseTran(covar_ind_tran)
  
  beta[n,1] <- log(1)
  beta[n,2] <- log(1)
  
  
  for (i in (n-1):1){
    
    tran <- Params2Tran(params_tran,i,tran_ind)
    
    #State 0 from 0
    bp_00 <- log(tran[1,1]) + logClassification(i+1,0,act,light,emit_act,emit_light,covar_ind_emis,wake_sleep_corr,act_light_binom) + beta[i+1,1]
    
    #State 1 from 0
    bp_01 <- log(tran[1,2]) + logClassification(i+1,1,act,light,emit_act,emit_light,covar_ind_emis,wake_sleep_corr,act_light_binom) + beta[i+1,2] 
    
    #State 0 from 1
    bp_10 <- log(tran[2,1]) + logClassification(i+1,0,act,light,emit_act,emit_light,covar_ind_emis,wake_sleep_corr,act_light_binom) + beta[i+1,1] 
    
    #State 1 from 1
    bp_11 <- log(tran[2,2]) + logClassification(i+1,1,act,light,emit_act,emit_light,covar_ind_emis,wake_sleep_corr,act_light_binom) + beta[i+1,2] 
    
    
    beta[i,1] <- logSumExp(c(bp_00,bp_01))
    beta[i,2] <- logSumExp(c(bp_10,bp_11))
  }
  
  return(beta)
  
}

Forward <- function(act,light,init,params_tran,emit_act,emit_light,covar_mat_emis,covar_mat_tran,wake_sleep_corr,act_light_binom){
  # alpha_list <- list()
  alpha_list <- foreach(ind = 1:dim(act)[2]) %dopar% {
    act_ind <- act[,ind]
    light_ind <- light[,ind]
  
    ForwardInd(act_ind,light_ind,init,params_tran,emit_act,emit_light,covar_mat_emis[ind,],covar_mat_tran[ind,],wake_sleep_corr,act_light_binom)
  }
  return(alpha_list)
}

ForwardInd <- function(act, light, init, params_tran, emit_act, emit_light,covar_ind_emis,covar_ind_tran,wake_sleep_corr,act_light_binom) {
  alpha <- matrix(0, ncol = 2, nrow=length(act))
  tran_ind <- ChooseTran(covar_ind_tran)
  
  alpha[1,1] <- log(init[1]) + logClassification(1,0,act,light,emit_act,emit_light,covar_ind_emis,wake_sleep_corr,act_light_binom)
  
  alpha[1,2] <- log(init[2]) + logClassification(1,1,act,light,emit_act,emit_light,covar_ind_emis,wake_sleep_corr,act_light_binom)
  
  for (i in 2:length(act)){
    
    tran <- Params2Tran(params_tran,i,tran_ind)
    
    #From state 0 to 0
    fp_00 <- alpha[i-1,1] + log(tran[1,1]) + logClassification(i,0,act,light,emit_act,emit_light,covar_ind_emis,wake_sleep_corr,act_light_binom)
    
    #From state 1 to 0
    fp_10 <- alpha[i-1,2] + log(tran[2,1]) + logClassification(i,0,act,light,emit_act,emit_light,covar_ind_emis,wake_sleep_corr,act_light_binom)
    
    #From state 0 to 1
    fp_01 <- alpha[i-1,1] + log(tran[1,2]) + logClassification(i,1,act,light,emit_act,emit_light,covar_ind_emis,wake_sleep_corr,act_light_binom)
    
    #From state 1 to 1
    fp_11 <- alpha[i-1,2] + log(tran[2,2]) + logClassification(i,1,act,light,emit_act,emit_light,covar_ind_emis,wake_sleep_corr,act_light_binom)
    
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

expit <- function(x){
  return(exp(x) / (1+exp(x)))
}

logit <- function(x){
  return(log(x/(1-x)))
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

CalcEmpiricTran <- function(mc,covar_mat_tran){
  tran <- matrix(0,2,2)
  trans <- list(tran,tran,tran,tran,tran,tran)
  for (ind in 1:dim(mc)[2]){
    tran_ind <- ChooseTran(covar_mat_tran[ind,])
    for (tim in 2:dim(mc)[1]){
      trans[[tran_ind]][mc[tim-1,ind]+1,mc[tim,ind]+1] <- trans[[tran_ind]][mc[tim-1,ind]+1,mc[tim,ind]+1]+1
    }
  }
  for (i in 1:6){
    trans[[i]] <- trans[[i]]/rowSums(trans[[i]])
  }
  return(trans)
}

CalcEmpiricLight <- function(mc,data,light_lod){
  emit <- matrix(0,2,2)
  mc_vec <- as.vector(mc)
  data_vec <- as.vector(data)
  lod_vec <- as.vector(light_lod)
  
  mc_vec <- mc_vec[lod_vec==0]
  data_vec <- data_vec[lod_vec==0]
  
  
  emit_param <- c(mean(data_vec[mc_vec == 0]),
                  sqrt(var(data_vec[mc_vec == 0])),
                  mean(data_vec[mc_vec == 1]),
                  sqrt(var(data_vec[mc_vec == 1])))
  
  emit_mat <- matrix(emit_param, byrow = T,2,2)
  return(emit_mat)
  
}

CalcEmpiricAct <- function(mc,act_array,covar_mat_emis,act_lod){
  
  state0_norm <- c()
  state0_covar1 <- c()
  state0_covar2 <- c()
  
  state1_norm <- c()
  state1_covar1 <- c()
  state1_covar2 <- c()
  
  
  for (ind in 1:dim(mc)[2]){
    current_act <- act_array[,,ind]
    current_mc <- mc[,ind]
    current_lod <- act_lod[,ind]
    
    current_act <- current_act[current_lod==0,]
    current_mc <- current_mc[current_lod==0]
    
    state0_norm <- c(state0_norm,current_act[current_mc == 0,1])
    state1_norm <- c(state1_norm,current_act[current_mc == 1,1])
    
    if (covar_mat_emis[ind,2] == 1){
      state0_covar1 <- c(state0_covar1,current_act[current_mc == 0,2])
      state1_covar1 <- c(state1_covar1,current_act[current_mc == 1,2])
    }
    
    if (covar_mat_emis[ind,3] == 1){
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
  return(sqrt(SSE/(n)))
}

unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

#Gets predicted act without covariates buy subtracting expected value of covariate 
#Can then calculate correlation since we simulate correlation not on covariates
RecoverAct <- function(act,emit_act,covar_mat_emis,weights_mat){
  recovered_act <- act
  for (ind in 1:dim(act)[2]){
    recovered_act[,ind] <- act[,ind] - 
      sum((covar_mat_emis[ind,] * emit_act[,1,1])[2:3]) * (1 - weights_mat[,ind]) -
      sum((covar_mat_emis[ind,] * emit_act[,1,2])[2:3]) * (weights_mat[,ind])
  }
  return(recovered_act)
}

LogLike <- function(params_tran){
  
  alpha <- Forward(act,light,init,params_tran,emit_act,emit_light,covar_mat_emis,covar_mat_tran,wake_sleep_corr,act_light_binom)
  return(-CalcLikelihood(alpha))
}

CalcWeightCorr <- function(act,light,weights_vec,lod_light_weight,lod_act_weight){
  
  act_wake_vec <- as.vector(act)
  act_sleep_vec <- as.vector(act)
  
  act_wake_vec <- act_wake_vec - as.matrix(act_cv.df[,2:3]) %*% emit_act[2:3,1,1]
  act_sleep_vec <- act_sleep_vec - as.matrix(act_cv.df[,2:3]) %*% emit_act[2:3,1,2]
  

  light_vec <- as.vector(light)
  
  inds_to_keep <- !is.na(light_vec) & !is.na(act_wake_vec)
  weights_vec <- weights_vec[inds_to_keep]
  lod_light_weight <- lod_light_weight[inds_to_keep]
  lod_act_weight <- lod_act_weight[inds_to_keep]
  act_wake_vec <- act_wake_vec[inds_to_keep]
  act_sleep_vec <- act_sleep_vec[inds_to_keep]
  light_vec <- light_vec[inds_to_keep]
  

  return(c(weightedCorr(act_wake_vec,light_vec,"Pearson",weights = (1-weights_vec)),
           weightedCorr(act_sleep_vec,light_vec,"Pearson",weights = ((weights_vec) * (1- lod_light_weight)* (1- lod_act_weight)))))
}


CalcTran <- function(alpha,beta,act,light,params_tran,emit_act,emit_light,covar_mat_emis,covar_mat_tran,wake_sleep_corr,act_light_binom, return_grad = F){
  
  len <- dim(act)[1]

  gradient <- matrix(0,2,8)
  hessian <- matrix(0,16,16)
  hessian_vec <- matrix(0,2,8)
  cos_part_vec <- matrix(0,2,6)
  sin_part_vec <- matrix(0,2,6)
  cos_sin_part <- numeric(2)
  
  for (init_state in 1:2){
    for (new_state in 1:2){
      
    tran_vals <- foreach(ind = 1:length(alpha), .combine = 'cbind') %dopar% {
      
        tran_ind <- ChooseTran(covar_mat_tran[ind,])
        tran_vec <- sapply(c(1:(len-1)),FUN = Params2Tran,params_tran = params_tran,index=tran_ind)
        
        #1,1->1 & 2,1->2 & 1,2->3 & 2,2->4
        tran_vec_ind <- init_state * new_state
        if (init_state == 1 & new_state == 2){tran_vec_ind <- 3}
        
        alpha_ind <- alpha[[ind]]
        beta_ind <- beta[[ind]]
        likelihood <- logSumExp(c(alpha_ind[len,]))
        
        act_ind <- act[,ind]
        light_ind <- light[,ind]
        
        exp(alpha_ind[1:(len-1),init_state] + 
          beta_ind[2:len,new_state] + 
          log(tran_vec[tran_vec_ind,]) + 
          unlist(lapply(c(2:len),
                        logClassification,current_state = new_state-1,act = act_ind,light=light_ind,
                        emit_act=emit_act,emit_light=emit_light,covar_ind=covar_mat_emis[ind,], 
                        wake_sleep_corr = wake_sleep_corr,act_light_binom=act_light_binom)) - 
          likelihood) 
      }
    
    
    for (ind in 1:length(alpha)){
      tran_ind <- ChooseTran(covar_mat_tran[ind,])
      tran_vec <- sapply(c(1:(len-1)),FUN = Params2Tran,params_tran = params_tran,index=tran_ind)
      
      if(init_state == 1 & new_state == 1){
        # tran_prime <- -tran[1,2]
        # tran_prime_prime <- -tran[1,1] * tran[1,2]
        tran_prime <- -tran_vec[3,]
        tran_prime_prime <- -tran_vec[3,]*tran_vec[1,]
        
      } else if(init_state == 1 & new_state == 2){ 
        # tran_prime <- tran[1,1]
        # tran_prime_prime <- -tran[1,1] * tran[1,2]
        tran_prime <- tran_vec[1,]
        tran_prime_prime <- -tran_vec[3,]*tran_vec[1,]
        
      } else if(init_state == 2 & new_state == 2){ 
        # tran_prime <- -tran[2,1]
        # tran_prime_prime <- -tran[2,1] * tran[2,2]
        tran_prime <- -tran_vec[2,]
        tran_prime_prime <- -tran_vec[2,] * tran_vec[4,]
        
      } else if(init_state == 2 & new_state == 1){ 
        # tran_prime <- tran[2,2]
        # tran_prime_prime <- -tran[2,1] * tran[2,2]
        tran_prime <- tran_vec[4,]
        tran_prime_prime <- -tran_vec[2,] * tran_vec[4,]
      }
      
      cos_vec <- cos(2*pi*c(1:(len-1))/96)
      sin_vec <- sin(2*pi*c(1:(len-1))/96)
      
      gradient[init_state,tran_ind] <- gradient[init_state,tran_ind] + sum(tran_vals[,ind]*tran_prime)
      if (tran_ind != 1){
        gradient[init_state,1] <- gradient[init_state,1] + sum(tran_vals[,ind]*tran_prime)
      }
      gradient[init_state,7] <- gradient[init_state,7] + sum(tran_vals[,ind]*tran_prime*cos_vec)
      gradient[init_state,8] <- gradient[init_state,8] + sum(tran_vals[,ind]*tran_prime*sin_vec)
      
      
      
      hessian_vec[init_state,tran_ind] <- hessian_vec[init_state,tran_ind] + sum(tran_vals[,ind]*tran_prime_prime)
      cos_part_vec[init_state,tran_ind] <- cos_part_vec[init_state,tran_ind] + sum(tran_vals[,ind]*tran_prime_prime*cos_vec)
      sin_part_vec[init_state,tran_ind] <- sin_part_vec[init_state,tran_ind] + sum(tran_vals[,ind]*tran_prime_prime*sin_vec)
      if (tran_ind != 1){
        hessian_vec[init_state,1] <- hessian_vec[init_state,1] + sum(tran_vals[,ind]*tran_prime_prime)
        cos_part_vec[init_state,1] <- cos_part_vec[init_state,1] + sum(tran_vals[,ind]*tran_prime_prime*cos_vec)
        sin_part_vec[init_state,1] <- sin_part_vec[init_state,1] + sum(tran_vals[,ind]*tran_prime_prime*sin_vec)
      }
      
      hessian_vec[init_state,7] <- hessian_vec[init_state,7] + sum(tran_vals[,ind]*tran_prime_prime*cos_vec^2)
      hessian_vec[init_state,8] <- hessian_vec[init_state,8] + sum(tran_vals[,ind]*tran_prime_prime*sin_vec^2)
      cos_sin_part[init_state] <- cos_sin_part[init_state] + sum(tran_vals[,ind]*tran_prime_prime*cos_vec*sin_vec)
      
      
    }
    
      
    }
  }
  
  gradient <- as.vector(t(gradient))
  
  diag(hessian) <- c(hessian_vec[1,],hessian_vec[2,])
  hessian[1,1:6] <- c(hessian_vec[1,1:6])
  hessian[1:6,1] <- c(hessian_vec[1,1:6])
  hessian[7,1:6] <- cos_part_vec[1,]
  hessian[1:6,7] <- cos_part_vec[1,]
  hessian[8,1:6] <- sin_part_vec[1,]
  hessian[1:6,8] <- sin_part_vec[1,]
  hessian[7,8] <- cos_sin_part[1]
  hessian[8,7] <- cos_sin_part[1]
  
  hessian[9:14,9] <- c(hessian_vec[2,1:6])
  hessian[9,9:14] <- c(hessian_vec[2,1:6])
  hessian[15,9:14] <- cos_part_vec[2,]
  hessian[9:14,15] <- cos_part_vec[2,]
  hessian[16,9:14] <- sin_part_vec[2,]
  hessian[9:14,16] <- sin_part_vec[2,]
  hessian[15,16] <- cos_sin_part[2]
  hessian[16,15] <- cos_sin_part[2]


  if(return_grad){
    return(-gradient)
  } else {
    return(list(-gradient,-hessian))
  }
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
library(numDeriv)
library(tictoc)
library(dplyr)
library(lme4)
# library(emdbook)


#### Set True Parameters ####
init_true <- c(.3,.7)


params_tran_true <- c(-3,.5,1,-1,-1,1,.9,-.8,
                      -2.2,.3,.6,.6,-.6,-.6,-.75,.8)

if (use_covar){
  emit_act_wake_true <- matrix(c(3,1,
                                 1/5,0,
                                 -1/4,0),byrow=T,3)
  
  emit_act_sleep_true <- matrix(c(-1/3,2,
                                  1/3,0,
                                  -1/2,0),byrow=T,3)
  
  
} else {
  emit_act_wake_true <- matrix(c(3,1,
                                 0,0,
                                 0,0),byrow=T,3)
  
  emit_act_sleep_true <- matrix(c(-1/3,2,
                                  0,0,
                                  0,0),byrow=T,3)
}
  

emit_act_true <- abind(emit_act_wake_true,emit_act_sleep_true,along=3)

dimnames(emit_act_true)[[1]] <- c("Norm Act","Covar 1", "Covar 2")
dimnames(emit_act_true)[[2]] <- c("Mean","Std Dev")
dimnames(emit_act_true)[[3]] <- c("Wake","Sleep")
#


emit_light_true <- matrix(c(4,2,
                            0,2),byrow = T,ncol = 2)

colnames(emit_light_true) <- c("Mean","Std Dev")
rownames(emit_light_true) <- c("Wake","Sleep")

#Correlation between activity and light for sleep and wake state
wake_sleep_corr_true <- c(.4,.3)

#Prob of being below detection limit
act_light_binom_true <- c(.05,.7)


#### Simulate True Data ####

if (!real_data){
  
  sim_num <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
  if (is.na(sim_num)){sim_num <- 13}
  if (set_seed){set.seed(sim_num)}
  
  if (large_sim){
    day_length <- 96 * 4
    num_of_people <- 3000
  } else {
    day_length <- 96 * 2
    num_of_people <- 500
  }
  
  
  n <- day_length * num_of_people
  id <- data.frame(SEQN = c(1:num_of_people))
  
  simulated_hmm <- SimulateHMM(day_length,num_of_people,2,init_true,params_tran_true,emit_act_true,emit_light_true,wake_sleep_corr_true,act_light_binom_true, use_covar)
  
  mc <- simulated_hmm[[1]]
  act <- simulated_hmm[[2]]
  light <- simulated_hmm[[3]]
  covar_mat_emis <- simulated_hmm[[4]]
  covar_mat_tran <- simulated_hmm[[5]]
  act_array <- simulated_hmm[[6]]
  act_lod <- simulated_hmm[[7]]
  light_lod <- simulated_hmm[[8]]
  
  
  trans_true_emp <- CalcEmpiricTran(mc,covar_mat_tran)
  init_true_emp <- c(sum(mc[1,] == 0),sum(mc[1,] == 1)) / dim(mc)[2]
  emit_light_true_emp <- CalcEmpiricLight(mc,light,light_lod)
  params_light_true_emp <- as.vector(t(emit_light_true_emp))
  emit_act_true_emp <- CalcEmpiricAct(mc,act_array,covar_mat_emis,act_lod)
  params_act_true_emp <- EmitAct2ParamsAct(emit_act_true_emp)
  wake_sleep_corr_true_emp <- c(cor(as.vector(light)[as.vector(mc) == 0],as.vector(act_array[,1,])[as.vector(mc) == 0]),
                                cor(as.vector(light)[as.vector(mc) == 1 & as.vector(light_lod) == 0 & as.vector(act_lod) == 0],
                                    as.vector(act_array[,1,])[as.vector(mc) == 1 & as.vector(light_lod) == 0 & as.vector(act_lod) == 0]))
  
  act_light_binom_true_emp <- c(sum(act_lod[mc==1]==1)/length(act_lod[mc==1]),
                                sum(light_lod[mc==1]==1)/length(light_lod[mc==1]))
  
  
  #CHECK MISSING AND CORRELATION BEFORE SIM
  act <- apply(act,2,InduceMissingVec, prob = .3)
  light <- apply(light,2,InduceMissingVec, prob = .3)
}



#### Initialize starting parameters ####
init  <- c(runif(1,.1,.5),0)
init[2] <- 1 - init[1]


params_tran <- params_tran_true + runif(16,-.3,.3)

emit_act_wake <- matrix(c(runif(1,2,4),runif(1,1/2,3/2),
                          runif(1,.1,.3),runif(1,0,0),
                          runif(1,-.3,-.2),runif(1,0,0)),byrow=T,3)

emit_act_sleep <- matrix(c(runif(1,-2/3,0),runif(1,1,3),
                           runif(1,.3,.4),runif(1,0,0),
                           runif(1,-.6,-.4),runif(1,0,0)),byrow=T,3)

emit_act <- abind(emit_act_wake,emit_act_sleep,along=3)

emit_light <- matrix(c(runif(1,3,5),runif(1,1,3),
                       runif(1,-1/2,1/2),runif(1,1,3)), byrow = T, 2)

params_act <- EmitAct2ParamsAct(emit_act)
params_light <- as.vector(t(emit_light))

wake_sleep_corr <- c(runif(1,.3,.5),runif(1,.1,.3))
act_light_binom <- c(runif(1,0,.1),runif(1,.4,.6))
#### Load in Real Data ####
if(real_data){
  load("WaveHdata.rda")
  
  sim_num <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
  
  id <- wave_data[[4]]
  
  #CHANGE RACE AND POV
  covar_mat_tran <- matrix(0,ncol = 6, nrow = length(id$poverty))
  covar_mat_tran[,1] <- 1
  for (i in 1:length(id$poverty)){
    covar_mat_tran[i,id$poverty[i]] <- 1
  }
  
  act <- log(wave_data[[1]] + epsilon)
  light <- log(wave_data[[2]] + epsilon)
  
  covar_mat_emis <- as.matrix(wave_data[[3]])
  
  #remove SEQN identifier
  act <- t(act)[2:865,]
  light <- t(light)[2:865,]
  
  day_length <- dim(act)[1]
  
  print("Loaded Data")
  
  
  
  
  init  <- c(1/3,2/3)
  
  params_tran <- c(-3,.01,.01,.01,.01,.01,.01,.01,-2.2,.01,.01,.01,.01,.01,.01,.01)
  
  emit_act_wake <- matrix(c(3,1,
                            0,0,
                            0,0),byrow=T,3)
  
  emit_act_sleep <- matrix(c(-1/3,2,
                             0,0,
                             0,0),byrow=T,3)
  
  emit_act <- abind(emit_act_wake,emit_act_sleep,along=3)
  
  emit_light <- matrix(c(4,2,
                         0,2), byrow = T, 2)
  
  params_act <- EmitAct2ParamsAct(emit_act)
  params_light <- as.vector(t(emit_light))
  
  wake_sleep_corr <- c(.4,.3)
  act_light_binom <- c(.05,.75)
}

#### EM #### 

# init <- init_true_emp
# params_tran <- params_tran_true
# emit_act <- emit_act_true_emp
# emit_light <- emit_light_true_emp
# wake_sleep_corr <- wake_sleep_corr_true_emp
# act_light_binom <- act_light_binom_true_emp

# init <- init_true
# emit_act <- emit_act_true
# emit_light <- emit_light_true
# wake_sleep_corr <- wake_sleep_corr_true
# act_light_binom <- act_light_binom_true

RepCovarInd <- function(covar_ind){return(t(replicate(day_length, covar_ind)))}
covar_mat_emis_rep <- apply(covar_mat_emis,2,RepCovarInd)

id_mat <- apply(as.matrix(id$SEQN),2,RepCovarInd)

act_cv.df <- data.frame(activity = as.vector(act),
                        cv1 = covar_mat_emis_rep[,2],
                        cv2 = covar_mat_emis_rep[,3],
                        SEQN = id_mat)

act_cv.df <- act_cv.df %>% 
  mutate(cv1 = ifelse(activity == log(epsilon),0,cv1)) %>%
  mutate(cv2 = ifelse(activity == log(epsilon),0,cv2))


print("PRE PAR")

if ((parallel::detectCores() - 1) > 8){
  n.cores <- 120
} else {
  n.cores <- 7
}

cl <- makeCluster(n.cores)

print("PRE REG")

registerDoParallel(cl)

print("POST REG")

clusterExport(cl,c('ForwardInd','BackwardInd','logClassification','logSumExp','dmvnorm','ChooseTran', 'epsilon',
                   'Params2Tran','Param2TranHelper','expit'))

foreach::getDoParRegistered()
foreach::getDoParWorkers()

print("PRE ALPHA")
alpha <- Forward(act,light,init,params_tran,emit_act,emit_light,covar_mat_emis,covar_mat_tran,wake_sleep_corr,act_light_binom)
print("POST ALPHA")
beta <- Backward(act,light,params_tran,emit_act,emit_light,covar_mat_emis,covar_mat_tran,wake_sleep_corr,act_light_binom)

break

new_likelihood <- CalcLikelihood(alpha)
likelihood_vec <- c(new_likelihood)
likelihood <- -Inf
like_diff <- new_likelihood - likelihood

# grad_num <- grad(LogLike,params_tran)

while(abs(like_diff) > .001){
  # tic()
  likelihood <- new_likelihood
  
  weights_mat <- ProbWeights(alpha,beta)
  weights_vec <- as.vector(weights_mat)
  
  init <- CalcInit(alpha,beta)
  
  gradhess <- CalcTran(alpha,beta,act,light,params_tran,emit_act,emit_light,covar_mat_emis,covar_mat_tran,wake_sleep_corr,act_light_binom)
  grad <- gradhess[[1]]
  hess <- gradhess[[2]]

  params_tran <- params_tran - solve(hess,grad)
  
  light_vec <- as.vector(light)
  lod_light_weight <- as.numeric(light_vec==log(epsilon))
  act_light_binom[2] <- sum(lod_light_weight,na.rm = T)/sum(weights_vec[!is.na(as.vector(light))])
  
  act_vec <- as.vector(act)
  lod_act_weight <- as.numeric(act_vec==log(epsilon))
  act_light_binom[1] <- sum(lod_act_weight,na.rm = T)/sum(weights_vec[!is.na(as.vector(act))])
  
  x <- lmer(activity ~cv1+cv2 + (1 | SEQN), data = act_cv.df, weights = (1-weights_vec))
  SSE <- resid(x)^2 %*% (1-weights_vec[!is.na(as.vector(act))])
  n <- sum(weights_vec)
  sqrt(SSE/(n))
  
  
  k=length(model$coefficients)
  SSE <- model$residuals^2 %*% weights_vec
  n <- sum(weights_vec)
  return(sqrt(SSE/(n)))
  
  ########
  #Est normal distributions
  if (use_covar){
    wake_act_lm <- lm(activity ~cv1+cv2,data = act_cv.df, weights = (1-weights_vec))
    sleep_act_lm <- lm(activity ~cv1+cv2,data = act_cv.df, weights = weights_vec * (1 - lod_act_weight))
  } else {
    wake_act_lm <- lm(activity ~1,data = act_cv.df, weights = (1-weights_vec))
    sleep_act_lm <- lm(activity ~1,data = act_cv.df, weights = weights_vec * (1 - lod_act_weight))
  }
  wake_act_sigma <- WeightedSE(wake_act_lm,1-weights_vec[!is.na(act_cv.df$activity)])
  sleep_act_sigma <- WeightedSE(sleep_act_lm,weights_vec[!is.na(act_cv.df$activity)] * (1- lod_act_weight[!is.na(act_cv.df$activity)]))

  wake_light_lm <- lm(as.vector(light) ~ 1, weights = (1-weights_vec))
  sleep_light_lm <- lm(as.vector(light) ~ 1, weights = (weights_vec) * (1- lod_light_weight))
  wake_light_sigma <- WeightedSE(wake_light_lm,1-weights_vec[!is.na(light)])
  sleep_light_sigma <- WeightedSE(sleep_light_lm,weights_vec[!is.na(light)]* (1- lod_light_weight[!is.na(light)]))
  
  
  params_light <- c(wake_light_lm[[1]],wake_light_sigma,sleep_light_lm[[1]],sleep_light_sigma)
  if (use_covar){
    params_act <- c(c(rbind(summary(wake_act_lm)$coefficients[,1],c(wake_act_sigma,0,0))),
                    c(rbind(summary(sleep_act_lm)$coefficients[,1],c(sleep_act_sigma,0,0))))
  } else {
    params_act <- c(c(rbind(c(summary(wake_act_lm)$coefficients[,1],0,0),c(wake_act_sigma,0,0))),
                    c(rbind(c(summary(sleep_act_lm)$coefficients[,1],0,0),c(sleep_act_sigma,0,0))))
  }
 

  emit_act <- ParamsAct2EmitAct(params_act)
  emit_light <- matrix(params_light,byrow = T,2)

  emit_act[,2,] <- abs(emit_act[,2,])
  emit_light[,2] <- abs(emit_light[,2])


  if(use_bivar){
    wake_sleep_corr <- CalcWeightCorr(act,light,weights_vec,lod_light_weight,lod_act_weight)
  }

  alpha <- Forward(act,light,init,params_tran,emit_act,emit_light,covar_mat_emis,covar_mat_tran,wake_sleep_corr,act_light_binom)
  beta <- Backward(act,light,params_tran,emit_act,emit_light,covar_mat_emis,covar_mat_tran,wake_sleep_corr,act_light_binom)
  
  new_likelihood <- CalcLikelihood(alpha)
  like_diff <- new_likelihood - likelihood
  print(like_diff)
  likelihood_vec <- c(likelihood_vec,new_likelihood)
}


if (!real_data){
  true_params <- list(init_true_emp,params_tran_true,emit_light_true_emp,emit_act_true_emp,wake_sleep_corr_true_emp,act_light_binom_true_emp)
  est_params <- list(init,params_tran,emit_light,emit_act,wake_sleep_corr,act_light_binom)
  params_to_save <- list(true_params,est_params,likelihood_vec)
} else {
  est_params <- list(init,params_tran,emit_light,emit_act,wake_sleep_corr,act_light_binom)
  params_to_save <- list(est_params,likelihood_vec)
}
  

  
if(large_sim){
  save(params_to_save,file = paste0("AHMM",sim_num,".rda"))
}
parallel::stopCluster(cl)
unregister_dopar()

