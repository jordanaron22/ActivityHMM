large_sim <- F
set_seed <- T

real_data <- F
epsilon <- .00001


#### Functions ####

SimulateHMM <- function(day_length,num_of_people,init,params_tran,emit_act,
                        act_light_binom, re_set,pi_l){
  
  re_vec <- sample(re_set,num_of_people,T,pi_l)
  
  covar_mat_tran <- t(rmultinom(num_of_people,1,rep(1/6,6)))
  covar_mat_tran <- cbind((numeric(num_of_people) + 1),covar_mat_tran[,2:6])
  
  for (ind in 1:num_of_people){
    hidden_states <- numeric(day_length)
    activity <- numeric(day_length)
    
    act_lod_vec <- numeric(day_length)
    
    for (i in 1:day_length){
      tran <- Params2Tran(params_tran,i,ChooseTran(covar_mat_tran[ind,]))
      
      if (i == 1) {
        hidden_states[1] <- rbinom(1,1,init[2])
      } else {
        hidden_states[i] <- rbinom(1,1,tran[hidden_states[i-1] + 1,2])
      }
      
      
      mu_act <- emit_act[hidden_states[i] + 1,1]
      sig_act <- emit_act[hidden_states[i] + 1,2]
      
      activity[i] <-rnorm(1,mu_act,sig_act) 
      
      if (hidden_states[i] == 1){
        
        if (rbinom(1,1,act_light_binom[1]) == 1){
          activity[i] <- log(epsilon)
          act_lod_vec[i] <- 1
        } 
        
      } else {
        activity[i] <- activity[i] + re_vec[ind]
      }
    }
    
    if (ind == 1){
      hidden_states_matrix <- hidden_states
      activity_matrix <- activity
      act_lod <- act_lod_vec
    } else {
      hidden_states_matrix <- cbind(hidden_states_matrix,hidden_states)
      activity_matrix <- cbind(activity_matrix,activity)
      act_lod <- cbind(act_lod,act_lod_vec)
    }
  }
  
  
  return(list(hidden_states_matrix,activity_matrix,covar_mat_tran,act_lod,re_vec))
}

logClassification <- function(time,current_state,act,emit_act,act_light_binom,rande){
  
  mu_act <- emit_act[current_state+1,1]
  sig_act <- emit_act[current_state+1,2]
  
  if (!is.na(act[time])) {
    
    if (current_state == 0){
      lognorm_dens <- log(dnorm(act[time]-rande,mu_act,sig_act)) 
    } else {
      lognorm_dens <- log(((1-act_light_binom[1]) * (act[time]!=log(epsilon)) * dnorm(act[time],mu_act,sig_act))+
              (act_light_binom[1] * (act[time]==log(epsilon))))
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

Backward <- function(act,params_tran,emit_act,covar_mat_tran,act_light_binom,re_set){
  beta_list <- foreach(ind = 1:dim(act)[2]) %:%
    foreach(rand_ind = 1:length(re_set)) %dopar% {
      act_ind <- act[,ind]
      BackwardInd(act_ind,params_tran,emit_act,covar_mat_tran[ind,],act_light_binom,re_set[rand_ind])
    }
  
  beta <- lapply(beta_list, simplify2array)
  
  return(beta)
}

BackwardInd <- function(act, params_tran, emit_act,covar_ind_tran,act_light_binom,rande) {
  
  n <- length(act)
  beta <- matrix(0, ncol = 2, nrow = n)
  tran_ind <- ChooseTran(covar_ind_tran)
  
  beta[n,1] <- log(1)
  beta[n,2] <- log(1)
  
  
  for (i in (n-1):1){
    
    tran <- Params2Tran(params_tran,i+1,tran_ind)
    
    #State 0 from 0
    bp_00 <- log(tran[1,1]) + logClassification(i+1,0,act,emit_act,act_light_binom,rande) + beta[i+1,1]
    
    #State 1 from 0
    bp_01 <- log(tran[1,2]) + logClassification(i+1,1,act,emit_act,act_light_binom,rande) + beta[i+1,2] 
    
    #State 0 from 1
    bp_10 <- log(tran[2,1]) + logClassification(i+1,0,act,emit_act,act_light_binom,rande) + beta[i+1,1] 
    
    #State 1 from 1
    bp_11 <- log(tran[2,2]) + logClassification(i+1,1,act,emit_act,act_light_binom,rande) + beta[i+1,2] 
    
    
    beta[i,1] <- logSumExp(c(bp_00,bp_01))
    beta[i,2] <- logSumExp(c(bp_10,bp_11))
  }
  
  return(beta)
  
}

Forward <- function(act,init,params_tran,emit_act,covar_mat_tran,act_light_binom,rande_vec){
  # alpha_list <- list()
  alpha_list <- foreach(ind = 1:dim(act)[2]) %:%
    foreach(rand_ind = 1:length(rande_vec)) %dopar% {
      act_ind <- act[,ind]
    
      ForwardInd(act_ind,init,params_tran,emit_act,covar_mat_tran[ind,],act_light_binom,rande_vec[rand_ind])
    }
  
  alpha <- lapply(alpha_list, simplify2array)
  
  return(alpha)
}


ForwardInd <- function(act, init, params_tran, emit_act,covar_ind_tran,act_light_binom,rande) {
  alpha <- matrix(0, ncol = 2, nrow=length(act))
  tran_ind <- ChooseTran(covar_ind_tran)
  
  alpha[1,1] <- log(init[1]) + logClassification(1,0,act,emit_act,act_light_binom,rande)
  
  alpha[1,2] <- log(init[2]) + logClassification(1,1,act,emit_act,act_light_binom,rande)
  
  for (i in 2:length(act)){
    
    tran <- Params2Tran(params_tran,i,tran_ind)
    
    #From state 0 to 0
    fp_00 <- alpha[i-1,1] + log(tran[1,1]) + logClassification(i,0,act,emit_act,act_light_binom,rande)
    
    #From state 1 to 0
    fp_10 <- alpha[i-1,2] + log(tran[2,1]) + logClassification(i,0,act,emit_act,act_light_binom,rande)
    
    #From state 0 to 1
    fp_01 <- alpha[i-1,1] + log(tran[1,2]) + logClassification(i,1,act,emit_act,act_light_binom,rande)
    
    #From state 1 to 1
    fp_11 <- alpha[i-1,2] + log(tran[2,2]) + logClassification(i,1,act,emit_act,act_light_binom,rande)
    
    alpha[i,1] <- logSumExp(c(fp_00,fp_10))
    alpha[i,2] <- logSumExp(c(fp_01,fp_11))
  }
  
  return(alpha)
}

CalcInit <- function(alpha, beta,re_mat){

  num_obs <- dim(alpha[[1]][,,1])[1]
  time <- 1
  init_0_vec <- c()
  init_1_vec <- c()
  
  for(ind in 1:length(alpha)){ 
    ind_like <- logSumExp(c(as.vector(t(t(alpha[[ind]][num_obs,,]) + log(re_mat[ind,])))))
    
    init_0_vec <- c(init_0_vec, alpha[[ind]][time,1,] + beta[[ind]][time,1,] + log(re_mat[ind,]) - ind_like)
    
    
    init_1_vec <- c(init_1_vec, alpha[[ind]][time,2,] + beta[[ind]][time,2,] + log(re_mat[ind,]) - ind_like)
  }
  
  
  init_0 <- logSumExp(init_0_vec)
  init_1 <- logSumExp(init_1_vec)
  
  return(exp(c(init_0,init_1) - logSumExp(c(init_0,init_1))))
  
}

expit <- function(x){
  return(exp(x) / (1+exp(x)))
}

logit <- function(x){
  return(log(x/(1-x)))
}

CalcLikelihood <- function(alpha,re_mat){
  num_obs <- dim(alpha[[1]][,,1])[1]
  like_vec <- c()
  #i is number of people
  for (i in 1:length(alpha)){
    ind_like <- logSumExp(c(SumOverREIndTime(alpha,re_mat,i,num_obs)))
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

CalcEmpiricAct <- function(mc,act,act_lod,re_vec_true){
  
  state0_norm <- c()
  
  state1_norm <- c()
  
  
  for (ind in 1:dim(mc)[2]){
    current_act <- act[,ind]
    current_mc <- mc[,ind]
    current_lod <- act_lod[,ind]
    
    current_act <- current_act[current_lod==0]
    current_mc <- current_mc[current_lod==0]
    
    state0_norm <- c(state0_norm,current_act[current_mc == 0]-re_vec_true[ind])
    state1_norm <- c(state1_norm,current_act[current_mc == 1])
  }
  
  emit_act_emp_true <- matrix(c(mean(state0_norm),sqrt(var(state0_norm)),
                               mean(state1_norm),sqrt(var(state1_norm))),byrow=T,ncol = 2)
  
  return(emit_act_emp_true)
  
}

InduceMissingVec <- function(vec,prob){
  for (i in 1:length(vec)){
    if (runif(1) < prob){vec[i] <- NA}
  }
  return(vec)
}

WeightedSE <- function(model, weights_vec){
  SSE <- model$residuals^2 %*% weights_vec
  n <- sum(weights_vec)
  return(sqrt(SSE/(n-1.5)))
}

unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

LogLike <- function(params_tran){
  
  alpha <- Forward(act,init,params_tran,emit_act,covar_mat_tran,act_light_binom)
  return(-CalcLikelihood(alpha))
}

CalcTran <- function(alpha,beta,act,params_tran,emit_act,covar_mat_tran,act_light_binom,re_mat,re_set, return_grad = F){
  
  len <- dim(act)[1]
  
  gradient <- matrix(0,2,8)
  hessian <- matrix(0,16,16)
  hessian_vec <- matrix(0,2,8)
  cos_part_vec <- matrix(0,2,6)
  sin_part_vec <- matrix(0,2,6)
  cos_sin_part <- numeric(2)
  
  for (init_state in 1:2){
    for (new_state in 1:2){
      
      tran_vals_re <- foreach(re_ind = 1:dim(re_mat)[2]) %:% 
        foreach(ind = 1:length(alpha), .combine = 'cbind')%dopar% {
          
          tran_ind <- ChooseTran(covar_mat_tran[ind,])
          tran_vec <- sapply(c(2:(len)),FUN = Params2Tran,params_tran = params_tran,index=tran_ind)
          
          #1,1->1 & 2,1->2 & 1,2->3 & 2,2->4
          tran_vec_ind <- init_state * new_state
          if (init_state == 1 & new_state == 2){tran_vec_ind <- 3}
          
          alpha_ind <- alpha[[ind]]
          beta_ind <- beta[[ind]]
          likelihood <- logSumExp(SumOverREIndTime(alpha,re_mat,ind,len))
          
          act_ind <- act[,ind]
          
          exp(alpha_ind[1:(len-1),init_state,re_ind] + 
                beta_ind[2:len,new_state,re_ind] + 
                log(tran_vec[tran_vec_ind,]) + 
                log(re_mat[ind,re_ind]) + 
                unlist(lapply(c(2:len),
                              logClassification,current_state = new_state-1,act = act_ind,
                              emit_act=emit_act, act_light_binom=act_light_binom,rande = re_set[re_ind])) - 
                likelihood) 
        }
      
      tran_vals_re <- simplify2array(tran_vals_re)
      
      tran_vals <- apply(tran_vals_re, c(1,2), sum)
      
      for (ind in 1:length(alpha)){
        tran_ind <- ChooseTran(covar_mat_tran[ind,])
        tran_vec <- sapply(c(2:(len)),FUN = Params2Tran,params_tran = params_tran,index=tran_ind)
        
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
        
        #Maybe change this back to 1:len-1
        cos_vec <- cos(2*pi*c(2:(len))/96)
        sin_vec <- sin(2*pi*c(2:(len))/96)
        
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

SumOverREIndTime <- function(fb,re_mat,ind,time, add_re = T){
  
  fb_ind <- fb[[ind]]
  re_vec <- re_mat[ind,]
  
  fb_sum <- numeric(2)
  if (add_re){
    fb_sum[1] <- logSumExp(c(fb_ind[time,1,] + log(re_vec)))
    fb_sum[2] <- logSumExp(c(fb_ind[time,2,] + log(re_vec)))
  } else {
    fb_sum[1] <- logSumExp(c(fb_ind[time,1,]))
    fb_sum[2] <- logSumExp(c(fb_ind[time,2,]))
  }
    
  return(fb_sum)
}


Marginalize <- function(alpha,beta){
  alpha_beta <- simplify2array(alpha) + simplify2array(beta)
  #instead of population level re_mat
  #should it be ind level re_weights
  
  for (ind in 1:dim(alpha_beta)[4]){
    for (re_ind in 1:dim(alpha_beta)[3]){
      # alpha_beta[,,re_ind,ind] <- alpha_beta[,,re_ind,ind] + log(re_weights[ind,re_ind])
      alpha_beta[,,re_ind,ind] <- alpha_beta[,,re_ind,ind] + log(re_mat[ind,re_ind])
    }
  }
  
  
  alpha_beta <- apply(alpha_beta,c(1,2,4),logSumExp)
  # alpha_beta <- lapply(seq(dim(alpha_beta)[3]), function(x) alpha_beta[ , , x])
  
  ind_like_mat <- apply(alpha_beta,c(1,3),logSumExp)
  
  alpha_beta[,1,] <- alpha_beta[,1,] - ind_like_mat
  alpha_beta[,2,] <- alpha_beta[,2,] - ind_like_mat
  alpha_beta <- exp(alpha_beta)
  alpha_beta <- alpha_beta[,2,]
  
  return(alpha_beta)
}

CalcProbRE <- function(alpha,re_mat){
  #if using double check this works
  
  len <- dim(alpha[[1]])[1]
  re_len <- dim(alpha[[1]])[3]
  re_weight_vec <- numeric(re_len)
  re_weights <- matrix(0,nrow = length(alpha),ncol = re_len)
  
  
  
  for (ind in 1:length(alpha)){
    for (re_ind in 1:re_len){
      re_weight_vec[re_ind] <- logSumExp(alpha[[ind]][len,,re_ind]) + log(re_mat[ind,re_ind])
    }
    
    for (re_ind in 1:re_len){
      re_weights[ind,re_ind] <- exp(re_weight_vec[re_ind] - logSumExp(c(re_weight_vec)))
    }
    
  }
  
  return(re_weights)
  
}

RepCovarInd <- function(covar_ind){return(t(replicate(day_length, covar_ind)))}
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
library(mclust)
# library(emdbook)


#### Set True Parameters ####
init_true <- c(.3,.7)


params_tran_true <- c(-3,.5,1,-1,-1,1,.9,-.8,
                      -2.2,.3,.6,.6,-.6,-.6,-.75,.8)


emit_act_true <- matrix(c(3,1,
                          -1/3,2),byrow = T,ncol = 2)

colnames(emit_act_true) <- c("Mean","Std Dev")
rownames(emit_act_true) <- c("Wake","Sleep")


#Prob of being below detection limit
act_light_binom_true <- c(.05)

# re_set_true <- c(0)
# re_set_true <- c(-1/2,2)
# re_set_true <- c(-1/2,1,2)
re_set_true <- c(-1,-1/2,1/2,1)

# pi_l_true <- c(1)
# pi_l_true <- c(.6,.4)
# pi_l_true <- c(.5,.4,.1)
pi_l_true <- c(.15,.35,.4,.1)

#### Simulate True Data ####

if (!real_data){
  
  sim_num <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
  if (is.na(sim_num)){sim_num <- 1}
  if (set_seed){set.seed(sim_num)}
  
  if (large_sim){
    day_length <- 96 * 3
    num_of_people <- 2500
  } else {
    day_length <- 96 * (96/96) 
    num_of_people <- 500
  }
  
  
  n <- day_length * num_of_people
  id <- data.frame(SEQN = c(1:num_of_people))
  
  simulated_hmm <- SimulateHMM(day_length,num_of_people,init_true,params_tran_true,emit_act_true,
                               act_light_binom_true, re_set_true, pi_l_true)
  
  mc <- simulated_hmm[[1]]
  act <- simulated_hmm[[2]]
  covar_mat_tran <- simulated_hmm[[3]]
  act_lod <- simulated_hmm[[4]]
  re_vec_true <- simulated_hmm[[5]]
  re_mat_true <- matrix(pi_l_true, 
                            ncol = length(re_set_true), nrow = num_of_people, byrow = T)
  
  trans_true_emp <- CalcEmpiricTran(mc,covar_mat_tran)
  init_true_emp <- c(sum(mc[1,] == 0),sum(mc[1,] == 1)) / dim(mc)[2]
  
  emit_act_true_emp <- CalcEmpiricAct(mc,act,act_lod,re_vec_true)
  
  act_light_binom_true_emp <- c(sum(act_lod[mc==1]==1)/length(act_lod[mc==1]))
  
  re_mat_true_emp <- matrix(table(re_vec_true)/ length(re_vec_true), 
                            ncol = length(re_set_true), nrow = num_of_people, byrow = T)
  
  
  
  #CHECK MISSING AND CORRELATION BEFORE SIM
  # act <- apply(act,2,InduceMissingVec, prob = .3)
}

#### Initialize starting parameters ####
init  <- c(runif(1,.1,.5),0)
init[2] <- 1 - init[1]


params_tran <- params_tran_true + runif(16,-.3,.3)

emit_act <- matrix(c(runif(1,1,5),runif(1,1/2,3/2),
                       runif(1,-2/3,0),runif(1,1,3)), byrow = T, 2)

act_light_binom <- c(runif(1,0,.1))

re_set <- re_set_true + runif(length(re_set_true),-.1,.1)

re_mat_temp <- re_mat_true[1,] + runif(length(re_set_true))
re_mat_temp <- re_mat_temp/sum(re_mat_temp)
re_mat <- matrix(re_mat_temp,ncol = length(re_set_true),nrow = num_of_people, byrow = T)

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
  
  #remove SEQN identifier
  act <- t(act)[2:865,]
  light <- t(light)[2:865,]
  
  day_length <- dim(act)[1]
  
  print("Loaded Data")
  
  init  <- c(1/3,2/3)
  
  params_tran <- c(-3,.01,.01,.01,.01,.01,.01,.01,-2.2,.01,.01,.01,.01,.01,.01,.01)
  
  emit_act <- matrix(c(4,2,
                       0,2), byrow = T, 2)
  
  act_light_binom <- c(.05)
}

#### EM #### 

# init <- init_true_emp
# params_tran <- params_tran_true
# emit_act <- emit_act_true_emp
# act_light_binom <- act_light_binom_true_emp
# re_set <- re_set_true
# re_mat <- re_mat_true_emp
# re_vec <- re_vec_true

# init <- init_true
# params_tran <- params_tran_true
# emit_act <- emit_act_true
# act_light_binom <- act_light_binom_true
# re_set <- re_set_true
# re_mat <- re_mat_true
# re_vec <- re_vec_true

id_mat <- apply(as.matrix(id$SEQN),2,RepCovarInd)

act_cv.df <- data.frame(activity = as.vector(act),
                        SEQN = id_mat,
                        re = apply(as.matrix(re_vec_true),2,RepCovarInd))



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

clusterExport(cl,c('ForwardInd','BackwardInd','logClassification','logSumExp','dnorm','ChooseTran', 'epsilon',
                   'Params2Tran','Param2TranHelper','expit','SumOverREIndTime'))

foreach::getDoParRegistered()
foreach::getDoParWorkers()



print("PRE ALPHA")
alpha <- Forward(act,init,params_tran,emit_act,covar_mat_tran,act_light_binom,re_set)
print("POST ALPHA")
beta <- Backward(act,params_tran,emit_act,covar_mat_tran,act_light_binom,re_set)

# apply(alpha[[2]][,,1]+beta[[2]][,,1],1,logSumExp)

new_likelihood <- CalcLikelihood(alpha,re_mat)
likelihood_vec <- c(new_likelihood)
likelihood <- -Inf
like_diff <- new_likelihood - likelihood

# grad_num <- grad(LogLike,params_tran)


while(abs(like_diff) > .0001){
  likelihood <- new_likelihood
  
  weights_mat <- Marginalize(alpha,beta)
  weights_vec <- as.vector(weights_mat)
  re_prob <- CalcProbRE(alpha,re_mat)
  # weights_vec <- as.vector(mc)
  
  
  ##################
  
  init <- CalcInit(alpha,beta,re_mat)

  gradhess <- CalcTran(alpha,beta,act,params_tran,emit_act,covar_mat_tran,act_light_binom,re_mat,re_set)
  grad <- gradhess[[1]]
  hess <- gradhess[[2]]
  params_tran <- params_tran - solve(hess,grad)
  
  
  ##################
  

  act_vec <- as.vector(act)
  lod_act_weight <- as.numeric(act_vec==log(epsilon))
  act_light_binom[1] <- sum(lod_act_weight,na.rm = T)/sum(weights_vec[!is.na(as.vector(act))])

  
  ##################

  act_cv_em.df <- act_cv.df %>% mutate(weights = (1-weights_vec))
  
  act_mean <- matrix(0,ncol = 3,nrow = length(unique(act_cv_em.df$SEQN)))
  for (i in 1:length(unique(act_cv_em.df$SEQN))){
    act_mean[i,1] <- unique(act_cv_em.df$SEQN)[i]
    act_cv_working <- act_cv_em.df %>% filter(SEQN == unique(act_cv_em.df$SEQN)[i])
    act_mean[i,2] <- weighted.mean(act_cv_working$activity,w = act_cv_working$weights)
    act_mean[i,3] <- sum(act_cv_working$weights)
  }

  # x <- Mclust(act_mean[,2],modelName="V")
  x <- me.weighted(data = act_mean[,2], modelName = "E",
                   z = re_prob, weights = act_mean[,3])

  pi_l <- x$parameters$pro
  re_mat <- matrix(pi_l,ncol = length(pi_l),nrow = length(alpha), byrow = T)



  re_set <- x$parameters$mean - emit_act[1,1]
  re_vec <- x$z %*% re_set
  
  act_mean[,2] <- act_mean[,2] - re_vec
  
  act_cv.df <- act_cv.df %>% mutate(re = apply(as.matrix(re_vec),2,RepCovarInd))
  act_cv_em.df <- act_cv_em.df %>% mutate(re = apply(as.matrix(re_vec),2,RepCovarInd)) 
  
  
  ##################
  
  
  sleep_act_lm <- lm(activity ~1,data = act_cv.df, weights = weights_vec * (1 - lod_act_weight))
  wake_act_sigma <- sqrt(sum((act_cv_em.df$activity - act_cv_em.df$re - weighted.mean(act_mean[,2],w = act_mean[,3]))^2 * 
                          act_cv_em.df$weights) / (sum(act_cv_em.df$weights)-1.5))
  sleep_act_sigma <- WeightedSE(sleep_act_lm,weights_vec[!is.na(act_cv.df$activity)] * (1- lod_act_weight[!is.na(act_cv.df$activity)]))

  
  ##################

  # emit_act[1,1] <- mean(act_mean[,2])
  emit_act[1,1] <- weighted.mean(act_mean[,2],w = act_mean[,3])
  # emit_act[1,1] <- emit_act_true[1,1]

  emit_act[1,2] <- wake_act_sigma
  # emit_act[1,2] <- emit_act_true[1,2]

  emit_act[2,1] <- summary(sleep_act_lm)$coefficients[1]
  # emit_act[2,1] <- emit_act_true[2,1]

  emit_act[2,2] <- sleep_act_sigma
  # emit_act[2,2] <- emit_act_true[2,2]

  emit_act[,2] <- abs(emit_act[,2])

  
  ##################


  alpha <- Forward(act,init,params_tran,emit_act,covar_mat_tran,act_light_binom,re_set)
  beta <- Backward(act,params_tran,emit_act,covar_mat_tran,act_light_binom,re_set)
  
  new_likelihood <- CalcLikelihood(alpha,re_mat)
  like_diff <- new_likelihood - likelihood
  print(like_diff)
  likelihood_vec <- c(likelihood_vec,new_likelihood)
  
}


if (!real_data){
  true_params <- list(init_true_emp,params_tran_true,emit_act_true_emp,act_light_binom_true_emp,re_set_true,re_mat_true_emp[1,])
  est_params <- list(init,params_tran,emit_act,act_light_binom,re_set,pi_l)
  params_to_save <- list(true_params,est_params,likelihood_vec)
} else {
  est_params <- list(init,params_tran,emit_act,act_light_binom,re_set,pi_l)
  params_to_save <- list(est_params,likelihood_vec)
}
  
  
if(large_sim){
  save(params_to_save,file = paste0("AHMM",sim_num,".rda"))
}
parallel::stopCluster(cl)
unregister_dopar()


#####################

