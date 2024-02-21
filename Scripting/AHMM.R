set_seed <- T
sim_num <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
# sim_num <- 1

real_data <- F
sim_covar <- T

epsilon <- 1e-15
# epsilon <- 1e-100
lepsilon <- log(epsilon)



wake_params <- c(2.562106,.5992697)
sleep_params <- c(-1.1387,1.9459)

obs_per_day <- 96

RE_num <- as.numeric(commandArgs(TRUE)[1])
sim_size <- as.numeric(commandArgs(TRUE)[2])
RE_type <- as.character(commandArgs(TRUE)[3])
print(paste("Sim Seed:",sim_num,"Size",sim_size,"RE type",RE_type,"Clust Num:",RE_num))


if(is.na(RE_num)){RE_num <- 8}
if(is.na(sim_size)){sim_size <- 4}
if(is.na(RE_type)){RE_type <- "gamma"}


#### Functions ####

SimulateHMM <- function(day_length,num_of_people,init,params_tran,emit_act,
                        act_light_binom,pi_l,re_set){
  
  if(RE_type == "norm"){re_vec <- rnorm(num_of_people,0,2)}
  if(RE_type == "student"){re_vec <- rt(num_of_people,2)}
  if(RE_type == "stud3nt"){re_vec <- rt(num_of_people,3)}
  if(RE_type == "gamma"){re_vec <- rgamma(num_of_people,2,1)-2}
  if(RE_type == "mix1"){re_vec <- c(rnorm(num_of_people/2,-2,2),rgamma(num_of_people/2,2,1))}
  if(RE_type == "mix2"){re_vec <- c(rnorm(num_of_people/2,1,2),-rgamma(num_of_people/4,4,2),runif(num_of_people/4,-2,2))}
  
  
  tran_covar_num <- 3
  
  covar_mat_tran <- t(rmultinom(num_of_people,1,rep(1/tran_covar_num,tran_covar_num)))
  covar_mat_tran <- cbind((numeric(num_of_people) + 1),covar_mat_tran[,2:tran_covar_num])
  
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

logClassification <- function(time,current_state,act,emit_act,clust_i){
  
  mu_act <- emit_act[current_state+1,1,clust_i]
  sig_act <- emit_act[current_state+1,2,clust_i]
  
  if (!is.na(act[time])) {
    
    if (current_state == 0){
      lognorm_dens <- log(dnorm(act[time],mu_act,sig_act)) 
    } else {
      if (act[time]==lepsilon){
        lognorm_dens <- log(1)
      } else{
        lognorm_dens <- dnorm(act[time],mu_act,sig_act,log = T)
      }
      
    }
    
    
  } else {lognorm_dens <- 0}
  
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

BackwardInd <- function(act_ind, tran_list, emit_act,tran_ind,clust_i) {
  
  n <- length(act_ind)
  beta <- matrix(0, ncol = 2, nrow = n)
  
  beta[n,1] <- log(1)
  beta[n,2] <- log(1)
  
  log_class_0 <- unlist(lapply(c(1:n),logClassification,current_state = 0,act= act_ind, emit_act=emit_act, clust_i = clust_i))
  log_class_1 <- unlist(lapply(c(1:n),logClassification,current_state = 1,act= act_ind, emit_act=emit_act, clust_i = clust_i))
  
  # log_class_0 <- logClassificationC(current_state = 0, act_obs = act_ind,
  #                                   mu = emit_act[1,1,clust_i], sig = emit_act[1,2,clust_i],
  #                                   act_binom = act_light_binom,lod = lepsilon)
  # 
  # log_class_1 <- logClassificationC(current_state = 1, act_obs = act_ind,
  #                                   mu = emit_act[2,1,clust_i], sig = emit_act[2,2,clust_i],
  #                                   act_binom = act_light_binom,lod = lepsilon)
  
  
  for (i in (n-1):1){
    
    #(i-1)%%96+1 is almost modulo 96 but multiples of 96 become 1
    tran <- tran_list[[tran_ind]][[(i+1-1)%%96+1]]
    
    #State 0 from 0
    # bp_00 <- log(tran[1,1]) + logClassification(i+1,0,act_ind,emit_act,act_light_binom,clust_i) + beta[i+1,1]
    bp_00 <- log(tran[1,1]) + log_class_0[i+1] + beta[i+1,1]
    
    #State 1 from 0
    # bp_01 <- log(tran[1,2]) + logClassification(i+1,1,act_ind,emit_act,act_light_binom,clust_i) + beta[i+1,2] 
    bp_01 <- log(tran[1,2]) + log_class_1[i+1] + beta[i+1,2] 
    
    #State 0 from 1
    # bp_10 <- log(tran[2,1]) + logClassification(i+1,0,act_ind,emit_act,act_light_binom,clust_i) + beta[i+1,1] 
    bp_10 <- log(tran[2,1]) + log_class_0[i+1] + beta[i+1,1] 
    
    #State 1 from 1
    # bp_11 <- log(tran[2,2]) + logClassification(i+1,1,act_ind,emit_act,act_light_binom,clust_i) + beta[i+1,2] 
    bp_11 <- log(tran[2,2]) + log_class_1[i+1] + beta[i+1,2] 
    
    
    beta[i,1] <- logSumExp(c(bp_00,bp_01))
    beta[i,2] <- logSumExp(c(bp_10,bp_11))
  }
  
  return(beta)
  
}

ForwardIndAll <- function(act,init,tran_list,emit_act_array,tran_ind_vec, lepsilon, log_sweights_vec){
  n <- dim(act)[2]
  fu <- dim(act)[1]
  alpha <- list()
  for (ind in 1:n){
    tran_ind <- tran_ind_vec[ind]
    emit_ind <- array(emit_act_array[,,ind],dim = c(2,2,1))
    alpha_ind <-ForwardIndC(act[,ind],init,tran_list,emit_ind,tran_ind,0,lepsilon,log_sweights_vec[ind])
    alpha_ind <- array(alpha_ind,c(fu,2,1))
    alpha[[ind]] <- alpha_ind
  }
  return(alpha)
}

BackwardIndAll <- function(act, tran_list, emit_act_array,tran_ind_vec,lepsilon){
  n <- dim(act)[2]
  fu <- dim(act)[1]
  beta <- list()
  for (ind in 1:n){
    tran_ind <- tran_ind_vec[ind]
    emit_ind <- array(emit_act_array[,,ind],dim = c(2,2,1))
    beta_ind <- BackwardIndC(act[,ind],tran_list,emit_ind,tran_ind,0,lepsilon)
    beta_ind <- array(beta_ind,c(fu,2,1))
    beta[[ind]] <- beta_ind
  }
  return(beta)
}

TranByTimeVec <- function(index, params_tran, time_vec){
  return(lapply(time_vec, Params2Tran, params_tran = params_tran,index=index))
}

ForwardInd <- function(act_ind, init, tran_list, emit_act,tran_ind,clust_i) {
  alpha <- matrix(0, ncol = 2, nrow=length(act_ind))
  
  log_class_0 <- unlist(lapply(c(1:length(act_ind)),logClassification,current_state = 0,act= act_ind, emit_act=emit_act, clust_i = clust_i))
  log_class_1 <- unlist(lapply(c(1:length(act_ind)),logClassification,current_state = 1,act= act_ind, emit_act=emit_act, clust_i = clust_i))
  
  # log_class_0 <- logClassificationC(current_state = 0, act_obs = act_ind,
  #                                   mu = emit_act[1,1,clust_i], sig = emit_act[1,2,clust_i],
  #                                   act_binom = act_light_binom,lod = lepsilon)
  # 
  # log_class_1 <- logClassificationC(current_state = 1, act_obs = act_ind, 
  #                                   mu = emit_act[2,1,clust_i], sig = emit_act[2,2,clust_i],
  #                                   act_binom = act_light_binom,lod = lepsilon)
  
  alpha[1,1] <- log(init[1]) + log_class_0[1]
  
  alpha[1,2] <- log(init[2]) + log_class_1[1]
  
  for (i in 2:length(act_ind)){
    
    #(i-1)%%96+1 is almost modulo 96 but multiples of 96 become 1
    tran <- tran_list[[tran_ind]][[(i-1)%%96+1]]
    
    #From state 0 to 0
    fp_00 <- alpha[i-1,1] + log(tran[1,1]) + log_class_0[i]
    
    #From state 1 to 0
    fp_10 <- alpha[i-1,2] + log(tran[2,1]) + log_class_0[i]
    
    #From state 0 to 1
    fp_01 <- alpha[i-1,1] + log(tran[1,2]) + log_class_1[i]
    
    #From state 1 to 1
    fp_11 <- alpha[i-1,2] + log(tran[2,2]) + log_class_1[i]
    
    alpha[i,1] <- logSumExp(c(fp_00,fp_10))
    alpha[i,2] <- logSumExp(c(fp_01,fp_11))
  }
  
  return(alpha)
}


CalcInit <- function(alpha, beta,pi_l,pop_pool = T){
  
  num_obs <- dim(alpha[[1]][,,1])[1]
  time <- 1
  init_0_vec <- matrix(0,length(alpha),length(pi_l))
  init_1_vec <- matrix(0,length(alpha),length(pi_l))
  init_vec <- matrix(0,length(alpha),2)
  
  ind_like_vec <- CalcLikelihoodIndVec(alpha,pi_l)
  
  
  for(ind in 1:length(alpha)){ 
    ind_like <- ind_like_vec[ind]
    
    init_0_vec[ind,] <- alpha[[ind]][time,1,] + beta[[ind]][time,1,] + log(pi_l) - ind_like + log_sweights_vec[ind]
    init_1_vec[ind,] <- alpha[[ind]][time,2,] + beta[[ind]][time,2,] + log(pi_l) - ind_like + log_sweights_vec[ind]
    
    init_vec[ind,] <- exp(c(logSumExp(c(init_0_vec[ind,])),logSumExp(c(init_1_vec[ind,]))) - 
                            logSumExp(c(init_0_vec[ind,],init_1_vec[ind,])))
  }
  
  
  init_0 <- logSumExp(init_0_vec)
  init_1 <- logSumExp(init_1_vec)
  
  if(!pop_pool){return(init_vec)}
  return(exp(c(init_0,init_1) - logSumExp(c(init_0,init_1))))
  
}

expit <- function(x){
  to_ret <- exp(x) / (1+exp(x))
  if (is.na(to_ret)){return(1)}
  return(to_ret)
}

logit <- function(x){
  return(log(x/(1-x)))
}

CalcLikelihood <- function(alpha,pi_l){
  return(sum(CalcLikelihoodIndVec(alpha,pi_l)))
}

CalcLikelihoodIndVec <- function(alpha,pi_l){
  num_obs <- dim(alpha[[1]][,,1])[1]
  like_vec <- c()
  #i is number of people
  for (i in 1:length(alpha)){
    ind_like <- logSumExp(c(SumOverREIndTime(alpha,pi_l,i,num_obs)))
    like_vec <- c(like_vec,ind_like)
  }
  return(like_vec)
}

CalcEmpiricTran <- function(mc,covar_mat_tran){
  
  tran_covar_num <- 3
  
  tran <- matrix(0,2,2)
  trans <- list(tran,tran,tran,tran,tran,tran)
  for (ind in 1:dim(mc)[2]){
    tran_ind <- ChooseTran(covar_mat_tran[ind,])
    for (tim in 2:dim(mc)[1]){
      trans[[tran_ind]][mc[tim-1,ind]+1,mc[tim,ind]+1] <- trans[[tran_ind]][mc[tim-1,ind]+1,mc[tim,ind]+1]+1
    }
  }
  for (i in 1:tran_covar_num){
    trans[[i]] <- trans[[i]]/rowSums(trans[[i]])
  }
  return(trans)
}

CalcEmpiricAct <- function(mc,act,act_lod,clust_ind_true){
  
  
  state0_list <- vector(mode = "list", length = length(unique(clust_ind_true)))
  state1_norm <- c()
  
  
  for (ind in 1:dim(mc)[2]){
    current_act <- act[,ind]
    current_mc <- mc[,ind]
    current_lod <- act_lod[,ind]
    clust <- clust_ind_true[ind]
    
    current_act <- current_act[current_lod==0]
    current_mc <- current_mc[current_lod==0]
    
    state0_list[[clust]] <- c(state0_list[[clust]],current_act[current_mc == 0])
    state1_norm <- c(state1_norm,current_act[current_mc == 1])
  }
  
  stnrd_state0 <- unlist(mapply('-',state0_list,unlist(lapply(state0_list,mean)), SIMPLIFY= FALSE))
  
  
  emit_act_emp_true <- array(dim = c(2,2,length(mean_set_true)))
  emit_act_emp_true[1,1,] <- unlist(lapply(state0_list,mean))
  emit_act_emp_true[1,2,] <- sqrt(var(stnrd_state0))
  emit_act_emp_true[2,1,] <- mean(state1_norm)
  emit_act_emp_true[2,2,] <- sqrt(var(state1_norm))
  
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
  
  alpha <- Forward(act,init,params_tran,emit_act,covar_mat_tran)
  return(-CalcLikelihood(alpha))
}

Params2TranVector <- function(index,len,params_tran){
  return(sapply(c(2:(len)),FUN = Params2Tran,params_tran = params_tran,index=index))
}

Params2TranVectorT <- function(index,len,params_tran){
  return(t(sapply(c(2:(len)),FUN = Params2Tran,params_tran = params_tran,index=index)))
}

IndLike <- function(alpha,pi_l,ind,len){
  likelihood <- logSumExp(SumOverREIndTime(alpha,pi_l,ind,len))
  return(likelihood)
}


CalcTran <- function(alpha,beta,act,params_tran,emit_act,covar_mat_tran,pi_l, return_grad = F){
  
  len <- dim(act)[1]
  
  gradient <- matrix(0,2,8)
  hessian <- matrix(0,16,16)
  hessian_vec <- matrix(0,2,8)
  cos_part_vec <- matrix(0,2,6)
  sin_part_vec <- matrix(0,2,6)
  cos_sin_part <- numeric(2)
  
  
  # tran_mat <- lapply(c(1:dim(covar_mat_tran)[2]),TranByTimeVec, params_tran = params_tran, time_vec = c(1:obs_per_day))
  tran_list_mat <- lapply(c(1:dim(covar_mat_tran)[2]),Params2TranVector, len = len, params_tran = params_tran)
  tran_ind_vec <- apply(covar_mat_tran,1,ChooseTran)
  ind_like_vec <- unlist(lapply(c(1:length(alpha)),IndLike,alpha = alpha, pi_l = pi_l, len = len))
  
  for (init_state in 1:2){
    for (new_state in 1:2){
      
      tran_vals_re <- foreach(re_ind = 1:dim(emit_act)[3]) %:% 
        foreach(ind = 1:length(alpha), .combine = 'cbind')%dopar% {
          
          tran_ind <- tran_ind_vec[ind]
          # tran_vec <- sapply(c(2:(len)),FUN = Params2Tran,params_tran = params_tran,index=tran_ind)
          tran_vec <- tran_list_mat[[tran_ind]]
          
          #1,1->1 & 2,1->2 & 1,2->3 & 2,2->4
          tran_vec_ind <- init_state * new_state
          if (init_state == 1 & new_state == 2){tran_vec_ind <- 3}
          
          alpha_ind <- alpha[[ind]]
          beta_ind <- beta[[ind]]
          likelihood <- ind_like_vec[ind]
          
          act_ind <- act[,ind]
          
          class_vec <- logClassificationC(current_state = new_state-1,act_obs = act_ind[-1],
                                          mu = emit_act[new_state,1,re_ind],sig = emit_act[new_state,2,re_ind],
                                          lod = lepsilon)
          
          exp(alpha_ind[1:(len-1),init_state,re_ind] + 
                beta_ind[2:len,new_state,re_ind] + 
                log(tran_vec[tran_vec_ind,]) + 
                log(pi_l[re_ind]) + 
                class_vec - 
                likelihood) 
        }
      
      
      
      
      
      
      
      tran_vals_re <- simplify2array(tran_vals_re)
      
      tran_vals <- apply(tran_vals_re, c(1,2), sum)
      
      for (ind in 1:length(alpha)){
        
        tran_ind <- tran_ind_vec[ind]
        tran_vec <- tran_list_mat[[tran_ind]]
        
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

CalcTranC <- function(alpha,beta,act,params_tran,emit_act,covar_mat_tran,pi_l){
  
  len <- dim(act)[1]
  
  gradient <- matrix(0,2,8)
  hessian <- matrix(0,16,16)
  hessian_vec <- matrix(0,2,8)
  cos_part_vec <- matrix(0,2,6)
  sin_part_vec <- matrix(0,2,6)
  cos_sin_part <- numeric(2)
  
  
  # tran_mat <- lapply(c(1:dim(covar_mat_tran)[2]),TranByTimeVec, params_tran = params_tran, time_vec = c(1:obs_per_day))
  tran_list_mat <- lapply(c(1:dim(covar_mat_tran)[2]),Params2TranVectorT, len = len, params_tran = params_tran)
  tran_ind_vec <- apply(covar_mat_tran,1,ChooseTran)
  ind_like_vec <- unlist(lapply(c(1:length(alpha)),IndLike,alpha = alpha, pi_l = pi_l, len = len))
  
  
  tran_vals_re_00 <- CalcTranHelperC(0,0,
                                     act,tran_list_mat, tran_ind_vec, emit_act, ind_like_vec, alpha, beta,lepsilon, pi_l)
  tran_vals_re_01 <- CalcTranHelperC(0,1,
                                     act,tran_list_mat, tran_ind_vec, emit_act, ind_like_vec, alpha, beta,lepsilon, pi_l)
  tran_vals_re_10 <- CalcTranHelperC(1,0,
                                     act,tran_list_mat, tran_ind_vec, emit_act, ind_like_vec, alpha, beta,lepsilon, pi_l)
  tran_vals_re_11 <- CalcTranHelperC(1,1,
                                     act,tran_list_mat, tran_ind_vec, emit_act, ind_like_vec, alpha, beta,lepsilon, pi_l)
  
  for (init_state in 1:2){
    for (new_state in 1:2){
      
      
      
      if (init_state == 1 & new_state == 1){tran_vals_re <- tran_vals_re_00}
      if (init_state == 1 & new_state == 2){tran_vals_re <- tran_vals_re_01}
      if (init_state == 2 & new_state == 1){tran_vals_re <- tran_vals_re_10}
      if (init_state == 2 & new_state == 2){tran_vals_re <- tran_vals_re_11}
      
      
      tran_vals <- apply(tran_vals_re, c(1,2), sum)
      
      for (ind in 1:length(alpha)){
        
        tran_ind <- tran_ind_vec[ind]
        tran_vec <- tran_list_mat[[tran_ind]]
        
        if(init_state == 1 & new_state == 1){
          # tran_prime <- -tran[1,2]
          # tran_prime_prime <- -tran[1,1] * tran[1,2]
          tran_prime <- -tran_vec[,3]
          tran_prime_prime <- -tran_vec[,3]*tran_vec[,1]
          
        } else if(init_state == 1 & new_state == 2){ 
          # tran_prime <- tran[1,1]
          # tran_prime_prime <- -tran[1,1] * tran[1,2]
          tran_prime <- tran_vec[,1]
          tran_prime_prime <- -tran_vec[,3]*tran_vec[,1]
          
        } else if(init_state == 2 & new_state == 2){ 
          # tran_prime <- -tran[2,1]
          # tran_prime_prime <- -tran[2,1] * tran[2,2]
          tran_prime <- -tran_vec[,2]
          tran_prime_prime <- -tran_vec[,2] * tran_vec[,4]
          
        } else if(init_state == 2 & new_state == 1){ 
          # tran_prime <- tran[2,2]
          # tran_prime_prime <- -tran[2,1] * tran[2,2]
          tran_prime <- tran_vec[,4]
          tran_prime_prime <- -tran_vec[,2] * tran_vec[,4]
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
  
  grad <- -gradient
  hess <- -hessian
  
  if (sim_covar == T){
    null_inds <- grad != 0
    params_tran[null_inds] <- params_tran[null_inds] - solve(hess[null_inds,null_inds],grad[null_inds])
  } else {
    params_tran <- params_tran - solve(hess,grad)
  }
  
  return(params_tran)
}

CalcTranInd <- function(alpha,beta,act,params_tran,emit_act_array,covar_mat_tran,pi_l){
  
  len <- dim(act)[1]
  
  gradient <- matrix(0,2,8)
  hessian <- matrix(0,16,16)
  hessian_vec <- matrix(0,2,8)
  cos_part_vec <- matrix(0,2,6)
  sin_part_vec <- matrix(0,2,6)
  cos_sin_part <- numeric(2)
  
  # tran_mat <- lapply(c(1:dim(covar_mat_tran)[2]),TranByTimeVec, params_tran = params_tran, time_vec = c(1:obs_per_day))
  tran_list_mat <- lapply(c(1:dim(covar_mat_tran)[2]),Params2TranVectorT, len = len, params_tran = params_tran)
  tran_ind_vec <- apply(covar_mat_tran,1,ChooseTran)
  ind_like_vec <- unlist(lapply(c(1:length(alpha)),IndLike,alpha = alpha, pi_l = pi_l, len = len))
  
  
  tran_vals_re_00 <- CalcTranIndHelperC(0,0,
                                        act,tran_list_mat, tran_ind_vec, emit_act_array, ind_like_vec, alpha, beta,lepsilon, pi_l)
  tran_vals_re_01 <- CalcTranIndHelperC(0,1,
                                        act,tran_list_mat, tran_ind_vec, emit_act_array, ind_like_vec, alpha, beta,lepsilon, pi_l)
  tran_vals_re_10 <- CalcTranIndHelperC(1,0,
                                        act,tran_list_mat, tran_ind_vec, emit_act_array, ind_like_vec, alpha, beta,lepsilon, pi_l)
  tran_vals_re_11 <- CalcTranIndHelperC(1,1,
                                        act,tran_list_mat, tran_ind_vec, emit_act_array, ind_like_vec, alpha, beta,lepsilon, pi_l)
  
  
  # # tran_mat <- lapply(c(1:dim(covar_mat_tran)[2]),TranByTimeVec, params_tran = params_tran, time_vec = c(1:obs_per_day))
  # tran_list_mat <- lapply(c(1:dim(covar_mat_tran)[2]),Params2TranVector, len = len, params_tran = params_tran)
  # tran_ind_vec <- apply(covar_mat_tran,1,ChooseTran)
  # ind_like_vec <- unlist(lapply(c(1:length(alpha)),IndLike,alpha = alpha, pi_l = pi_l, len = len))
  
  for (init_state in 1:2){
    for (new_state in 1:2){
      
      
      
      if (init_state == 1 & new_state == 1){tran_vals_re <- tran_vals_re_00}
      if (init_state == 1 & new_state == 2){tran_vals_re <- tran_vals_re_01}
      if (init_state == 2 & new_state == 1){tran_vals_re <- tran_vals_re_10}
      if (init_state == 2 & new_state == 2){tran_vals_re <- tran_vals_re_11}
      
      
      
      tran_vals <- apply(tran_vals_re, c(1,2), sum)
      
      for (ind in 1:length(alpha)){
        
        tran_ind <- tran_ind_vec[ind]
        tran_vec <- tran_list_mat[[tran_ind]]
        
        if(init_state == 1 & new_state == 1){
          # tran_prime <- -tran[1,2]
          # tran_prime_prime <- -tran[1,1] * tran[1,2]
          tran_prime <- -tran_vec[,3]
          tran_prime_prime <- -tran_vec[,3]*tran_vec[,1]
          
        } else if(init_state == 1 & new_state == 2){ 
          # tran_prime <- tran[1,1]
          # tran_prime_prime <- -tran[1,1] * tran[1,2]
          tran_prime <- tran_vec[,1]
          tran_prime_prime <- -tran_vec[,3]*tran_vec[,1]
          
        } else if(init_state == 2 & new_state == 2){ 
          # tran_prime <- -tran[2,1]
          # tran_prime_prime <- -tran[2,1] * tran[2,2]
          tran_prime <- -tran_vec[,2]
          tran_prime_prime <- -tran_vec[,2] * tran_vec[,4]
          
        } else if(init_state == 2 & new_state == 1){ 
          # tran_prime <- tran[2,2]
          # tran_prime_prime <- -tran[2,1] * tran[2,2]
          tran_prime <- tran_vec[,4]
          tran_prime_prime <- -tran_vec[,2] * tran_vec[,4]
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
  
  
  grad <- -gradient
  hess <- -hessian
  
  if (sim_covar == T){
    null_inds <- grad != 0
    params_tran[null_inds] <- params_tran[null_inds] - solve(hess[null_inds,null_inds],grad[null_inds])
  } else {
    params_tran <- params_tran - solve(hess,grad)
  }
  
  return(params_tran)
  
  
}

CalcTranBothC <- function(alpha,beta,act,params_tran,emit_act,covar_mat_tran,pi_l,pop_pool){
  
  len <- dim(act)[1]
  
  gradient <- matrix(0,2,8)
  hessian <- matrix(0,16,16)
  hessian_vec <- matrix(0,2,8)
  cos_part_vec <- matrix(0,2,6)
  sin_part_vec <- matrix(0,2,6)
  cos_sin_part <- numeric(2)
  
  
  # tran_mat <- lapply(c(1:dim(covar_mat_tran)[2]),TranByTimeVec, params_tran = params_tran, time_vec = c(1:obs_per_day))
  tran_list_mat <- lapply(c(1:dim(covar_mat_tran)[2]),Params2TranVectorT, len = len, params_tran = params_tran)
  tran_ind_vec <- apply(covar_mat_tran,1,ChooseTran)
  ind_like_vec <- unlist(lapply(c(1:length(alpha)),IndLike,alpha = alpha, pi_l = pi_l, len = len))
  
  if (pop_pool){
    tran_vals_re_00 <- CalcTranHelperC(0,0,
                                       act,tran_list_mat, tran_ind_vec, emit_act, ind_like_vec, alpha, beta,lepsilon, pi_l)
    tran_vals_re_01 <- CalcTranHelperC(0,1,
                                       act,tran_list_mat, tran_ind_vec, emit_act, ind_like_vec, alpha, beta,lepsilon, pi_l)
    tran_vals_re_10 <- CalcTranHelperC(1,0,
                                       act,tran_list_mat, tran_ind_vec, emit_act, ind_like_vec, alpha, beta,lepsilon, pi_l)
    tran_vals_re_11 <- CalcTranHelperC(1,1,
                                       act,tran_list_mat, tran_ind_vec, emit_act, ind_like_vec, alpha, beta,lepsilon, pi_l)
  } else {
    tran_vals_re_00 <- CalcTranIndHelperC(0,0,
                                          act,tran_list_mat, tran_ind_vec, emit_act, ind_like_vec, alpha, beta,lepsilon, pi_l)
    tran_vals_re_01 <- CalcTranIndHelperC(0,1,
                                          act,tran_list_mat, tran_ind_vec, emit_act, ind_like_vec, alpha, beta,lepsilon, pi_l)
    tran_vals_re_10 <- CalcTranIndHelperC(1,0,
                                          act,tran_list_mat, tran_ind_vec, emit_act, ind_like_vec, alpha, beta,lepsilon, pi_l)
    tran_vals_re_11 <- CalcTranIndHelperC(1,1,
                                          act,tran_list_mat, tran_ind_vec, emit_act, ind_like_vec, alpha, beta,lepsilon, pi_l)
    
  }
  
  for (init_state in 1:2){
    for (new_state in 1:2){
      
      
      
      if (init_state == 1 & new_state == 1){tran_vals_re <- tran_vals_re_00}
      if (init_state == 1 & new_state == 2){tran_vals_re <- tran_vals_re_01}
      if (init_state == 2 & new_state == 1){tran_vals_re <- tran_vals_re_10}
      if (init_state == 2 & new_state == 2){tran_vals_re <- tran_vals_re_11}
      
      
      tran_vals <- apply(tran_vals_re, c(1,2), sum)
      
      for (ind in 1:length(alpha)){
        
        tran_ind <- tran_ind_vec[ind]
        tran_vec <- tran_list_mat[[tran_ind]]
        
        if(init_state == 1 & new_state == 1){
          # tran_prime <- -tran[1,2]
          # tran_prime_prime <- -tran[1,1] * tran[1,2]
          tran_prime <- -tran_vec[,3]
          tran_prime_prime <- -tran_vec[,3]*tran_vec[,1]
          
        } else if(init_state == 1 & new_state == 2){ 
          # tran_prime <- tran[1,1]
          # tran_prime_prime <- -tran[1,1] * tran[1,2]
          tran_prime <- tran_vec[,1]
          tran_prime_prime <- -tran_vec[,3]*tran_vec[,1]
          
        } else if(init_state == 2 & new_state == 2){ 
          # tran_prime <- -tran[2,1]
          # tran_prime_prime <- -tran[2,1] * tran[2,2]
          tran_prime <- -tran_vec[,2]
          tran_prime_prime <- -tran_vec[,2] * tran_vec[,4]
          
        } else if(init_state == 2 & new_state == 1){ 
          # tran_prime <- tran[2,2]
          # tran_prime_prime <- -tran[2,1] * tran[2,2]
          tran_prime <- tran_vec[,4]
          tran_prime_prime <- -tran_vec[,2] * tran_vec[,4]
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
  
  grad <- -gradient
  hess <- -hessian
  
  if (sim_covar == T){
    null_inds <- grad != 0
    params_tran[null_inds] <- params_tran[null_inds] - solve(hess[null_inds,null_inds],grad[null_inds])
  } else {
    params_tran <- params_tran - solve(hess,grad)
  }
  
  return(params_tran)
}



SolveCatch <- function(block_ind_hess,block_ind_grad) {
  tryCatch(
    {
      solve(block_ind_hess,block_ind_grad)
    },
    error = function(cond) {
      # message("Non-Invertible Matrix")
      numeric(length(block_ind_grad))
    },
    warning = function(cond) {
      NULL
    },
    finally = {}
  )
}

SumOverREIndTime <- function(fb,pi_l,ind,time, add_re = T){
  
  fb_ind <- fb[[ind]]
  
  fb_sum <- numeric(2)
  if (add_re){
    fb_sum[1] <- logSumExp(c(fb_ind[time,1,] + log(pi_l)))
    fb_sum[2] <- logSumExp(c(fb_ind[time,2,] + log(pi_l)))
  } else {
    fb_sum[1] <- logSumExp(c(fb_ind[time,1,]))
    fb_sum[2] <- logSumExp(c(fb_ind[time,2,]))
  }
  
  return(fb_sum)
}

CondMarginalize <- function(alpha,beta,pi_l){
  alpha_beta <- simplify2array(alpha) + simplify2array(beta)
  
  
  for (ind in 1:dim(alpha_beta)[4]){
    for (re_ind in 1:dim(alpha_beta)[3]){
      # alpha_beta[,,re_ind,ind] <- alpha_beta[,,re_ind,ind] + log(re_weights[ind,re_ind])
      alpha_beta[,,re_ind,ind] <- alpha_beta[,,re_ind,ind] + log(pi_l[re_ind])
    }
  }
  
  ind_like_mat <- apply(alpha_beta,c(1,4),logSumExp)
  
  weight_array <- array(0, dim = c(dim(alpha_beta)[1],dim(alpha_beta)[4],dim(alpha_beta)[3]))
  for (ind in 1:dim(alpha_beta)[4]){
    for (t in 1:dim(alpha_beta)[1]){
      weight_array[t,ind,] <- alpha_beta[t,1,,ind] - ind_like_mat[t,ind]
    }
  }
  
  return(weight_array)
}


CalcProbRE <- function(alpha,pi_l){
  
  len <- dim(alpha[[1]])[1]
  re_len <- dim(alpha[[1]])[3]
  re_weight_vec <- numeric(re_len)
  re_weights <- matrix(0,nrow = length(alpha),ncol = re_len)
  
  
  
  for (ind in 1:length(alpha)){
    for (re_ind in 1:re_len){
      re_weights[ind,re_ind] <- logSumExp(alpha[[ind]][len,,re_ind]) + log(pi_l[re_ind])
    }
    re_weights[ind,] <- exp(re_weights[ind,] - logSumExp(c(re_weights[ind,])))
    
  }
  
  return(re_weights)
  
}

RepCovarInd <- function(covar_ind){return(t(replicate(day_length, covar_ind)))}

ViterbiInd <- function(ind, RE_num){
  
  params_tran_ind <- params_tran
  if (RE_num == 0){
    emit_ind <- array(emit_act_array[,,ind],dim = c(2,2,1))
  } else {
    emit_ind <- emit_act
  }
  
  clust_ind <- which.max(re_prob[ind,])
  
  
  viterbi_mat <- matrix(NA,2,day_length)
  viterbi_mat[1,1] <- log(init[1]) + logClassification(1,0,act[,ind],emit_ind,clust_ind)
  viterbi_mat[2,1] <- log(init[2]) + logClassification(1,1,act[,ind],emit_ind,clust_ind)
  
  viterbi_ind_mat <- matrix(NA,2,day_length)
  
  tran_ind <- ChooseTran(covar_mat_tran[ind,])
  if(RE_num == 0) {tran_ind <- 1}
  
  for (time in 2:day_length){
    
    tran <- Params2Tran(params_tran_ind,time,tran_ind)
    
    viterbi_mat[1,time] <- logClassification(time,0,act[,ind],emit_ind,clust_ind)+ 
      max(viterbi_mat[1,time-1] + log(tran[1,1]),
          viterbi_mat[2,time-1] + log(tran[2,1]))
    
    
    viterbi_mat[2,time] <- logClassification(time,1,act[,ind],emit_ind,clust_ind) + 
      max(viterbi_mat[1,time-1] + log(tran[1,2]),
          viterbi_mat[2,time-1] + log(tran[2,2]))
    
    
    viterbi_ind_mat[1,time] <-  which.max(c(viterbi_mat[1,time-1] + log(tran[1,1]),
                                            viterbi_mat[2,time-1] + log(tran[2,1])))
    
    
    viterbi_ind_mat[2,time] <- which.max(c(viterbi_mat[1,time-1] + log(tran[1,2]),
                                           viterbi_mat[2,time-1] + log(tran[2,2])))
    
    
  }
  
  decoded_mc <- c(which.max(viterbi_mat[,time]))
  for(time in day_length:2){
    decoded_mc <- c(viterbi_ind_mat[decoded_mc[1],time],decoded_mc)
  }
  
  return(decoded_mc-1)
} 

if(sim_covar){tran_covar_num <- 3}
if(!sim_covar){tran_covar_num <- 6}


Tran2DF <- function(params_tran){
  race_list <- c("0","1","2","3","4","5") 
  
  
  tran_vec <- sapply(c(1:96),FUN = Params2Tran,params_tran = params_tran,index=1)
  tran_df <- data.frame(prob = c(tran_vec[3,],tran_vec[2,]),
                        type = rep(c("Sleep", "Wake"),each= 96),
                        time = rep(c(1:96)/4,2),
                        race = race_list[1])
  
  for (i in 2:tran_covar_num){
    tran_vec <- sapply(c(1:96),FUN = Params2Tran,params_tran = params_tran,index=i)
    tran_df_working <- data.frame(prob = c(tran_vec[3,],tran_vec[2,]),
                                  type = rep(c("Sleep", "Wake"),each= 96),
                                  time = rep(c(1:96)/4,2),
                                  race = race_list[i])
    
    tran_df <- rbind(tran_df,tran_df_working)
  }
  
  return(tran_df)
}

readCpp <- function(path) {
  tryCatch(
    {
      sourceCpp(file = path)
    },
    error = function(cond) {
      message("Wrong environment")
      # Choose a return value in case of error
      NA
    },
    warning = function(cond) {
      message("Wrong environment")
      # Choose a return value in case of warning
      NULL
    },
    finally = {
      message("Done")
    }
  )
}

CalcMeansWake <- function(act,weights_array){
  n <- dim(weights_array)[2]
  re_num <- dim(weights_array)[3]
  mean_vec <- numeric(re_num)
  for (re_index in 1:re_num){
    num <- 0
    denom <- 0
    for (ind in 1:n){
      act_ind0 <- act[,ind] 
      inds_keep <- !is.na(act_ind0)
      
      num <- num + (weights_array[inds_keep,ind,re_index]) %*% act_ind0[inds_keep]
      denom <- denom + sum(weights_array[inds_keep,ind,re_index])
    }
    mean_vec[re_index] <- num/denom
  }
  return(mean_vec)
}

CalcMeansWakeInd <- function(act,weights_array){
  n <- dim(weights_array)[2]
  mean_vec <- numeric(n)
  for (ind in 1:n){
    act_ind0 <- act[,ind] 
    inds_keep <- !is.na(act_ind0)
    
    num <- (weights_array[inds_keep,ind,1]) %*% act_ind0[inds_keep]
    denom <- sum(weights_array[inds_keep,ind,1])
    
    mean_vec[ind] <- num/denom
    if(is.na(mean_vec[ind])){mean_vec[ind] <- 2}
  }
  
  return(mean_vec)
}

CalcSigmaWake <- function(act,weights_array,mean_vec){
  n <- dim(weights_array)[2]
  re_num <- dim(weights_array)[3]
  num <- 0
  denom <- 0
  
  for (re_index in 1:re_num){
    for (ind in 1:n){
      
      resid0 <- (act[,ind]-mean_vec[re_index])
      inds_keep <- !is.na(resid0)
      
      num <- num + (weights_array[inds_keep,ind,re_index]) %*% (resid0[inds_keep]^2)
      denom <- denom + sum(weights_array[inds_keep,ind,re_index])
    }
  }
  return(sqrt(num/(denom)))
}

CalcSigmaWakeInd <- function(act,weights_array,wake_means_ind){
  n <- dim(weights_array)[2]
  sigma_vec <- numeric(n)
  for (ind in 1:n){
    resid0 <- (act[,ind]-wake_means_ind[ind])
    inds_keep <- !is.na(resid0)
    
    num <- (weights_array[inds_keep,ind,1]) %*% (resid0[inds_keep]^2)
    denom <- sum(weights_array[inds_keep,ind,1])
    
    # if (num == 0){num <- .0001}
    
    sigma_vec[ind] <- sqrt(num/denom)
    
    if(sigma_vec[ind] == Inf | is.na(sigma_vec[ind]) | sigma_vec[ind] == 0){
      print("AAAAAAAAA")
      sigma_vec[ind] <- .01
    }
  }
  
  return(sigma_vec)
}

CalcMeanSleep <- function(act,weights_mat,lepsilon,pop_pool = T){
  n <- dim(weights_mat)[2]
  num <- 0
  denom <- 0
  mean_vec <- numeric(n)
  
  leps_indicator <- (act==lepsilon)
  leps_indicator[is.na(leps_indicator)] <- F
  for (ind in 1:n){
    
    act_ind0 <- act[,ind] 
    inds_keep <- !is.na(act_ind0)
    
    leps_ind <- leps_indicator[,ind]
    inds_keep <- (inds_keep & !leps_ind)
    
    num <- num + (1-weights_mat[inds_keep,ind]) %*% act_ind0[inds_keep]
    denom <- denom + sum(1-weights_mat[inds_keep,ind])
    mean_vec[ind] <- ((1-weights_mat[inds_keep,ind]) %*% act_ind0[inds_keep])/sum(1-weights_mat[inds_keep,ind])
    
    if(is.na(mean_vec[ind])){mean_vec[ind] <- -1}
  }
  if (!pop_pool){return(mean_vec)}
  return(num/denom)
}

CalcSigmaSleep <- function(act,weights_mat,sleep_mean, lepsilon,pop_pool=T){
  n <- dim(weights_mat)[2]
  num <- 0
  denom <- 0
  sigma_vec <- numeric(n)
  
  leps_indicator <- (act==lepsilon)
  leps_indicator[is.na(leps_indicator)] <- F
  
  for (ind in 1:n){
    
    if (pop_pool){resid0 <- c(act[,ind] - sleep_mean)
    } else{resid0 <- c(act[,ind] - sleep_mean[ind])}
    
    inds_keep <- !is.na(resid0)
    
    leps_ind <- leps_indicator[,ind]
    inds_keep <- (inds_keep & !leps_ind)
    
    num <- num + (1-weights_mat[inds_keep,ind]) %*% resid0[inds_keep]^2
    denom <- denom + sum(1-weights_mat[inds_keep,ind])
    sigma_vec[ind] <- sqrt(((1-weights_mat[inds_keep,ind]) %*% resid0[inds_keep]^2)/sum(1-weights_mat[inds_keep,ind]))
    
    if(sigma_vec[ind] == Inf | is.na(sigma_vec[ind]) | sigma_vec[ind] == 0){
      print("BBBBBBBB")
      sigma_vec[ind] <- .01
    }
  }
  if (!pop_pool){return(sigma_vec)}
  return(sqrt(num/denom))
}

TranIndList <- function(params_tran_mat,obs_per_day,num_of_people){
  tran_ind_list <- list()
  
  for (ind in 1:num_of_people){
    tran_ind_list[[ind]] <- TranByTimeVec(1,params_tran_mat[ind,],c(1:obs_per_day))
  }
  
  return(tran_ind_list)
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
library(mclust)
library(Rcpp)
library(RcppArmadillo)
library(profvis)
# library(emdbook)

readCpp("Scripting/cFunctions.cpp")
readCpp("/panfs/jay/groups/29/mfiecas/aron0064/ActHMM/Rcode/cFunctions.cpp")

#### Set True Parameters ####
central_mean <- wake_params[1]


var_factor <- RE_num %/% 2 * runif(1,.8,1.2)
re_set <- seq(-2*var_factor,2*var_factor,length.out = RE_num) * runif(RE_num,.65,.85)

pi_l_true <- rep(1/RE_num,RE_num)
if(RE_num == 0){pi_l_true <- c(1)}

if (RE_num <= 1){
  re_set <- c(0)
} 



mean_set_true <- central_mean + re_set


init_true <- c(0.3347129,0.6652871)

if (!sim_covar){
  params_tran_true <- c(-3.23,.05,.1,.07,.09,.09,1.17,.5,
                        -2.69,.04,.08,.13,.15,.21,-1.26,-.58)
} else{
  params_tran_true <- c(-3,.5,1,0,0,0,.9,-.8,
                        -2.2,.3,.6,0,0,0,-.75,.8)
}




emit_act_true <- array(dim = c(2,2,length(mean_set_true)))
emit_act_true[1,1,] <- mean_set_true
emit_act_true[1,2,] <- wake_params[2]
emit_act_true[2,1,] <- sleep_params[1]
emit_act_true[2,2,] <- sleep_params[2]


colnames(emit_act_true) <- c("Mean","Std Dev")
rownames(emit_act_true) <- c("Wake","Sleep")

emit_act_true_sim <- emit_act_true[,,1]
emit_act_true_sim[1,1] <- central_mean

#Prob of being below detection limit
act_light_binom_true <- c(.1)
# act_light_binom_true <- c(0)

#### Simulate True Data ####

if (!real_data){
  if (is.na(sim_num)){sim_num <- 99}
  if (set_seed){set.seed(sim_num)}
  
  if (sim_size == 0){
    day_length <- 96 
    num_of_people <- 1000
  } else if (sim_size == 1){
    day_length <- 96 
    num_of_people <- 1000
  } else if (sim_size == 2){
    day_length <- 96
    num_of_people <- 1000 * 5
  } else if (sim_size == 3){
    day_length <- 96 * 7  
    num_of_people <- 1000
  } else if (sim_size == 4){
    day_length <- 96 * 7  
    num_of_people <- 1000 * 5
  }
  
  log_sweights_vec <- numeric(num_of_people)
  
  
  n <- day_length * num_of_people
  id <- data.frame(SEQN = c(1:num_of_people))
  
  simulated_hmm <- SimulateHMM(day_length,num_of_people,init_true,params_tran_true,emit_act_true_sim,
                               act_light_binom_true, pi_l_true,re_set)
  
  mc <- simulated_hmm[[1]]
  act <- simulated_hmm[[2]]
  covar_mat_tran <- simulated_hmm[[3]]
  tran_ind_vec <- apply(covar_mat_tran,1,ChooseTran)
  act_lod <- simulated_hmm[[4]]
  clust_ind_true <- simulated_hmm[[5]]
  
  trans_true_emp <- CalcEmpiricTran(mc,covar_mat_tran)
  init_true_emp <- c(sum(mc[1,] == 0),sum(mc[1,] == 1)) / dim(mc)[2]
  
  # emit_act_true_emp <- CalcEmpiricAct(mc,act,act_lod,clust_ind_true)
  ###NEED TO THINK ABOUT THIS MORE
  emit_act_true_emp <- emit_act_true
  
  act_light_binom_true_emp <- c(sum(act_lod[mc==1]==1)/length(act_lod[mc==1]))
  
  pi_l_true_emp <- table(clust_ind_true)/ length(clust_ind_true)
  
  
  
  #CHECK MISSING AND CORRELATION BEFORE SIM
  act <- apply(act,2,InduceMissingVec, prob = .3)
  # act <- apply(act,2,InduceMissingVec, prob = 0)
}

#### Initialize starting parameters ####
init  <- c(runif(1,.1,.5),0)
init[2] <- 1 - init[1]


params_tran <- params_tran_true + runif(16,-.2,.2)
if (sim_covar){params_tran[c(4:6,12:14)] <- 0} 


emit_act <- emit_act_true
emit_act[1,1,]  <- emit_act[1,1,] + runif(length(emit_act[1,1,]),-1/2,1/2)
emit_act[1,2,] <- emit_act[1,2,] + runif(1,-1/4,1/4)
emit_act[2,1,] <- emit_act[2,1,] + runif(1,-1/4,1/4)
emit_act[2,2,] <- emit_act[2,2,] + runif(1,-1/2,1/2)

emit_act_array <- array(emit_act,dim=c(2,2,num_of_people))

# act_light_binom <- c(runif(1,0,.1))
act_light_binom <- .1
act_light_binom_ind <- rep(.1, num_of_people)

pi_l <- pi_l_true + runif(length(pi_l_true),0,.25)
pi_l <- pi_l/sum(pi_l)

#### Load in Real Data ####
if(real_data){
  
  load("Wavedata_G.rda")
  load("Wavedata_H.rda")
  
  sim_num <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
  
  
  id_G <- wave_data_G[[4]]
  id_H <- wave_data_H[[4]]
  id <- rbind(id_G,id_H)
  
  log_sweights_vec_G <- log(id_G[,3])
  log_sweights_vec_H <- log(id_H[,3])
  
  log_sweights_vec <- c(log_sweights_vec_G/2,log_sweights_vec_H/2)
  
  #CHANGE RACE AND POV
  covar_mat_tran <- matrix(0,ncol = 6, nrow = length(id$poverty))
  covar_mat_tran[,1] <- 1
  for (i in 1:length(id$poverty)){
    covar_mat_tran[i,id$poverty[i]] <- 1
  } 
  
  tran_ind_vec <- apply(covar_mat_tran,1,ChooseTran)
  
  act_G <- log(wave_data_G[[1]] + epsilon)
  light_G <- log(wave_data_G[[2]] + epsilon)
  act_H <- log(wave_data_H[[1]] + epsilon)
  light_H <- log(wave_data_H[[2]] + epsilon)
  
  #remove SEQN identifier
  act_G <- t(act_G)[2:865,]
  light_G <- t(light_G)[2:865,]
  act_H <- t(act_H)[2:865,]
  light_H <- t(light_H)[2:865,]
  
  act <- cbind(act_G,act_H)
  
  day_length <- dim(act)[1]
  num_of_people <- dim(act)[2]
  
  print("Loaded Data")
  
  init  <- c(1/3,2/3)
  
  params_tran <- c(-3,.01,.01,.01,.01,.01,.01,.01,-2.2,.01,.01,.01,.01,.01,.01,.01)
  
  #need to fix this
  # emit_act <- matrix(c(4,2,
  #                      0,2), byrow = T, 2)
  
  act_light_binom <- c(.05)
}


#### EM #### 


# init <- init_true_emp
# params_tran <- params_tran_true
# emit_act <- emit_act_true_emp
# pi_l <- pi_l_true_emp

init <- init_true
params_tran <- params_tran_true
emit_act <- emit_act_true
pi_l <- pi_l_true
emit_act_array <- array(matrix(c(2.5,2,-1,2),2,2,byrow = T),dim=c(2,2,num_of_people))

tran_list <- lapply(c(1:dim(covar_mat_tran)[2]),TranByTimeVec, params_tran = params_tran, time_vec = c(1:obs_per_day))
# tran_ind_list <- TranIndList(params_tran_mat,obs_per_day,num_of_people)

time_vec <- c()

print("PRE ALPHA")
if (RE_num==0){
  alpha <- ForwardIndAll(act,init,tran_list,emit_act_array,tran_ind_vec,lepsilon, log_sweights_vec)
  beta <- BackwardIndAll(act,tran_list,emit_act_array,tran_ind_vec,lepsilon)
} else {
  alpha <- ForwardC(act,init,tran_list,emit_act,tran_ind_vec,lepsilon, log_sweights_vec)
  beta <- BackwardC(act,tran_list,emit_act,tran_ind_vec,lepsilon)
}

# apply(alpha[[2]][,,1]+beta[[2]][,,1],1,logSumExp)

new_likelihood <- CalcLikelihood(alpha,pi_l)
likelihood_vec <- c(new_likelihood)
likelihood <- -Inf
like_diff <- new_likelihood - likelihood

# grad_num <- grad(LogLike,params_tran)

while(abs(like_diff) > 1e-3){
  # for(alsdkjf in 1:1){ 
  
  tic()
  start_time <- Sys.time()
  likelihood <- new_likelihood
  
  
  ################## MC Parameters
  
  
  init <- CalcInit(alpha,beta,pi_l,T)
  
  if (RE_num == 0){
    params_tran <- CalcTranBothC(alpha,beta,act,params_tran,emit_act_array,covar_mat_tran,pi_l,F)
  } else {
    params_tran <- CalcTranBothC(alpha,beta,act,params_tran,emit_act,covar_mat_tran,pi_l,T)
    
  }
  
  tran_list <- lapply(c(1:dim(covar_mat_tran)[2]),TranByTimeVec, params_tran = params_tran, time_vec = c(1:obs_per_day))
  # tran_ind_list <- TranIndList(params_tran_mat,obs_per_day,num_of_people)
  
  ################## Weights
  #Weights are prob currently in the wake state
  weights_array <- exp(CondMarginalize(alpha,beta,pi_l))
  weights_mat <-rowSums(weights_array, dims = 2)
  weights_vec <- as.vector(weights_mat)
  
  
  ################## Binom param (need to implement)
  # act_vec <- as.vector(act)
  # lod_act_weight <- as.numeric(act_vec==log(epsilon))
  # act_light_binom[1] <- sum(lod_act_weight,na.rm = T)/sum(1-weights_vec[!is.na(as.vector(act))])
  # act_light_binom_ind <- colSums(act==lepsilon,na.rm = T)/colSums(1-weights_mat*!is.na(act),na.rm = T)
  # act_light_binom_ind[act_light_binom_ind==0] <- .00001
  # act_light_binom_ind[act_light_binom_ind>=1] <- .99999
  
  
  
  ################## Emission Dist Param
  
  
  re_prob <- CalcProbRE(alpha,pi_l)
  pi_l <- colSums(re_prob)/num_of_people
  
  pi_l[pi_l<1e-100] <- 1e-100
  
  if (RE_num <= 1){pi_l <- c(1)}
  
  if (RE_num != 0){
    wake_means <- CalcMeansWake(act,weights_array)
    wake_sigma <- CalcSigmaWake(act,weights_array,wake_means)
    sleep_mean <- CalcMeanSleep(act,weights_mat,lepsilon)[[1]]
    sleep_sigma <- CalcSigmaSleep(act,weights_mat,sleep_mean,lepsilon,T)
    
    emit_act[1,1,] <- wake_means
    emit_act[1,2,] <- wake_sigma
    emit_act[2,1,] <- sleep_mean
    emit_act[2,2,] <- sleep_sigma
    emit_act[,2,] <- abs(emit_act[,2,])
    
    if (any(is.na(wake_means))){
      print("NA Wake Means")
      break
    }
    
  } else {
    wake_means_ind <- CalcMeansWakeInd(act,weights_array)
    wake_sigma_ind <- CalcSigmaWakeInd(act,weights_array,wake_means_ind)
    sleep_mean_ind <- CalcMeanSleep(act,weights_mat,lepsilon,F)
    sleep_sigma_ind <- CalcSigmaSleep(act,weights_mat,sleep_mean_ind,lepsilon,F)
    
    emit_act_array[1,1,] <- wake_means_ind
    emit_act_array[1,2,] <- wake_sigma_ind
    emit_act_array[2,1,] <- sleep_mean_ind
    emit_act_array[2,2,] <- sleep_sigma_ind
    # emit_act_array[,2,] <- abs(emit_act[,2,])
  }
  
  ##################
  
  #Reorder to avoid label switching
  #Cluster means go from small to large
  if (length(pi_l)>1){
    reord_inds <- order(emit_act[1,1,])
    emit_act <- emit_act[,,reord_inds]
    pi_l <- pi_l[reord_inds]
  }
  ##################
  
  
  if (RE_num==0){
    alpha <- ForwardIndAll(act,init,tran_list,emit_act_array,tran_ind_vec,lepsilon, log_sweights_vec)
    beta <- BackwardIndAll(act,tran_list,emit_act_array,tran_ind_vec,lepsilon)
  } else {
    alpha <- ForwardC(act,init,tran_list,emit_act,tran_ind_vec,lepsilon, log_sweights_vec)
    beta <- BackwardC(act,tran_list,emit_act,tran_ind_vec,lepsilon)
  }
  
  new_likelihood <- CalcLikelihood(alpha,pi_l)
  like_diff <- new_likelihood - likelihood
  print(paste("RE num:",RE_num,"Like:",round(like_diff,6)))
  likelihood_vec <- c(likelihood_vec,new_likelihood)
  
  end_time <- Sys.time()
  time_vec <- c(time_vec,as.numeric(difftime(end_time, start_time, units = "secs")))
  # print(paste("RE num",RE_num, "memory",sum(gc(T)[,6])))
  toc()
}

if (RE_num == 0){emit_act <- array(apply(emit_act_array,c(1,2),median),dim = c(2,2,1))}

print(paste("Sim Num:",sim_num,"RE Num:",RE_num,"Ending"))

decoded_mat <- sapply(c(1:num_of_people), ViterbiInd, RE_num=RE_num)

starting_conditions <- list(wake_params,
                            sleep_params,
                            c(num_of_people,day_length),
                            RE_type,
                            RE_num,
                            sim_num,
                            time_vec)

if (!real_data){
  total_acc <- (sum(decoded_mat == 0 & mc ==0) + sum(decoded_mat == 1 & mc ==1))/(num_of_people*day_length)
  wake_acc <- sum(decoded_mat == 0 & mc ==0)/sum(mc ==0)
  sleep_acc <- sum(decoded_mat == 1 & mc ==1)/sum(mc ==1)
  
  perc_pred_wake <- sum(decoded_mat == 0 & mc ==0)/sum(decoded_mat ==0)
  perc_pred_sleep <- sum(decoded_mat == 1 & mc ==1)/sum(decoded_mat == 1) 
  
  
  
  tran_df <- Tran2DF(params_tran) %>% 
    mutate(truth = Tran2DF(params_tran_true)[,1]) %>% 
    mutate(resid = prob - truth)
  
  Q1 <- function(x){return(quantile(x,probs = c(.01)))}
  Q2.5 <- function(x){return(quantile(x,probs = c(.025)))}
  Q10 <- function(x){return(quantile(x,probs = c(.1)))}
  Q25 <- function(x){return(quantile(x,probs = c(.25)))}
  
  
  Q75 <- function(x){return(quantile(x,probs = c(.75)))}
  Q90 <- function(x){return(quantile(x,probs = c(.9)))}
  Q975 <- function(x){return(quantile(x,probs = c(.975)))}
  Q99 <- function(x){return(quantile(x,probs = c(.99)))}
  
  tran_sum_df <- tran_df %>% group_by(type,race) %>% summarise(across(resid, 
                                                                      list(min = min, Q1 = Q1, Q10 = Q10, Q25 = Q25, med = median,
                                                                           Q75 = Q75, Q90 = Q90, Q99 = Q99, max = max)))
  
  tran_list <- list(tran_df,tran_sum_df)
}

if (!real_data){
  true_params <- list(init_true_emp,params_tran_true,emit_act_true_emp,act_light_binom_true_emp,pi_l_true_emp,clust_ind_true)
  est_params <- list(init,params_tran,emit_act,act_light_binom_true_emp,pi_l,re_prob)
  mc_list <- list(mc,decoded_mat,total_acc,wake_acc,sleep_acc,perc_pred_wake,perc_pred_sleep)
  params_to_save <- list(true_params,est_params,likelihood_vec,mc_list,starting_conditions,tran_list)
  
  if(sim_size != 0){
    save(params_to_save,file = paste0(RE_type,RE_num,"Size",sim_size,"Seed",sim_num,".rda"))
  }
  
} else {
  est_params <- list(init,params_tran,emit_act,.05,pi_l,re_prob)
  params_to_save <- list(est_params,likelihood_vec,decoded_mat,starting_conditions,Tran2DF(params_tran))
  
  save(params_to_save,file = paste0("MHMM",RE_num,".rda"))
}


#####################

# readCpp("Scripting/cFunctions.cpp")
# logClassificationC(current_state = 1,act_obs = act[,1],
#                    mu = emit_act[2,1,1],sig = emit_act[2,2,1],
#                    act_binom = act_light_binom,lod = lepsilon)
# 
# ForwardIndC(act[,190],init,tran_list,emit_act,1,act_light_binom,0,lepsilon)

# for (i in 1:length(alpha)){x <- ViterbiInd(i,RE_num)}
