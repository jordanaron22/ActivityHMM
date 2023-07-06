library(foreign)
library(tidyverse)
library(data.table)


varName_wave = Sys.getenv('varNameWave')

var_name <- strsplit(varName_wave," ")[[1]][1]
wave <- strsplit(varName_wave," ")[[1]][2]
print(varName_wave)
print(var_name)
print(wave)

if (wave == "H"){
  acth_full <- read.xport("../data/PAXMIN_H.XPT")
} else if (wave == "G"){
  acth_full <- read.xport("../data/PAXMIN_G.XPT")
}

gc()

acth <- acth_full %>% mutate(minind = (PAXSSNMP %/% 4800) + 1)
acth <- acth %>% subset(minind != 11530)

seqn_list <- unique(acth$SEQN)

rm(acth_full)
gc()

NHANESWideToLong <- function(acth,var_name){
  seqn_list <- unique(acth$SEQN)
  
  working_wide <- NHANESWideToLongHelper(acth,var_name,seqn_list[1])
  
  for (i in 2:length(seqn_list)){
    working_wide_ind <- NHANESWideToLongHelper(acth,var_name,seqn_list[i])
    
    working_wide <- rbind(working_wide,working_wide_ind)
    
    if (i == 10){print(i)}
    if (i %% 200 == 0){print(i)}
  }
  return(working_wide)
}

NHANESWideToLongHelper <- function(acth, var_name, seqn){
  # act_sub <- acth %>% select(SEQN, PAXDAYM, PAXDAYWM, all_of(var_name),minind)
  
  working_act_sub <- acth %>% select(SEQN, PAXDAYM, PAXDAYWM, all_of(var_name),minind) %>% filter(SEQN == seqn)
  days_rec <- unique(working_act_sub$PAXDAYM)
  
  first_day_time <- sum(working_act_sub$PAXDAYM == 1)
  last_day_time <- sum(working_act_sub$PAXDAYM == days_rec[length(days_rec)])
  
  if (length(days_rec) > 1){
    min_index <- c(seq(from = (1440 - first_day_time + 1), to = 1440),
                   rep(seq(1,1440),(length(days_rec)-2)),
                   seq(1,last_day_time))
  } else if (length(days_rec) == 1) {
    min_index <- c(seq(from = (1440 - first_day_time + 1), to = 1440))
  }
  
  working_act_sub <- working_act_sub %>% mutate(min_index = min_index)
  
  long_act_ind <- pivot_wider(working_act_sub, id_cols = c(SEQN,PAXDAYM,PAXDAYWM), names_from = c(min_index), values_from = all_of(var_name))
  
  if (length(days_rec) == 2){
    cols_to_add <- seq(last_day_time + 1, 1440 - first_day_time)
    df_to_add <- data.frame(V1 = rep(NA,length(cols_to_add)),
                            V2 = rep(NA,length(cols_to_add)))
    rownames(df_to_add) <- cols_to_add
    long_act_ind <- cbind(long_act_ind,t(df_to_add))
    
  } else if (length(days_rec) == 1){
    cols_to_add <- seq(1, 1440 - first_day_time)
    df_to_add <- data.frame(V1 = rep(NA,length(cols_to_add)))
    rownames(df_to_add) <- cols_to_add
    long_act_ind <- cbind(long_act_ind,t(df_to_add))
  }
  
  long_act_ind <- long_act_ind[,c(colnames(long_act_ind)[1:3],paste0(seq(1,1440)))]
  
  
  return(long_act_ind)
}

## Light - PAXLXSDM
## Activity - PAXMTSM
## Wake/sleep/non-wear pred - PAXPREDM
## Quality Flag - PAXQFM 
## Seconds - PAXTSM
## Sleep mode - PAXAISMM
long_data <- NHANESWideToLong(acth,var_name)

save(long_data, file = paste0(var_name,wave,"_long.rda"))



