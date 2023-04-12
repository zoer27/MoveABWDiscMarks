####Simulation Testing for estimation power of ABW movement models
#12/9/22 
#Zoe Rand

library(tidyverse)
library(popbio)
library(cmdstanr)
library(MASS)
library(posterior)
library(bayesplot)
library(ggridges)
options(mc.cores = parallel::detectCores())

# Data ----------------------------------------------------------
Abund_Dat<-read_csv("Custom_Model/Data/basin_abundances_withMat.csv") #for dates and basins of abundance estimates
AtYears<-Abund_Dat %>% filter(Basin == "Atlantic") %>% dplyr::select(Year) %>% as.matrix()
IndYears<-Abund_Dat %>% filter(Basin == "Indian") %>% dplyr::select(Year) %>% as.matrix()
PacYears<-Abund_Dat %>% filter(Basin == "Pacific") %>% dplyr::select(Year) %>% as.matrix()
#last year is 2008

#release data
Rel_dat<-read_csv("Custom_Model/Data/mark_releases.csv")

#catch data
Catch_dat<-read_csv("Custom_Model/Data/catchcorrected.csv")

# Simulation Function ---------------------------------------------

sim_recs_nbs<-function(K, m_mat, tl, catch_dat, release_dat, r, s, cv_dat, theta, nyearPop, seed = 123){
  #This function takes the following parameters:
  #K = SO-wide carrying capacity/initial population size
  #m_mat = movement matrix between basins (3x3)
  #tl = mark loss paramter
  #catch_dat = matrix of catches in each basin and year (3xNyear)
  #release_dat = mark release data in matrix form
  #r = intrinsic growth
  #s = natural survival
  #cv_dat = vector of cv's for each year of abundnace estimates
  #theta = overdispersion parameter
  #nyearPop = number of years to run the population model
  #seed = random number seed to use when simulating data
  
  #This function returns as a list: 
  #Abund = simulated abundance data
  #Recs = simulated recovery data
  #Npop = predicted population in each year and basin from population model (for checking sims)
  #stat_dist = stationary distribution from the movement matrix (for checking)
  
  #setting seed
  set.seed(seed)
  
  #setting up parameters
  rvec<-eigen(m_mat)
  lvec<-ginv(rvec$vectors)
  stat_dist<-lvec[1,]/sum(lvec[1,])
  K_vect<-stat_dist*K
  MSY<-(r*K_vect*2.39)/((3.39)^((1/2.39) + 1))
  r_star <- r/(1-s)
  
  #projecting population
  
  Npop<-matrix(data = NA, nrow = nyearPop, ncol = 3)
  if(nrow(catch_dat < nyearPop)){
    catch0<-matrix(0, ncol = 3, nrow = nyearPop-nrow(catch_dat))
    catch_dat<-rbind(catch_dat, catch0)
  }
  
  Npop[1,] = K_vect #first year is carrying capacity
  
  for (i in 2:nyearPop){
    for(I in 1:3){
      Npop[i,I] = s*Npop[i-1,I] + (1-s)*Npop[i-1,I]*(1+r_star*(1-(((1-s)*r_star*2.39*Npop[i-1,I])/(MSY[I]*(3.39^(1/2.39 +1))))^2.39))
      
      if(catch_dat[i-1,I] <= (Npop[i,I])){ #limiting harvest rate so can't be greater than max set
        Npop[i,I] = Npop[i,I] - catch_dat[i-1,I]
      }else{
        Npop[i,I] = 1; #if catch is greater than abundance turns pop into 1
      }
      Npop[i,I] = max(Npop[i,I], 1)#if population goes below 1 this becomes 1
    }
    Npop[i,] = Npop[i,] %*% m_mat
    #print(Npop[i,])
  }
  
  #harvest rate calculation
  h<-matrix(data = NA, nrow = nrow(release_dat), ncol = 3)
  for (i in 1:(nrow(release_dat))){ #harvest starts in 1926 and ends in 1972
    for (I in 1:3){
      h[i,I] = catch_dat[i+13,I]/Npop[i+13,I]
      h[i,I] = min(h[i,I], 1) #harvest rate can't be more than 1 (if pop is 1 then all whales are harvested and harvest rate is 1)
    }
  }
  
  PredH = h[2:nrow(release_dat),] #no harvest rates before 1927
  #mark model
  NT<-matrix(data = NA, nrow = nrow(release_dat), ncol = 9)
  NTR<-NT[-1,] #recovered tags don't include first year
  
  NT[1,1] = release_dat[1,1] #first year releases
  NT[1, 2:4] = rep(0,3)
  NT[1,5] = release_dat[1,2]
  NT[1,6:8] = rep(0,3)
  NT[1,9] = release_dat[1,3]
  m<-m_mat #so don't have to rewrite in the model
  
  for(i in 2:(nrow(release_dat))){ #recoveries for all years
    new_tags = release_dat[i,]
    NT[i,1] = (NT[i-1, 1]*(1-h[i-1,1])*s*m[1,1]) + (NT[i-1,2]*(1-h[i-1,2])*s*m[2,1]) + (NT[i-1,3]*(1-h[i-1,3])*s*m[3,1])
    NT[i,2] = (NT[i-1, 1]*(1-h[i-1,1])*s*m[1,2]) + (NT[i-1,2]*(1-h[i-1,2])*s*m[2,2]) + (NT[i-1,3]*(1-h[i-1,3])*s*m[3,2])
    NT[i,3] = (NT[i-1, 1]*(1-h[i-1,1])*s*m[1,3]) + (NT[i-1,2]*(1-h[i-1,2])*s*m[2,3]) + (NT[i-1,3]*(1-h[i-1,3])*s*m[3,3])
    NT[i,4] = (NT[i-1, 4]*(1-h[i-1,1])*s*m[1,1]) + (NT[i-1,5]*(1-h[i-1,2])*s*m[2,1]) + (NT[i-1,6]*(1-h[i-1,3])*s*m[3,1])
    NT[i,5] = (NT[i-1, 4]*(1-h[i-1,1])*s*m[1,2]) + (NT[i-1,5]*(1-h[i-1,2])*s*m[2,2]) + (NT[i-1,6]*(1-h[i-1,3])*s*m[3,2])
    NT[i,6] = (NT[i-1, 4]*(1-h[i-1,1])*s*m[1,3]) + (NT[i-1,5]*(1-h[i-1,2])*s*m[2,3]) + (NT[i-1,6]*(1-h[i-1,3])*s*m[3,3])
    NT[i,7] = (NT[i-1, 7]*(1-h[i-1,1])*s*m[1,1]) + (NT[i-1,8]*(1-h[i-1,2])*s*m[2,1]) + (NT[i-1,9]*(1-h[i-1,3])*s*m[3,1])
    NT[i,8] = (NT[i-1, 7]*(1-h[i-1,1])*s*m[1,2]) + (NT[i-1,8]*(1-h[i-1,2])*s*m[2,2]) + (NT[i-1,9]*(1-h[i-1,3])*s*m[3,2])
    NT[i,9] = (NT[i-1, 7]*(1-h[i-1,1])*s*m[1,3]) + (NT[i-1,8]*(1-h[i-1,2])*s*m[2,3]) + (NT[i-1,9]*(1-h[i-1,3])*s*m[3,3])
    for(I in 1:3){
      NTR[i-1,I] = NT[i,I] * PredH[i-1, I];
      NTR[i-1, I+3] = NT[i,I + 3] * PredH[i-1, I];
      NTR[i-1, I + 6] = NT[i,I + 3] * PredH[i-1, I];
    }
    NT[i,1] = NT[i, 1] + new_tags[1];
    NT[i,5] = NT[i, 5] + new_tags[2];
    NT[i, 9] = NT[i,9] + new_tags[3];
  }
  
  
  NTR = NTR*(1-tl) #recovered tags are not lost
  
  for(i in 1:(nrow(release_dat)-1)){
    for(k in 1:9){
      NTR[i,k] = max(NTR[i,k], 0.001) #if 0 is in predictions replace with 0.001 so likelihood works
    }
  }
  
  #simulating population data 
  
  Abund<-tibble("Years" = c(AtYears, IndYears, PacYears), 
                "Basin" = c(rep("Atlantic", nrow(AtYears)), rep("Indian", nrow(IndYears)), rep("Pacific", nrow(PacYears))),
                "Npop" = NA)

  CV<-cv_dat
  sig<-sqrt(log(CV^2 +1))
  
  Years<-seq(min(Catch_dat$Year), min(Catch_dat$Year) + nyearPop, by = 1)
  At_index<-which(Years %in% AtYears)
  In_index<-which(Years %in% IndYears)
  Pac_index<-which(Years %in% PacYears)
  
  
  Abund$Npop[Abund$Basin == "Atlantic"]<-rlnorm(nrow(AtYears), meanlog = log(Npop[At_index,1]), sdlog = sig)
  Abund$Npop[Abund$Basin == "Indian"]<-rlnorm(nrow(IndYears), meanlog = log(Npop[In_index,2]), sdlog = sig)
  Abund$Npop[Abund$Basin == "Pacific"]<-rlnorm(nrow(PacYears), meanlog = log(Npop[Pac_index,3]), sdlog = sig)
  

  #simulating recovery data
  
  #negative binomial simulation--overdispersed data
  Recs<-rnbinom(length(NTR), size = theta, mu = NTR)
  
  dim(Recs)<-dim(NTR)
  
  return(list(Abund = Abund, Recs = Recs, Npop = Npop, stat_dist = stat_dist))
}

# Functions To Run Simulations Multiple Times -----------------------------

#compiling model
file<-"Code/Stan Files/SimulationModel.stan"
mod<-cmdstan_model(file)

#MSY fun for initial values
MSYfun<-function(K, r){
  s = 0.96
  r_star = r/(1-s)
  MSY = ((1-s)*r_star*K*2.39)/((3.39)^((1/2.39) + 1))
  return(MSY)
}

#tells you whether the "true" value of the parameter is included in the 95% credible interval
true_val<-function(K, m, tl, r, Kbsn, theta, sumstats){
  included<-rep(NA, 12)
  if(K >= exp(sumstats$q025[1]) & K <= exp(sumstats$q975[1])){ 
    included[1]<-1
  } else{
    included[1]<-0
  }
  for(i in 1:9){
    if(as.vector(m)[i] >= sumstats$q025[i+1] & as.vector(m)[i] <= sumstats$q975[i+1]){
      included[i+1] <-1
    } else{
      included[i+1]<-0
    }
  }
  if(tl >= sumstats$q025[11] & tl <= sumstats$q975[11]){
    included[11]<-1
  } else{
    included[11]<-0
  }
  if(r >= sumstats$q025[12] & r <= sumstats$q975[12]){
    included[12]<-1
  } else{
    included[12]<-0
  }
  if(Kbsn[1] >= sumstats$q025[13] & Kbsn[1] <= sumstats$q975[13]){
    included[13]<-1
  } else{
    included[13]<-0
  }
  if(Kbsn[2] >= sumstats$q025[14] & Kbsn[2] <= sumstats$q975[14]){
    included[14]<-1
  } else{
    included[14]<-0
  }
  if(Kbsn[3] >= sumstats$q025[15] & Kbsn[3] <= sumstats$q975[15]){
    included[15]<-1
  } else{
    included[15]<-0
  }
  if(theta >= sumstats$q025[16] & theta <= sumstats$q975[16]){
    included[16]<-1
  } else{
    included[16]<-0
  }
  return(included)
}

simfunc<-function(nint, K, m, tl, r, s, cv_dat, theta, catch_dat, release_dat){ 
  #This function takes a number of iterations, a set of "true" parameter values, 
    #and the catch and mark release data and simulates abundance and mark recovery datasets
    #Then it fits the model using these simulated datasets
    #And calculates wether the "true" value of the parameter is included in the credible
    #intervals returned by the model
  
  #vector to store the seeds that were used
  seeds<-rep(NA, nint)
  
  
  #Information for Stan
  N_basin<-3
  Basin_index<-c(1,2,3)
  nyeartag<-nrow(release_dat)
  N_yearcatch<-nrow(catch_dat)
  
  #initial values
  init_fun <- function() {
    rinit = runif(1, 0.06, 0.08)
    out<-list(
      MAB=runif(1, 0.2, 0.33), MAC=runif(1, 0.01, 0.33), MBA=runif(1, 0.2, 0.33),
      MBC=runif(1, 0.2, 0.33), MCA=runif(1, 0.05, 0.33), MCB=runif(1, 0.2, 0.33),
      tl = runif(1,0.93,1),
      r = rinit,
      lnMSY = log(MSYfun(250000, rinit)),
      inv_theta = 1/runif(1, 5, 10))
    return(out)
  }
  
  #to save results
  results<-tibble(par = c("lnK", "m[1,1]","m[2,1]","m[3,1]","m[1,2]","m[2,2]","m[3,2]","m[1,3]","m[2,3]","m[3,3]", "tl", "r", "K[1]", "K[2]", "K[3]", "theta"), 
                  trth = c(log(K), m[1,1], m[2,1], m[3,1], m[1,2], m[2,2], m[3,2], m[1,3], m[2,3], m[3,3], tl, r, NA, NA, NA, theta),
                  nincl = rep(0,16), prop = rep(NA,16))
  particular_results_list<-list()
  
  #running simulations (change seed each time)
  for(i in 1:nint){
    print(paste(i, "/", nint))
    seeds[i]<-round(runif(1,1,1000), 0)
    simdat<-sim_recs_nbs(K = K, m_mat = m, tl = tl, cv_dat = cv_dat, theta = theta, catch_dat = catch_mat, release_dat = Rel_mat, nyearPop = Nyear_pop, r=r, s=s, seed = seeds[i])
    
    Kbasin<-simdat$Npop[1,]
    results$trth[13:15]<-Kbasin
    
    N_abunddat<-nrow(simdat$Abund)
    Abund_Est<-simdat$Abund$Npop
    Abund_CV<-rep(cv_dat, N_abunddat)
    
    Years<-seq(min(Catch_dat$Year), max(simdat$Abund$Years), by = 1)
    At_index<-which(Years %in% AtYears)
    In_index<-which(Years %in% IndYears)
    Pac_index<-which(Years %in% PacYears)
    
    NyrA = length(At_index)
    NyrI = length(In_index)
    NyrP = length(Pac_index)
    
    Rec_mat<-simdat$Recs
    
    the_data <- list(ABr=N_abunddat, Nbasin = N_basin, Ibasin = Basin_index,  NyearAbund = Nyear_pop,
                     Nyearcatch = N_yearcatch, Nyeartag = nyeartag, AEstBr = Abund_Est, ACVBr = Abund_CV, NAYearBrA = NyrA,
                     NAYearBrI = NyrI, NAYearBrP = NyrP,
                     AYear = At_index, IYearBr = In_index, PYearBr = Pac_index, 
                     Catch = catch_mat, Rec = Rec_mat, Rel = Rel_mat, ub = 0.34, s = s)
    
    fit <- mod$sample(data = the_data, seed = 123, refresh = 100,
                      iter_warmup = 800, iter_sampling = 800, chains = 2, init=init_fun, adapt_delta = 0.98, max_treedepth = 20)
    
    #after model fit save results
    parsofint<-c("lnK", "m", "tl", "r", "K", "theta")
    #summary statistics 
    sumstats<-fit$summary(parsofint) %>% left_join(fit$summary(parsofint, ~quantile(.x, probs = c(0.025, 0.975)))) %>% rename(q025 = `2.5%`, q975 = `97.5%`)
    #is the true value
    included<-true_val(K = K, m = m, tl = tl, r = r, Kbsn = Kbasin, theta = theta, sumstats = sumstats)
    #increment number included for each parameter
    results$nincl<- results$nincl + included
    #calculate proportion of simulations that include the true value of the parameter
    results$prop<-results$nincl/i
    #save draws and other information about the fit for diagnostics
    particular_results_list$fits[[i]]<-fit$draws(variables = c("MAA","MAB", "MAC", "MBA", "MBB", "MBC", "MCA", "MCB", "MCC", "tl", "r", "lnK", "K", "theta", "lp__"))
    particular_results_list$rhats[[i]]<-rhat(fit) #to check for convergence
    particular_results_list$np[[i]]<-nuts_params(fit) #to check for convergence
    particular_results_list$lps[[i]]<-log_posterior(fit)
  }
  #save seeds
  particular_results_list$seed<-seeds
  return(list(sum_results = results, part_results = particular_results_list))
}

# Simulate and Plot ---------------------------------
#catch data
Catch_dat<-Catch_dat %>% add_row(Year = 1973, Atlantic = 0, Indian = 0, Pacific = 0)
catch_mat<-as.matrix(Catch_dat[,2:4])

#release data
Rel_mat<-Rel_dat %>%
  arrange(Year) %>%
  #select(c("A","B","C")) %>%
  as.matrix()
Rel_mat<-Rel_mat[,2:4]

#Low Movement "true" parameters
K<-200000
m<-matrix(data = NA, nrow = 3, ncol = 3)

m[1,]<-c(0.98, 0.01, 0.01)
m[2,]<-c(0.01, 0.98, 0.01)
m[3,]<-c(0.01, 0.01, 0.98)
Nyear_pop<-96
tl<-0.96
s<-0.96
r<-0.073
theta<-0.8
LowMove<-simfunc(nint = 5, K = K, m = m, tl = tl, r=r, s =s, theta = theta, cv_dat = 0.2, catch_dat = catch_mat, release_dat = Rel_mat)
#save results 
#saveRDS(LowMove, file = "Results/simulationLowMove.RDS")

#High Movement
K<-200000
m<-matrix(data = NA, nrow = 3, ncol = 3)
m[1,]<-c(0.5, 0.25, 0.25)
m[2,]<-c(0.25, 0.5, 0.25)
m[3,]<-c(0.25, 0.25, 0.5)
tl<-0.96
s<-0.96
r<-0.073
HighMove<-simfunc(nint = 5, K = K, m = m, tl = tl, s = s, r = r, cv_dat = 0.2, catch_dat = catch_mat, release_dat = Rel_mat, ntl = 1, nb = F)
#saveRDS(HighMove, file = "Results/simulationHighMove.RDS")


# Simulation Results ------------------------------------------------------
#loading in results (if necessary)
#HighMove<-readRDS("Results/simulationHighMove.RDS")
#LowMove<-readRDS("Results/simulationLowMove.RDS")

#plot movement rates and parameter estimates

#1) High Movement

y<-vector()
y1<-vector()
y2<-vector()
y3<-vector()
y4<-vector()

for(a in 1:9){
  y<-c(y,as.vector(as_draws_matrix(HighMove$part_results$fits[[1]])[,a]))
  y1<-c(y1,as.vector(as_draws_matrix(HighMove$part_results$fits[[2]])[,a]))
  y2<-c(y2,as.vector(as_draws_matrix(HighMove$part_results$fits[[3]])[,a]))
  y3<-c(y3,as.vector(as_draws_matrix(HighMove$part_results$fits[[4]])[,a]))
  y4<-c(y4,as.vector(as_draws_matrix(HighMove$part_results$fits[[5]])[,a]))
}


grp<-rep(c("MAA", "MAB", "MAC", "MBA", "MBB", "MBC", "MCA", "MCB", "MCC"), each = (1600))

tibble_highmove<-tibble(sim1 = y, sim2 = y1, sim3 = y2, sim4 = y3, sim5 = y4, group = grp)
tibble_highmove_long<-tibble_highmove %>% pivot_longer(c(sim1,sim2,sim3,sim4, sim5), names_to = "sim", values_to = "vals")


#check convergence
for(i in 1:5){
  print(mcmc_rhat_hist(ConstHighMove$part_results$rhats[[i]]))
} 

truth_tab<-tibble(grp = as.factor(c("MAA", "MAB", "MAC", "MBA", "MBB", "MBC", "MCA", "MCB", "MCC")), truth = c(0.5, 0.25, 0.25, 0.25, 0.5, 0.25, 0.25, 0.25, 0.5), grpend = c("MAB", "MAC", "MBA", "MBB", "MBC", "MCA", "MCB", "MCC", "MCd"))

plot_sim_high_move<-ggplot(tibble_highmove) + geom_density_ridges(aes(x = sim1, y = group), rel_min_height = 0.02, alpha = 0.4, scale = 1, fill = "cadetblue3", color = "lightsteelblue4") +  
  geom_density_ridges(aes(x = sim2, y = group), rel_min_height = 0.02, alpha = 0.4,scale = 1, fill = "cadetblue3",color = "lightsteelblue4") + 
  geom_density_ridges(aes(x = sim3, y = group), rel_min_height = 0.02, alpha = 0.4, scale = 1, fill = "cadetblue3", color = "lightsteelblue4") +  
  geom_density_ridges(aes(x = sim4, y = group), rel_min_height = 0.02, alpha = 0.4, scale = 1, fill = "cadetblue3", color = "lightsteelblue4") +  
  geom_density_ridges(aes(x = sim5, y = group), rel_min_height = 0.02, alpha = 0.4, scale = 1, fill = "cadetblue3", color = "lightsteelblue4") + 
  scale_y_discrete(breaks = c("MAA", "MAB", "MAC", "MBA", "MBB", "MBC", "MCA", "MCB", "MCC"), labels = c("MAA" = "Atl_Atl", "MAB" = "Atl_Ind", "MAC" = "Atl_Pac", "MBA" = "Ind_Atl", "MBB" = "Ind_Ind", "MBC" = "Ind_Pac", "MCA" = "Pac_Atl", "MCB" = "Pac_Ind", "MCC" = "Pac_Pac"), expand = c(0,0.05)) +
  scale_x_continuous(expand = c(0,0)) +
  labs(x = "Probability", y = "Movement Parameter (from_to)") + 
  geom_segment(data = truth_tab, aes(x = truth, xend = truth, y = grp, yend = grpend)) + 
  theme_classic() + ggtitle("(a) Constant Mark Loss High Movement")
plot_sim_high_move




#2) Low Movement
y<-vector()
y1<-vector()
y2<-vector()
y3<-vector()
y4<-vector()

#extracting results (full draws)
for(a in 1:9){
  y<-c(y,as.vector(as_draws_matrix(LowMove$part_results$fits[[1]])[,a]))
  y1<-c(y1,as.vector(as_draws_matrix(LowMove$part_results$fits[[2]])[,a]))
  y2<-c(y2,as.vector(as_draws_matrix(LowMove$part_results$fits[[3]])[,a]))
  y3<-c(y3,as.vector(as_draws_matrix(LowMove$part_results$fits[[4]])[,a]))
  y4<-c(y4,as.vector(as_draws_matrix(LowMove$part_results$fits[[5]])[,a]))
}


grp<-rep(c("MAA", "MAB", "MAC", "MBA", "MBB", "MBC", "MCA", "MCB", "MCC"), each = (2000)) #2000 is number of iterations for the chains, so change if necessary 

tibble_lowmove<-tibble(sim1 = y, sim2 = y1, sim3 = y2, sim4 = y3, sim5 = y4, group = grp)
tibble_lowmove_long<-tibble_lowmove %>% pivot_longer(c(sim1,sim2,sim3,sim4, sim5), names_to = "sim", values_to = "vals")

#checking convergence
for(i in 1:5){
  print(mcmc_rhat_hist(LowMove$part_results$rhats[[i]]))
} 

#looking at specific simulations
mcmc_rhat_hist(LowMove$part_results$rhats[[1]])
mcmc_trace(LowMove$part_results$fits[[1]], facet_args = list(ncol = 3))


truth_tab_low<-tibble(grp = as.factor(c("MAA", "MAB", "MAC", "MBA", "MBB", "MBC", "MCA", "MCB", "MCC")), truth = c(0.98, 0.01, 0.01, 0.01, 0.98, 0.01, 0.01, 0.01, 0.98), grpend = c("MAB", "MAC", "MBA", "MBB", "MBC", "MCA", "MCB", "MCC", "MCd"))

plot_sim_low_move<-ggplot(tibble_lowmove) + geom_density_ridges(aes(x = sim1, y = group), rel_min_height = 0.02, alpha = 0.4, scale = 1, fill = "cadetblue3", color = "lightsteelblue4") +  
  geom_density_ridges(aes(x = sim2, y = group), rel_min_height = 0.02, alpha = 0.4,scale = 1, fill = "cadetblue3",color = "lightsteelblue4") + 
  geom_density_ridges(aes(x = sim3, y = group), rel_min_height = 0.02, alpha = 0.4, scale = 1, fill = "cadetblue3", color = "lightsteelblue4") +  
  geom_density_ridges(aes(x = sim4, y = group), rel_min_height = 0.02, alpha = 0.4, scale = 1, fill = "cadetblue3", color = "lightsteelblue4") +  
  geom_density_ridges(aes(x = sim5, y = group), rel_min_height = 0.02, alpha = 0.4, scale = 1, fill = "cadetblue3", color = "lightsteelblue4") + 
  scale_y_discrete(breaks = c("MAA", "MAB", "MAC", "MBA", "MBB", "MBC", "MCA", "MCB", "MCC"), labels = c("MAA" = "Atl_Atl", "MAB" = "Atl_Ind", "MAC" = "Atl_Pac", "MBA" = "Ind_Atl", "MBB" = "Ind_Ind", "MBC" = "Ind_Pac", "MCA" = "Pac_Atl", "MCB" = "Pac_Ind", "MCC" = "Pac_Pac"), expand = c(0,0.05)) +
  scale_x_continuous(expand = c(0,0)) +
  labs(x = "Probability", y = "Movement Parameter (from_to)") + 
  geom_segment(data = truth_tab_low, aes(x = truth, xend = truth, y = grp, yend = grpend)) + 
  theme_classic() +  ggtitle("(b) Constant Mark Loss Low Movement")
plot_sim_low_move

#summary results
HighMove$sum_results 
LowMove$sum_results 