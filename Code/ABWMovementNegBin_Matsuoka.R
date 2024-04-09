###ABW Movement Model fitting to Matsuoka and Hakamada 2014 Abundance Data
##Negative binomial likelihood for tags, fixed survival, estimated r
###K in each basin determined by the stationary distribution of m
#Zoe Rand
#2/22/23 
library(tidyverse)
library(cmdstanr)
library(coda)
library(extraDistr)
library(posterior)
library(bayesplot)

#for parallelization
options(mc.cores = parallel::detectCores())

#set_cmdstan_path("/Users/zrand/cmdstan") #to use specific version of stan
cmdstan_version()

# Importing Data ----------------------------------------------------------
Catch_dat<-read_csv("FINAL Code/Data/catchcorrected.csv")
Rec_dat<-read_csv("FINAl Code/Data/mark_recoveries_3groups.csv")
Rel_dat<-read_csv("FINAL Code/Data/mark_releases.csv")
Abund_dat<-read_csv("FINAL Code/Data/basin_abundances_Mat.csv")

# Functions ---------------------------------------------------------------
source("FINAL Code/Functions/Data_Manip_Functions.R") #for stay function
# Data--------------------------------------------

#Basin Information 
N_basin<-3
Basin_index<-c(1,2,3)

#Catch Data
#adding 1973 with 0 catch
Catch_dat<-Catch_dat %>% add_row(Year = 1973, Atlantic = 0, Indian = 0, Pacific = 0)
N_yearcatch<-nrow(Catch_dat)
catch_mat<-as.matrix(Catch_dat[,2:4])


#recovery data as a matrix

Rec_mat<- Rec_dat %>%
  arrange(Yr_rec) %>% #this is actually seasons, 1926 has been removed already
  #select(-1) %>%
  as.matrix()
Rec_mat<-Rec_mat[,-1]

#recovery years
Rec_years<-Rec_dat$Yr_rec
Rec_years<-sort(Rec_years)

#release data as a matrix
Rel_mat<-Rel_dat %>%
  arrange(Year) %>%
  #select(-Year) %>%
  as.matrix()
Rel_mat<-Rel_mat[,-1]
nyeartag<-nrow(Rel_mat)



#Abundance Estimates
Abund_Est_Br<-Abund_dat$PopSize[which(Abund_dat$Source == "Branch")]
Abund_Est_Mat_Ind<-Abund_dat$PopSize[which(Abund_dat$Source == "Matsouka" & Abund_dat$Basin == "Indian")]
Abund_Est_Mat_Pac<-Abund_dat$PopSize[which(Abund_dat$Source == "Matsouka" & Abund_dat$Basin == "Pacific")]
Abund_CV_Br<-Abund_dat$CV[which(Abund_dat$Source == "Branch")]
Abund_CV_Mat<-Abund_dat$CV[which(Abund_dat$Source == "Matsouka")]
At_year<-Abund_dat$Year[which(Abund_dat$Basin == "Atlantic")] #survey years in the Atlantic
In_year_Br<-Abund_dat$Year[which(Abund_dat$Basin == "Indian" & Abund_dat$Source == "Branch")] #survey years in the Indian from Branch
Pac_year_Br<-Abund_dat$Year[which(Abund_dat$Basin == "Pacific" & Abund_dat$Source == "Branch")] #survey yeas in the Pacific 
In_year_Mat<-Abund_dat$Year[which(Abund_dat$Basin == "Indian" & Abund_dat$Source == "Matsouka")]
Pac_year_Mat<-Abund_dat$Year[which(Abund_dat$Basin == "Pacific" & Abund_dat$Source == "Matsouka")] #survey yeas in the Pacific 

N_abunddat_Br<-nrow(filter(Abund_dat, Source == "Branch")) #number of Abundance observations
N_abunddat_Mat<-nrow(filter(Abund_dat, Source == "Matsouka"))

#Year information and indexes
Years<-seq(min(Catch_dat$Year), max(Abund_dat$Year), by = 1)
Nyear_pop<-length(Years) 
At_index<-which(Years %in% At_year)
In_index_Br<-which(Years %in% In_year_Br)
Pac_index_Br<-which(Years %in% Pac_year_Br)
In_index_Mat<-which(Years %in% In_year_Mat)
Pac_index_Mat<-which(Years %in% Pac_year_Mat)



# Running Stan Model ------------------------------------------------------

the_data <- list(ABr=N_abunddat_Br, AMat = N_abunddat_Mat, AMat2 = N_abunddat_Mat/2, Nbasin = N_basin, Ibasin = Basin_index,  NyearAbund = Nyear_pop,
                 Nyearcatch = N_yearcatch, Nyeartag = nyeartag, AEstBr = Abund_Est_Br, ACVBr = Abund_CV_Br, NAYearBr = 3, AEstMatInd = Abund_Est_Mat_Ind, 
                 AEstMatPac = Abund_Est_Mat_Pac, ACVMat = Abund_CV_Mat, NAYearMat = 7, AYear = At_index, IYearBr = In_index_Br, PYearBr = Pac_index_Br, IYearMat = In_index_Mat, PYearMat = Pac_index_Mat, 
                 Catch = catch_mat, Rec = Rec_mat, Rel = Rel_mat,ub = 0.499, s = 0.96)

#initial value function 
#MSY fun gives you MSY based on the r and K you draw--helps get reasonable initial values
MSYfun<-function(K, r){
  s = 0.96
  r_star = r/(1-s)
  MSY = ((1-s)*r_star*K*2.39)/((3.39)^((1/2.39) + 1))
  return(MSY)
}
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

#compiling model

file<-"FINAL Code/ABWMovementNegBinomial.stan"
mod<-cmdstan_model(file) 

#sampling from model
fit <- mod$sample(data = the_data, seed = 400, refresh = 200,
                  iter_warmup = 1000, iter_sampling = 1000, chains = 4, init=init_fun, adapt_delta = 0.98, max_treedepth = 20)





#saving model fit--note with cmdstan need to save using save_object, which saves all of the draws otherwise you might not get all of them

fit$save_object(file = "FINAL Code/Results/FitNB_Matsuoka_ub05_112323.RDS")



parsofint<-c("lnK", "K", "m", "tl", "r", "q", "lnMSY", "theta")
sumstats<-fit$summary(parsofint) %>% left_join(fit$summary(parsofint, ~quantile(.x, probs = c(0.025, 0.975))))
sumstats
#write_csv(sumstats, "Results/sumstats22323.csv")
fit$cmdstan_diagnose()


#converting to stanfit--need to do this for model comparison
library(rstan)
fit2<-rstan::read_stan_csv(fit$output_files())
#saveRDS(fit2, "Results/FitNB22223Stanfit.RDS")




# Diagnostic Plots --------------------------------------------------------
#library(bayesplot)
draws_array<-fit$draws()
#lets you know if any chains didn't make it to the optimal space
apply(draws_array[,,"lp__"], 2, summary) 


lp_fit<-log_posterior(fit)
np_fit<-nuts_params(fit)
rhats<-rhat(fit)

#traceplots
mcmc_trace(fit$draws(c("lnK", "m", "tl", "r","q", "lnMSY", "theta")), facet_args= list(ncol = 5))

#divergences
mcmc_nuts_divergence(np_fit, lp_fit) 
mcmc_nuts_divergence(np_fit, lp_fit, chain = 1) #if there are concerning chains can highlight specific ones

#autocorrelation
mcmc_acf(draws_array, pars = vars(lnK, tl, r), lags = 10)
mcmc_acf(draws_array, pars = vars(MAA, MAB, MAC, MBA, MBB, MBC, MCA, MCB, MCC), lags = 10)

#rhat
mcmc_rhat_hist(rhats) 


#paramter densities
plot1<-mcmc_areas(fit$draws(c("lnK")))
plot1
plot2<-mcmc_areas(fit$draws(c("tl")))
plot2
plot3<-mcmc_areas(fit$draws(c("r")))
plot3
plot4<-mcmc_areas(fit$draws(c("m")))
plot4

#library(patchwork)
plot1 + plot2/plot3 + plot4





