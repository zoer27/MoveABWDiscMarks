#Model Diagnostics
library(tidyverse)
library(cmdstanr)
library(bayesplot)
library(patchwork)

options(mc.cores = parallel::detectCores())

#read in model fit
fit<-readRDS("FINAl Code/Results/FitNB_Matsuoka_ub05_112323.RDS")
#fit<-readRDS("FINAl Code/Results/FitNB_Hamabe_ub05_112323.RDS")

#parsofint<-c("lnK", "K", "m", "tl", "r", "q", "lnMSY") #for constant tag loss
#parsofint<-c("lnK", "K", "m", "tl1930", "tl1944", "tl1957", "r", "q", "lnMSY") #for time varying tag loss
parsofint<-c("lnK", "K", "m", "tl", "r", "q", "lnMSY", "theta") #for negative binomial constant tag loss
#parsofint<-c("lnK", "K", "m", "tl1930", "tl1944", "tl1957", "r", "q", "lnMSY", "theta") #for time varying tag loss
#checking lp values
draws_array<-fit$draws()
apply(draws_array[,,"lp__"], 2, summary) #make sure they all reach similar spaces
#draws_array<-draws_array[,c(1,2,4),] #to subset draws if one chain doesn't converge


#visual diagnostics
lp_fit<-log_posterior(fit)
np_fit<-nuts_params(fit)
rhats<-rhat(fit)

#plots
#par(mfrow=c(4,4),mar=c(3,4,2,1))
#only looking at chains that converged[,c(1,2,4),]
mcmc_trace(fit$draws(c("lnK", "m", "tl", "r","q", "lnMSY", "theta")), facet_args= list(ncol = 5))

#divergences
mcmc_nuts_divergence(np_fit, lp_fit) 
#mcmc_nuts_divergence(np_fit, lp_fit) 
mcmc_parcoord(fit$draws(c("lnK", "m", "tl", "r","q", "lnMSY", "theta")), np = np_fit)

#autocorrelation
mcmc_acf(draws_array, pars = vars(lnK, tl, r, theta), lags = 10)
mcmc_acf(draws_array, pars = vars(MAA, MAB, MAC, MBA, MBB, MBC, MCA, MCB, MCC), lags = 10)

#rhat
mcmc_rhat_hist(rhats) 


#paramter densities
plot1<-mcmc_areas(fit$draws(c("lnK")))
plot1
plot2<-mcmc_areas(fit$draws(c("tl")))
plot2
#plot2<-mcmc_areas(fit$draws(c("tl1930", "tl1944", "tl1957")))
#plot2
plot3<-mcmc_areas(fit$draws(c("r")))
plot3
plot4<-mcmc_areas(fit$draws(c("m")))
plot4

#time varying
(plot1/plot3) | (plot2) | plot4
#for negative binomial
plot5<-mcmc_areas(fit$draws(c("theta")))
plot5
(plot1/plot5) | (plot2/plot3) | plot4

#summary results
sumstats<-fit$summary(parsofint) %>% left_join(fit$summary(parsofint, ~quantile(.x, probs = c(0.025, 0.975))))
sumstats
#write_csv(sumstats, "FINAL Code/Results/sumstats_Hamabe_112323.csv")

#write_csv(sumstats, "FINAL Code/Results/sumstats_Matsuoka_112323.csv")
