#Running models with new groups
#8/21/23
#Zoe Rand
#NOTE: Almost all the variables for option 1 are overwritten by stuff for option 2, so need to be careful when running code
#Option 1:
#I/II vs III/IV vs V/VI
opt1<-c("1" = "A", "2" = "A", "3" = "B", "4" = "B", "5" = "C", "6" = "C")
#A is 120W to 0
#B is 0 to 130E
#C is 130E to 20W

#Option 2: 
#II/III, vs IV/V, vs VI/I
opt2<-c("1" = "C", "2" = "A", "3" = "A", "4" = "B", "5" = "B", "6" = "C")
#A is 60W to 70E
#B is 70E to 170W
#C is 170W to 60W

#IWC areas
cutpoints_IWC<-c(-170, -120, -60, 0, 70, 130)

# Libraries ---------------------------------------------------------------
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


# Loading Data ------------------------------------------------------------
#Option 1
Catch_dat<-read_csv("Data/catchcorrected_opt1.csv")
Rec_dat<-read_csv("Data/mark_recoveries_opt1.csv")
Rel_dat<-read_csv("Data/mark_releases_opt1.csv")
Abund_dat<-read_csv("Data/Abund_NewGrps_1.csv")
#using average CV for the missing data point
Abund_dat$CV[is.na(Abund_dat$CV)]<-mean(Abund_dat$CV, na.rm = T)

#Option 2
Catch_dat<-read_csv("Data/catchcorrected_opt2.csv")
Rec_dat<-read_csv("Data/mark_recoveries_opt2.csv")
Rel_dat<-read_csv("Data/mark_releases_opt2.csv")
Abund_dat<-read_csv("Data/Abund_NewGrps_2.csv")

# Functions ---------------------------------------------------------------
source("Functions/Data_Manip_Functions.R") #for stay function

# Data--------------------------------------------

#Basin Information 
N_basin<-3
Basin_index<-c(1,2,3)

#Catch Data
#adding 1973 with 0 catch
Catch_dat<-Catch_dat %>% add_row(Year = 1973, A = 0, B = 0, C = 0)
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


#option 1
#Abundance Estimates
Abund_Est_Br<-Abund_dat$Total[which(Abund_dat$Source == "Branch")]
Abund_Est_Mat_B<-Abund_dat$Total[which(Abund_dat$Source == "Hamabe" & Abund_dat$Area == "B" & Abund_dat$Total != 0)]
Abund_Est_Mat_C<-Abund_dat$Total[which(Abund_dat$Source == "Hamabe" & Abund_dat$Area == "C")]
Abund_CV_Br<-Abund_dat$CV[which(Abund_dat$Source == "Branch")]
Abund_CV_Mat<-Abund_dat$CV[which(Abund_dat$Source == "Hamabe")]
#adding 0 to the end of B
#Abund_Est_Mat_B<-c(Abund_Est_Mat_B, 0)

#Index for year with 0 estimate in Hamabe
Norm_year<-Abund_dat$Year[which(Abund_dat$Total == 0)]


#Year indexes
At_year<-Abund_dat$Year[which(Abund_dat$Area == "A")] #survey years in each area
In_year_Br<-Abund_dat$Year[which(Abund_dat$Area == "B" & Abund_dat$Source == "Branch")] 
Pac_year_Br<-Abund_dat$Year[which(Abund_dat$Area == "C" & Abund_dat$Source == "Branch")] 
In_year_Mat<-Abund_dat$Year[which(Abund_dat$Area == "B" & Abund_dat$Source == "Hamabe" & Abund_dat$Total != 0)]
Pac_year_Mat<-Abund_dat$Year[which(Abund_dat$Area == "C" & Abund_dat$Source == "Hamabe")] 
#Index for year with 0 estimate in Hamabe
Norm_year<-Abund_dat$Year[which(Abund_dat$Total == 0)]
In_year_Mat<-In_year_Mat[In_year_Mat != Norm_year]

Years<-seq(min(Catch_dat$Year), max(Abund_dat$Year), by = 1)
Nyear_pop<-length(Years) 
At_index<-which(Years %in% At_year) #A
In_index_Br<-which(Years %in% In_year_Br) #B
Pac_index_Br<-which(Years %in% Pac_year_Br) #C
In_index_Mat<-which(Years %in% In_year_Mat)#B
Pac_index_Mat<-which(Years %in% Pac_year_Mat) #C
Zero_index_Mat<-which(Years %in% Norm_year) #0



#CV indexes
IndCVIdx<-which(Abund_dat$Source == "Hamabe" & Abund_dat$Area == "B" & Abund_dat$Total !=0 )
PacCVIdx<-which(Abund_dat$Source == "Hamabe" & Abund_dat$Area == "C")
ZeroCVIdx<-which(Abund_dat$Total == 0)
NCVMat<-length(Abund_CV_Mat)

#number of abundance observations 
N_abunddat_Br<-nrow(filter(Abund_dat, Source == "Branch")) 
N_abunddat_Mat_1<-length(Abund_Est_Mat_B) 
N_abunddat_Mat_2<-length(Abund_Est_Mat_C)

#option 2: 
#Abundance Estimates
Abund_Est_Br<-Abund_dat$Total[which(Abund_dat$Source == "Branch")]
Abund_Est_Mat_B<-Abund_dat$Total[which(Abund_dat$Source == "Hamabe" & Abund_dat$Area == "B")]
Abund_CV_Br<-Abund_dat$CV[which(Abund_dat$Source == "Branch")]
Abund_CV_Mat<-Abund_dat$CV[which(Abund_dat$Source == "Hamabe")]
At_year<-Abund_dat$Year[which(Abund_dat$Area == "A")] #survey years in each area
In_year_Br<-Abund_dat$Year[which(Abund_dat$Area == "B" & Abund_dat$Source == "Branch")] 
Pac_year_Br<-Abund_dat$Year[which(Abund_dat$Area == "C" & Abund_dat$Source == "Branch")] 
In_year_Mat<-Abund_dat$Year[which(Abund_dat$Area == "B" & Abund_dat$Source == "Hamabe")]

#number of Abundance observations
N_abunddat_Br<-nrow(filter(Abund_dat, Source == "Branch")) 
N_abunddat_Mat<-nrow(filter(Abund_dat, Source == "Hamabe"))

#CV indexes
IndCVIdx<-which(Abund_dat$Source == "Hamabe" & Abund_dat$Area == "B")
NCVMat<-length(Abund_CV_Mat)
#Year information and indexes
Years<-seq(min(Catch_dat$Year), max(Abund_dat$Year), by = 1)
Nyear_pop<-length(Years) 
At_index<-which(Years %in% At_year) #A
In_index_Br<-which(Years %in% In_year_Br) #B
Pac_index_Br<-which(Years %in% Pac_year_Br) #C
In_index_Mat<-which(Years %in% In_year_Mat)#B
#Pac_index_Mat<-which(Years %in% Pac_year_Mat)


# Running Stan Model ------------------------------------------------------
#option 1


#initial value function 
#MSY fun gives you MSY based on the r and K you draw--helps get reasonable initial values
MSYfun<-function(K, r){
  s = 0.96
  r_star = r/(1-s)
  MSY = ((1-s)*r_star*K*2.39)/((3.39)^((1/2.39) + 1))
  return(MSY)
}
init_fun <- function() {
  rinit = runif(1, 0.06, 0.09)
  out<-list(
    MAB=runif(1, 0.2, 0.33), MAC=runif(1, 0.01, 0.33), MBA=runif(1, 0.2, 0.33),
    MBC=runif(1, 0.2, 0.33), MCA=runif(1, 0.05, 0.33), MCB=runif(1, 0.2, 0.33), 
    tl = runif(1,0.93,1),
    r = rinit,
    #q =runif(2, 0.00001, 1), #adding in q because it's a parameter now
    lnMSY = log(MSYfun(300000, rinit)),
    inv_theta = 1/runif(1, 5, 10))
  return(out)
}
#option 1
#compiling model

file<-"Stan Files/ABWNegBin_newgroups_opt1.stan"
mod<-cmdstan_model(file) 

#estimates in C and B, but there's a 0 in B
the_data <- list(ABr=N_abunddat_Br, AMat = N_abunddat_Mat_1, AMat2 = N_abunddat_Mat_2, AMat3 = 1, 
                 Nbasin = N_basin, Ibasin = Basin_index,  NyearAbund = Nyear_pop,Nyearcatch = N_yearcatch, 
                 Nyeartag = nyeartag, AEstBr = Abund_Est_Br, ACVBr = Abund_CV_Br, NAYearBr = 3, 
                 AEstMatInd = Abund_Est_Mat_B, AEstMatPac = Abund_Est_Mat_C, 
                 NCV = NCVMat, ACVMat = Abund_CV_Mat, IndCVIdx = IndCVIdx, PacCVIdx = PacCVIdx, ZeroCVIdx = ZeroCVIdx, 
                 NAYearMatI = N_abunddat_Mat_1, NAYearMatP = N_abunddat_Mat_2, 
                 AYear = At_index, IYearBr = In_index_Br, PYearBr = Pac_index_Br, IYearMat = In_index_Mat, PYearMat = Pac_index_Mat, 
                 NormYear = Zero_index_Mat, Catch = catch_mat, Rec = Rec_mat, Rel = Rel_mat,ub = 0.499, s = 0.96)

#sampling from model
fit <- mod$sample(data = the_data, seed = 123, refresh = 200,
                  iter_warmup = 2000, iter_sampling = 1000, chains = 4, init=init_fun, adapt_delta = 0.98, max_treedepth = 20)


#fit$save_object(file = "Results/NewGroupsOpt1.RDS")

#option 2
file<-"Stan Files/ABWNegBin_newgroups_opt2.stan"
mod<-cmdstan_model(file) 

#no pacific estimates, no 0 abundances
the_data <- list(ABr=N_abunddat_Br, AMat = N_abunddat_Mat, Nbasin = N_basin, Ibasin = Basin_index,  NyearAbund = Nyear_pop,
                 Nyearcatch = N_yearcatch, Nyeartag = nyeartag, AEstBr = Abund_Est_Br, ACVBr = Abund_CV_Br, NAYearBr = 3, AEstMatInd = Abund_Est_Mat_B, 
                 NCV = NCVMat, ACVMat = Abund_CV_Mat, IndCVIdx = IndCVIdx, 
                 NAYearMatI = 10, AYear = At_index, IYearBr = In_index_Br, PYearBr = Pac_index_Br, IYearMat = In_index_Mat,
                 NormYear = 0, Catch = catch_mat, Rec = Rec_mat, Rel = Rel_mat,ub = 0.5, s = 0.96)

#sampling from model
fit <- mod$sample(data = the_data, seed = 400, refresh = 200,
                  iter_warmup = 500, iter_sampling = 500, chains = 4, init=init_fun, adapt_delta = 0.98, max_treedepth = 20)


#fit$save_object(file = "Results/NewGroupsOpt2.RDS")



# Results -----------------------------------------------------------------
#color palettes:
pal1<-clrs<-c('#1b9e77','#d95f02','#7570b3')
pal2<-c("#ffd700", "#ff1493", "#87cefa")
#map of IWC areas
library(maptools)
data(wrld_simpl)
zlim <- c(0.295,0.705) 
## select out just the continent
worldmap <- wrld_simpl[wrld_simpl$NAME == "Antarctica", ] %>%
  ggplot(aes(x=long, y=lat, group=group)) +
  geom_polygon(fill=gray(0.6)) + 
  scale_y_continuous(limits=c(-90,-42.5))
g0 <- worldmap + coord_map("ortho", orientation=c(-90, 0, 0)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        ## remove space among facets
        panel.spacing.x=unit(0, "lines"),
        panel.spacing.y=unit(0, "lines"),
        ## Change border color
        panel.border=element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid.major.y = element_line(colour = "grey85", size=0.3),
        ## change background color of facet label
        strip.background =element_rect(fill='white', color='white'),
        ##removing legend                               
        #legend.position = "none",
        text = element_text(size = rel(1.5), family = "Arial"))
g0

IWC_dat<-tibble(LatFrom = rep(-90, length(cutpoints_IWC)), LatTo = rep(-50,length(cutpoints_IWC)),Lon = cutpoints_IWC)
#basin model
#Basins: -67.267, 20, 146.9167
g1 <- g0 + annotate("rect", xmin = -67, xmax = 20, ymin = -90, ymax = -50,
                   alpha = .2, fill = pal1[1]) +
  annotate("rect", xmin = 20, xmax = 147, ymin = -90, ymax = -50,
           alpha = .2, fill = pal1[2]) +
  annotate("rect", xmin = 147, xmax = 293, ymin = -90, ymax = -50,
           alpha = .2, fill = pal1[3]) +
  geom_segment(data = IWC_dat, aes(x = Lon, xend = Lon, y = LatFrom, yend = LatTo), 
               color = "gray40", linetype = "dashed", inherit.aes = F) + 
  annotate("text", x = c(0, 130, -120, -60, 70, -170, 173,173), 
           y = c(-50,-48,-48, -48, -48, -48, -50,-60), 
           hjust=c(0.5,0.1,0.9,0.7, 0.1, 0.5, 0.5,0.3),
           vjust=c(-0.4,0.9,0.9,0.7, 0.9,0.9,0.5,0.5),
           label=c("0°","130°E", "120°W","60°W", "70°E", "170°W","50°S", "60°S"), 
           color='gray40', size=rel(5), family = "Arial") + 
  annotate("text", x = c(-35, 45, 100, 170, -145, -95), 
           y = c(-80, -80, -80, -80, -80, -80),
           hjust = c(0.5,0.5, 0.5, 0.5, 0.5, 0.5),
           vjust = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
           label = c("II", "III", "IV", "V", "VI", "I"),
           color = "black", size = rel(7), family = "Arial") +
 labs(tag = "A) Basin model") + theme(plot.tag = element_text(size = rel(10)), 
                                      plot.tag.position = c(0.3, 1))

g1
#Alternate
g2 <- g0 + annotate("rect", xmin = -60, xmax = 70, ymin = -90, ymax = -50,
                    alpha = .2, fill = pal2[1]) +
  annotate("rect", xmin = 70, xmax = 190, ymin = -90, ymax = -50,
           alpha = .2, fill = pal2[2]) +
  annotate("rect", xmin = 190, xmax = 300, ymin = -90, ymax = -50,
           alpha = .2, fill = pal2[3]) +
  geom_segment(data = IWC_dat, aes(x = Lon, xend = Lon, y = LatFrom, yend = LatTo), 
               color = "gray40", linetype = "dashed", inherit.aes = F) + 
  annotate("text", x = c(0, 130, -120, -60, 70, -170, 173,173), 
           y = c(-50,-48,-48, -48, -48, -48, -50,-60), 
           hjust=c(0.5,0.1,0.9,0.7, 0.1, 0.5, 0.5,0.3),
           vjust=c(-0.4,0.9,0.9,0.7, 0.9,0.9,0.5,0.5),
           label=c("0°","130°E", "120°W","60°W", "70°E", "170°W","50°S", "60°S"), 
           color='gray40', size=rel(5), family = "Arial") + 
  annotate("text", x = c(-35, 45, 100, 170, -145, -95), 
           y = c(-80, -80, -80, -80, -80, -80),
           hjust = c(0.5,0.5, 0.5, 0.5, 0.5, 0.5),
           vjust = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
           label = c("II", "III", "IV", "V", "VI", "I"),
           color = "black", size = rel(7), family = "Arial") + 
  labs(tag = "B) Alternative model") + theme(plot.tag = element_text(size = rel(10)), 
                                             plot.tag.position = c(0.3, 1))

g2

library(patchwork)
Figure1<-g1 + g2
Figure1
ggsave("Figures/Newgroups_maps.png", Figure1, dpi = 600)


#option 2 Results
fit_opt2<-readRDS("FINAL Code/Results/NewGroupsOpt2.RDS")
draws_array2<-fit_opt2$draws()
#removing "log_like" variable because I didn't update the code for this so they are NAs
draws_array2<-draws_array2[,,!grepl("log_like", dimnames(draws_array2)$variable)]

#diagnostics
mcmc_trace(fit_opt2$draws(c("lnK", "m", "tl", "r","q", "lnMSY", "theta")), facet_args= list(ncol = 5))
rhats<-rhat(fit_opt2)
mcmc_rhat_hist(rhats) 

#parameter estimates
parsofint<-c("lnK", "K", "m", "tl", "r", "q", "lnMSY", "theta")
sumstats<-fit_opt2$summary(parsofint) %>% left_join(fit_opt2$summary(parsofint, ~quantile(.x, probs = c(0.025, 0.975))))
sumstats
write_csv(sumstats, "FINAL Code/Results/sumstats_newgroups.csv")

#plot movement rates
cols<-c("#65BCB6", "#65BCB6",
        "#65BCB6", "#65BCB6", 
        "#015b58", "#015b58")
color_scheme_set(cols)

m1<-mcmc_areas(
  draws_array2, 
  pars = c("MAA","MAB", "MAC", "MBA", "MBB", "MBC", "MCA", "MCB", "MCC"),
  prob = 0.95, 
  point_est = "none",
  area_method = "equal height",
) + 
  scale_x_continuous(limits = c(0, 1.02),expand = c(0,0)) + 
  scale_y_discrete(limits = c("MCC", "MCB", "MCA", "MBC", "MBB", "MBA", "MAC","MAB", "MAA"), 
                   labels = c("MAA" = "II&III_II&III", "MAB" = "II&III_IV&V", "MAC" = "II&III_VI&I", 
                              "MBA" = "IV&V_II&III", "MBB" = "IV&V_IV&V", 
                              "MBC" = "IV&V_VI&I", "MCA" = "VI&I_II&III", "MCB" = "VI&I_IV&V", "MCC" = "VI&I_VI&I")) + 
  #geom_vline(aes(xintercept = 0.5)) +
  labs(y = "Parameter (from_to)", x = "Movement probability") + 
  theme_classic() + 
  theme(text = element_text(size = 18, family = "Arial"))
m1

#ggsave("Figures/MovementRates_NewGroups.png", m1, dpi = 600)
