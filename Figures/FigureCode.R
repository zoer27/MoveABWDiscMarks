#### Figures for ABW Movement Manuscript
#3/3/23
#Zoe Rand


# Libraries ---------------------------------------------------------------
library(tidyverse)
library(PNWColors)
library(maptools)
library(maps)
library(mapproj)
library(geomtextpath)
library(patchwork)
# Data --------------------------------------------------------------------
ANT_BWCatch<-read_csv("Data/catches_with_locations.csv") #for map
Catch<-read_csv("Data/catchcorrected.csv")
ABW_norecaps<-read_csv("Data/marks_norecaps_clean.csv") #with locations for map (never recaptured)
ABW_recaps<-read_csv("Data/marks_recaps_diffseasons.csv") #with locations for map (recaptured, excluding same-seaon recoveries)
Abund_dat<-read_csv("Data/basin_abundances_withHamabe.csv")
Rec_dat<-read_csv("Data/mark_recoveries_3groups.csv")

#adding basin names
ABW_norecaps<-ABW_norecaps %>%
  mutate(MARK_basin = recode(MARK_basin, A = "Atlantic", B = "Indian", C = "Pacific"))
#years for plotting
Years<-seq(min(Catch$Year), max(Abund_dat$Year), by = 1)

#converting to long format for plotting
Catch_long<-Catch %>% pivot_longer(c("Atlantic", "Indian", "Pacific"), names_to = "basin", values_to = "catch")
tag_groups<-c("AA", "AB", "AC", "BA", "BB","BC", "CA", "CB","CC")
Rec_long<-Rec_dat %>% pivot_longer(all_of(tag_groups), names_to = "group", values_to = "recoveries") 


# Model Fits ----------------------------------------------------------------
#Model fit
#note: code to run model and simulations must be run before these will be available
fit<-readRDS("Results/FitNB_Hamabe_22323.RDS") #replace with name of saved model fit
draws_array<-fit$draws()
draws_mat<-as_draws_matrix(draws_array)

#thinned--makes it easier to plot
nthin <- 20 #takes every nthin draw
StanThinned<-thin_draws(draws_array, nthin)

#Simulations
LowMove<-readRDS("Results/simulationLowMove.RDS")
ConstHighMove<-readRDS("Results/simulationHighMove.RDS")


# Colors and Functions ----------------------------------------------------------------
#for Basins:
clrs<-c('#1b9e77','#d95f02','#7570b3') #Atlantic, Indian, Pacific

#for other parameters (not basin specific)
pal<-pnw_palette("Starfish", 5)

#to plot credible intervals
source("Functions/plotcredibleintervalfunc.R") #function from https://medium.com/@jireh/a-clever-use-of-ggplot-internals-bbb168133909

# Figure 2 ----------------------------------------------------------------

Figure2<-function(){
  #Atlantic is -67.26 to 20 Longitudes
  #Indian is 20 to 146.9167 Longitude
  #Pacific is everything else (146.9 to -67) 
  cutpoints<-c(-67.267, 20, 146.9167)
  cutpoints_IWC<-c(-170, -120, -60, 0, 70, 130)
  ##World Map
  data(wrld_simpl)
  zlim <- c(0.295,0.705) 
  ## select out just the continent
  worldmap <- wrld_simpl[wrld_simpl$NAME == "Antarctica", ] %>%
    ggplot(aes(x=long, y=lat, group=group)) +
    geom_polygon(fill=gray(0.6)) + 
    scale_y_continuous(limits=c(-90,-42.5))
  g0 <- worldmap + coord_map("ortho", orientation=c(-90, 0, 0)) +
    scale_color_discrete()+
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
          #panel.border=element_rect(color='white'),
          panel.border=element_blank(),
          panel.background = element_rect(fill = "white"),
          panel.grid.major.y = element_line(colour = "grey85", size=0.3),
          ## change background color of facet label
          strip.background =element_rect(fill='white', color='white'),
          ##removing legend                               
          legend.position = "none",
          text = element_text(size = rel(1.5), family = "Arial"))
  g0
  
  #lats and longs for lines
  IWC_dat<-tibble(LatFrom = rep(-90, length(cutpoints_IWC)), LatTo = rep(-70,length(cutpoints_IWC)),Lon = cutpoints_IWC)
  Basin_dat<-tibble(LatFrom = rep(-90, length(cutpoints)), LatTo = rep(-50,length(cutpoints)),Lon = cutpoints)
  
  #adding catch locations--using jitter to represent uncertainty
  #also adding basin labels
  basin_labs<-tibble(x = c(-60, 60,-120), xend = c(5,120,-170), y = c(-42.5, -42.5, -42.5), yend = c(-43, -43, -43), label = c("Atlantic", "Indian", "Pacific"))
  
  g1 <- g0+ geom_jitter(data=ANT_BWCatch, aes(x=Lon, y=Lat, group=NULL), pch = 21, color = "gray70", fill = "gray90",
                        size=0.8, alpha = 0.3, width = 0.5, height = 0.5) + geom_segment(data = IWC_dat, aes(x = Lon, xend = Lon, y = LatFrom, yend = LatTo), inherit.aes = F, linetype = "dashed", color = "gray40") + 
    geom_segment(data = Basin_dat, aes(x = Lon, xend = Lon, y = LatFrom, yend = LatTo), inherit.aes = F) + 
    annotate("text", x = c(0, 130, -120, -60, 70, -170, 173,173), 
             y = c(-50,-48,-48, -48, -48, -48, -50,-60), 
             hjust=c(0.5,0.1,0.9,0.9, 0.1, 0.5, 0.5,0.3),
             vjust=c(-0.4,0.9,0.9,0.9, 0.9,0.9,0.5,0.5),
             label=c("0°","130°E", "120°W","60°W", "70°E", "170°W","50°S", "60°S"), 
             color='gray40', size=rel(5), family = "Arial")
  g1
  
  #adding mark locations
  g2<-g1 + geom_point(data=ABW_norecaps, aes(x=Longitude, y=Latitude, group=NULL, color = MARK_basin, fill = MARK_basin),
                      size=0.8, shape = 24, alpha = 0.6) + scale_color_manual(values = clrs) + scale_fill_manual(values = clrs) +
    geom_textsegment(data = basin_labs, aes(x = x, y = y, xend = xend, yend = yend, label = label, color = as.factor(label), group = NULL), 
                     linecolor = "transparent", linetype = "dotted",vjust = c(0.9, 0.9, 0.9), family = "Arial", size = rel(8)) +
    annotate("text", x = c(-35, 45, 100, 170, -145, -95), 
             y = c(-80, -80, -80, -80, -80, -80),
             hjust = c(0.5,0.5, 0.5, 0.5, 0.5, 0.5),
             vjust = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
             label = c("II", "III", "IV", "V", "VI", "I"),
             color = "gray20", size = rel(6), family = "Arial")
  g2
  
  return(g2)
}

Figure2()


# Figure 3 ----------------------------------------------------------------

Figure3<-function(){
  scales_y <- list(
    Atlantic = scale_y_continuous(limits = c(0,13000), breaks = seq(0,13000, by = 2000)),
    Indian = scale_y_continuous(limits = c(0,16000), breaks = seq(0,16000, by = 2000)),
    Pacific = scale_y_continuous(limits = c(0,5000), breaks = seq(0,5000, by = 1500))
  )
  catchplot<-ggplot(data = Catch_long) + 
    geom_bar(aes(x = Year, y = catch, fill = basin),stat = "identity") + facet_grid_sc(rows = vars(basin), scales = list(y = scales_y), space = "free_y") +
    scale_x_continuous(breaks = seq(1913, 1972, by = 10), labels = as.character(seq(1913, 1972, by = 10)), expand = c(0,0)) + scale_fill_manual(values = clrs)+
    labs(y = "Catches", x = "Season") + theme_minimal() + theme(legend.position = "none",
                                                                panel.grid = element_blank(), 
                                                                text = element_text(size = 18, family = "Arial"),
                                                                strip.text = element_text(size = rel(1.1)),
                                                                panel.background = element_rect(color = bckgrd, fill = bckgrd), 
                                                                plot.background = element_rect(color = bckgrd,fill = bckgrd)
                                                                #axis.line = element_line(color = "gray70")
                                                                )
  
  return(catchplot)
}

Figure3()


# Figure4 -----------------------------------------------------------------

Figure4<-function(){
  #Total Carrying Capacity (K)
  prior_lnK<-matrix(data = runif(1000,9,13), nrow = 1000, ncol = 1)
  colnames(prior_lnK)<-"lnK"
  priorK<-as_tibble(prior_lnK)
  priorK$y<-dunif(priorK$lnK, min = 9, max = 13)
  priorK$K<-exp(priorK$lnK)
  
  lnkpost<-subset_draws(draws_array, "lnK")
  lnKub<-summarise_draws(lnkpost, ~quantile(.x, probs = c(0.025, 0.975)))[,3] %>% rename("q975" = "97.5%")
  lnKlb<-summarise_draws(lnkpost, ~quantile(.x, probs = c(0.025, 0.975)))[,2] %>% rename("q025" = "2.5%")
  lnkmed<-summarise_draws(lnkpost, "median")
  lnkposttab<-lnkpost %>% as_draws_matrix() %>% as_tibble()
  
  kpostdraws<-exp(lnkpost)
  kpostub<-summarise_draws(kpostdraws, ~quantile(.x, probs = c(0.025, 0.975)))[,3] %>% rename("q975" = "97.5%")
  kpostlb<-summarise_draws(kpostdraws, ~quantile(.x, probs = c(0.025, 0.975)))[,2] %>% rename("q025" = "2.5%")
  kposttab<-kpostdraws %>% as_draws_matrix() %>% as_tibble()
  kposttab1<-as.matrix(kposttab)
  kposttab2<-as_tibble(kposttab1)
  
  lnkpostplot<-ggplot() + geom_density(data = lnkposttab, aes(x = lnK, y = after_stat(density)),alpha = 1) + 
    theme_classic() + labs(x = "Carrying capacity (K)", y = "density", tag = "(a)") + 
    geom_density(data = priorK, aes(x = lnK, y = y), stat = "identity", fill = NA, linetype = "dashed") + 
    scale_x_continuous(expand = c(0,0), limits  = c(9,13.1), breaks = c(9, 10, 11, log(193000), 13), 
                       labels = c("8100","22000", "60000", "193000", "442000")) + 
    scale_y_continuous(expand = c(0,0)) +
    theme(axis.text.x = element_text(angle = 90),
          text = element_text(size = 18, family = "Arial"),
          axis.title.x = element_text(vjust=2),
          plot.tag.position = c(0, 1), 
          plot.tag = element_text(vjust = 1.3, hjust = -2.5))
  
  
  lnkpostplot
  
  lnkpostplot1<-lnkpostplot %>% plot_credible_interval(lnKlb$q025, lnKub$q975, fillcol = pal[4])
  lnkpostplot1
  
  #plotting in real space instead of log space
  prior_K<-matrix(data = rlunif(100000,exp(9),exp(13)), nrow = 100000, ncol = 1) #draws from lognormal uniform
  colnames(prior_K)<-"K"
  priorK<-as_tibble(prior_K)
  priorK$y<-dlunif(priorK$K, min = exp(9), max = exp(13))
  hist(priorK$K)
  
  kpostdraws<-exp(lnkpost)
  kpostub<-summarise_draws(kpostdraws, ~quantile(.x, probs = c(0.025, 0.975)))[,3] %>% rename("q975" = "97.5%")
  kpostlb<-summarise_draws(kpostdraws, ~quantile(.x, probs = c(0.025, 0.975)))[,2] %>% rename("q025" = "2.5%")
  kposttab<-kpostdraws %>% as_draws_matrix() %>% as_tibble()
  kposttab1<-as.matrix(kposttab)
  kposttab2<-as_tibble(kposttab1)
  
  kpostplot<-ggplot() + geom_density(data = kposttab2, aes(x = lnK, y = after_stat(density)),alpha = 1) + 
    theme_classic() + labs(x = "Carrying capacity (K)", y = "density", tag = "(a)") + 
    geom_density(data = priorK, aes(x = K, y = after_stat(density)), fill = NA, linetype = "dashed") + 
    #geom_density(data = priorK, aes(x = K, y = y), stat = "identity", fill = NA, linetype = "dashed") + 
    scale_x_continuous(expand = c(0,0)) + 
    coord_cartesian(xlim = c(160000, 220000)) +
    scale_y_continuous(expand = c(0,0)) +
    theme(axis.text.x = element_text(angle = 90),
          text = element_text(size = 18, family = "Arial"),
          axis.title.x = element_text(vjust=2),
          plot.tag.position = c(0, 1), 
          plot.tag = element_text(vjust = 1.3, hjust = -1), 
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank())
  
  kpostplot
  
  kpost1<-kpostplot %>% plot_credible_interval(kpostlb$q025, kpostub$q975, fillcol = pal[4])
  kpost1
  
  prior_r<-matrix(data = runif(1000,0,0.114), nrow = 1000, ncol = 1)
  colnames(prior_r)<-"r"
  priorr<-as_tibble(prior_r)
  hist(priorr$r)
  priorr$y<-dunif(priorr$r, min = 0, max = 0.114)
  
  r_post<-subset_draws(draws_array, "r")
  rub<-summarise_draws(r_post,~quantile(.x, probs = c(0.025, 0.975)))[,3] %>% rename("q975" = "97.5%")
  rlb<-summarise_draws(r_post,~quantile(.x, probs = c(0.025, 0.975)))[,2] %>% rename("q025" = "2.5%")
  rmed<-summarise_draws(r_post, "median")
  rposttab<-r_post %>% as_draws_matrix() %>% as_tibble()
  
  rpost<-ggplot() + geom_density(data = rposttab, aes(x = r, y = after_stat(density)),alpha = 1) + 
    theme_classic() + labs(x = "Intrinsic growth (r)", y = "density", tag = "(b)") + 
    geom_density(data = priorr, aes(x = r, y = y), stat = "identity", fill = NA, linetype = "dashed") + 
    #geom_density(data = priorr, aes(x = r, y = after_stat(density)), fill = NA, linetype = "dashed") + #this is the same as above but the lines are curved strangely
    #geom_vline(aes(xintercept = rmed$median), color = "#3cadbf", size = 2)+
    scale_x_continuous(expand = c(0,0), limits  = c(0,0.115), breaks = c(0.01, 0.03, 0.05, 0.07, 0.09, 0.11)) + 
    scale_y_continuous(expand = c(0,0)) + theme(panel.background = element_rect(fill = bckgrd), 
                                                plot.background = element_rect(fill = bckgrd),
                                                #axis.text = element_text(colour = "#0f3040", size = rel(0.8)), 
                                                #axis.ticks = element_line(colour = "#0f3040"), 
                                                #axis.line = element_line(colour = "#0f3040", size = rel(1)), 
                                                #axis.title = element_text(colour = "#0f3040"), 
                                                panel.border = element_blank(),
                                                text = element_text(size = 18, family = "Arial"),
                                                plot.tag.position = c(0, 1), 
                                                plot.tag = element_text(vjust = 1.3, hjust = -1), 
                                                axis.text.y = element_blank(), 
                                                axis.ticks.y = element_blank())
  
  
  rpost
  
  rpost1<-rpost %>% plot_credible_interval(rlb$q025, rub$q975, fillcol = pal[3])
  rpost1
  
  #tag loss
  prior_tl<-matrix(data = rbeta(1000, 1,1), nrow = 1000, ncol = 1)
  colnames(prior_tl)<-"tl"
  priortl<-as_tibble(prior_tl)
  priortl$yval<-dbeta(priortl$tl, 1, 1)
  hist(priortl$tl)
  #tl
  tl1post<-subset_draws(draws_array, "tl")
  tl1ub<-summarise_draws(tl1post, ~quantile(.x, probs = c(0.025, 0.975)))[,3] %>% rename("q975" = "97.5%")
  tl1lb<-summarise_draws(tl1post, ~quantile(.x, probs = c(0.025, 0.975)))[,2] %>% rename("q025" = "2.5%")
  tl1med<-summarise_draws(tl1post, "median")
  tl1posttab<-tl1post %>% as_draws_matrix() %>% as_tibble()
  
  tl1postplot<-ggplot() + geom_density(data = tl1posttab, aes(x = tl, y = after_stat(density)),alpha = 1) + 
    theme_classic() + labs(x = "Probability of mark loss (l)", y = "density", tag = "(c)")+ 
    geom_density(data = priortl, aes(x = tl, y = yval), stat = "identity", fill = NA, linetype = "dashed") + scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0), limits = c(0.7, 1.01)) + theme(text = element_text(size = 18, family = "Arial"),
                                                                       plot.tag.position = c(0, 1), 
                                                                       plot.tag = element_text(vjust = 1.3, hjust = -1), 
                                                                       axis.text.y = element_blank(), 
                                                                       axis.ticks.y = element_blank())
  
  tl1postplot
  
  tl1postplot1<-tl1postplot %>% plot_credible_interval(tl1lb$q025, tl1ub$q975, fillcol = pal[2])
  tl1postplot1
  
  
  combKtlandr<-kpost1 / (rpost1 + tl1postplot1)
  combKtlandr
}

Figure4()


# Figure 5 ----------------------------------------------------------------
Figure5<-function(){
  #abundance data
  #turn CVS into sigmas
  Abund_dat$C<-exp(1.96 * sqrt(log(1+Abund_dat$CV^2)))
  Abund_dat$Lower<-Abund_dat$PopSize/Abund_dat$C
  Abund_dat$Upper<-Abund_dat$PopSize*Abund_dat$C
  
  bckgrd<-"transparent"
  #library(grid)
  #library(gtable)
  ThinnedQ<-subset_draws(draws_array, variable = "q")
  ThinnedPPBr<-subset_draws(draws_array, variable = "PostPredBr")
  
  #3/3/23 note these are actually hamabe but didn't want to change variable names
  ThinnedPPMatInd<-subset_draws(draws_array, variable = "PostPredMatInd")
  ThinnedPPMatPac<-subset_draws(draws_array, variable = "PostPredMatPac")
  
  
  #PPMatIndQ<-rep(ThinnedQ[,, variable = "q[1]"], each = 10) * ThinnedPPMatInd #don't need this because automatically accounted for in rng
  #PPMatPacQ<-rep(ThinnedQ[,, variable = "q[2]"], each = 10) * ThinnedPPMatPac #don't need this because automatically accounted for in rng
  Quants_PPBr<-summarise_draws(ThinnedPPBr, "median", ~quantile(.x, probs = c(0.025, 0.975)))[c("median", "2.5%", "97.5%")]
  QuantsMIQ<-summarise_draws(ThinnedPPMatInd, "median", ~quantile(.x, probs = c(0.025, 0.975)))[c("median", "2.5%", "97.5%")]
  QuantsMPQ<-summarise_draws(ThinnedPPMatPac, "median", ~quantile(.x, probs = c(0.025, 0.975)))[c("median", "2.5%", "97.5%")]
  
  
  Quants_PPBr_violin<-ThinnedPPBr %>% as_draws_matrix %>% as_tibble %>% pivot_longer(everything(), names_to = "Par", values_to = "PopSize") %>% add_column(Year = rep(Abund_dat$Year[1:9], 4000), Basin = rep(Abund_dat$Basin[1:9], 4000))
  Quants_MIQ_violin<-ThinnedPPMatInd %>% as_draws_matrix %>% as_tibble %>% pivot_longer(everything(), names_to = "Par", values_to = "PopSize") %>% add_column(Year = rep(Abund_dat$Year[10:19], 4000), Basin = rep(Abund_dat$Basin[10:19], 4000))
  Quants_MPQ_violin<-ThinnedPPMatPac %>% as_draws_matrix %>% as_tibble %>% pivot_longer(everything(), names_to = "Par", values_to = "PopSize") %>% add_column(Year = rep(Abund_dat$Year[20:29], 4000), Basin = rep(Abund_dat$Basin[20:29], 4000))
  
  
  
  PopQuantsViolin<-bind_rows(Quants_PPBr_violin, Quants_MIQ_violin, Quants_MPQ_violin) %>% add_column(Source = c(rep(Abund_dat$Source[1:9], 4000), rep(Abund_dat$Source[10:19], 8000)))
  source.labs<-c("Branch 2007", "Hamabe et al. 2023")
  names(source.labs)<-c("Branch", "Hamabe")
  
  PPAbund<-ggplot() + geom_violin(data = PopQuantsViolin, aes(x = as.factor(Year), y = PopSize), fill = "gray60", color = "gray60", alpha = 0.5) + 
    geom_point(data = Abund_dat, aes(x = as.factor(Year), y = PopSize, color = Basin)) + 
    geom_linerange(data = Abund_dat, aes(x = as.factor(Year), ymin = Lower, ymax = Upper, color = Basin), size = 0.5) +
    facet_grid(Source ~ Basin, scales = "free_x", labeller = labeller(Source = source.labs)) + 
    scale_color_manual(values = clrs) + #scale_fill_manual(values = clrs)+
    scale_x_discrete(expand = c(0, 1)) + 
    scale_y_continuous(limits = c(0, 2600), breaks = seq(0, 2400, by  = 400)) + 
    labs(y = "Population Size", x = "Year")  + 
    theme_bw() + theme(legend.position = "none",
                       legend.background = element_rect(fill = bckgrd),
                       panel.background = element_rect(fill = bckgrd), 
                       panel.border = element_blank(),
                       panel.grid = element_blank(),
                       plot.background = element_rect(fill = bckgrd, color = NA),
                       #axis.text = element_text(colour = txt, size = rel(0.8)), 
                       axis.text.x = element_text(angle =90, size = rel(1.2), vjust = 0.5, hjust=1), 
                       axis.text.y = element_text(size = rel(1.2)),
                       axis.ticks = element_line(colour = "gray70"), 
                       axis.line = element_line(colour = "gray70", size = rel(1)), 
                       #axis.title = element_text(size = rel(1.7)),
                       axis.title = element_blank(),
                       strip.background = element_blank(), 
                       strip.text.x = element_text(size = rel(1.5)),
                       strip.text.y = element_blank(),
                       plot.tag.position = "topleft", 
                       plot.title = element_text(size = rel(1.7)), 
                       #text = element_text(family = "Arial")
    ) + 
    ggtitle("Posterior predictive distribution")
  
  
  PPAbund
  
  #getting rid of empty panel
  PPAbundgrob<-ggplotGrob(PPAbund)
  
  #gtable_show_layout(PPAbundgrob) 
  PPAbundgrob$layout$name
  idx <- which(PPAbundgrob$layout$name %in% c("panel-2-1"));
  for (i in idx) PPAbundgrob$grobs[[i]] <- nullGrob();
  #grid.draw(PPAbundgrob)
  
  #move y axis over
  
  idx <- which(PPAbundgrob$layout$name %in% c("axis-l-2"));
  PPAbundgrob$layout[idx, c("l", "r")] <- PPAbundgrob$layout[idx, c("l", "r")] + c(2);
  #grid.draw(PPAbundgrob)
  
  #move x axis up
  idx <- which(PPAbundgrob$layout$name %in% c("axis-b-1"));
  PPAbundgrob$layout[idx, c("t", "b")] <- PPAbundgrob$layout[idx, c("t", "b")] - c(2);
  grid.newpage();
  grid.draw(PPAbundgrob)
  PPAbund_violin<-as.ggplot(~plot(PPAbundgrob))
  PPAbund_violin
  #ggsave("Custom_Model/Figures/PPAbund_bar.png",PPAbundgrob, dpi = 600)
  
  #population trajectory
  bckgrd<-"#FFFFFF"
  txt<-"#000000"
  nyears<-length(Years)
  
  
  ThinnedNPop<-subset_draws(StanThinned, variable = "Npop")
  NpopAt<-summarise_draws(ThinnedNPop, ~quantile(.x, probs = c(0.025, 0.5, 0.975)))[1:nyears, c("50%", "2.5%", "97.5%")] %>% rename(median = "50%", q2.5 = "2.5%", q97.5 = "97.5%")
  NpopIn<-summarise_draws(ThinnedNPop, ~quantile(.x, probs = c(0.025, 0.5, 0.975)))[(nyears+1):(2*nyears), c("50%", "2.5%", "97.5%")] %>% rename(median = "50%", q2.5 = "2.5%", q97.5 = "97.5%")
  NpopPac<-summarise_draws(ThinnedNPop, ~quantile(.x, probs = c(0.025, 0.5, 0.975)))[(2*nyears + 1):(3*nyears), c("50%", "2.5%", "97.5%")] %>% rename(median = "50%", q2.5 = "2.5%", q97.5 = "97.5%")
  
  NpopAt$Basin<-rep("Atlantic", nyears)
  NpopIn$Basin<-rep("Indian", nyears)
  NpopPac$Basin<-rep("Pacific", nyears)
  
  #combine:
  NpopStanAllYr<-bind_rows(NpopAt, NpopIn, NpopPac)
  NpopStanAllYr$Year<-rep(Years,3)
  
  
  
  popplot<-ggplot(data = NpopStanAllYr) + 
    geom_ribbon(aes(x = Year, ymin = q2.5, ymax = q97.5, fill = Basin), alpha = 0.3) + 
    geom_line(aes(x = Year, y = median, color = Basin)) + 
    #geom_point(data = Abund_dat, aes(x = Year, y = PopSize, color = Basin)) + 
    #geom_linerange(data = Abund_dat, aes(x = Year, ymin = Lower, ymax = Upper, color = Basin), linetype = "twodash") +
    labs(y = "Population size")+ 
    scale_y_continuous(labels = comma, breaks = c(0,5000, 10000,50000,75000,100000), expand = c(0,0))  + 
    scale_x_continuous(expand = c(0,0), breaks = c(seq(1910, 2000, by =10 ))) + expand_limits(y = c(0,100000), x = c(1913, 2001)) +scale_color_manual(values = clrs) + scale_fill_manual(values = clrs) +  
    theme_classic() + theme(legend.background = element_rect(fill = bckgrd),
                            panel.background = element_rect(fill = bckgrd), 
                            plot.background = element_rect(fill = bckgrd),
                            axis.text = element_text(size = rel(1.5)), 
                            axis.ticks = element_line(colour = "gray70"), 
                            axis.line = element_line(colour = "gray70", size = rel(1)), 
                            axis.title = element_text(size = rel(1.7)), 
                            panel.border = element_blank(), 
                            legend.position = "bottom", 
                            plot.title = element_text(size = rel(1.7)), 
                            #text = element_text(family = "arial"),
                            legend.text=element_text(size= rel(1.3)), 
                            legend.title = element_text(size = rel(1.3))
    ) 
  popplot
  
  
  popplot3<-popplot + ggtitle("Posterior distribution")
  
  abundfig2<-ggdraw(popplot3 ) +
    draw_plot(PPAbundgrob, 0.33, 0.2, 0.65, 0.75) + 
    draw_label(source.labs[1], x = c(0.98), y = c(0.75), angle = c(270), size = 12, fontfamily = "Arial") + 
    draw_label(source.labs[2], x = c(0.98), y = c(0.40), angle = c(270), size = 12, fontfamily = "Arial")
  return(abundfig2)
}


# Figure 6 ----------------------------------------------------------------



