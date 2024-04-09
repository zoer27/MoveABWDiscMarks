#Figures for Manuscript (revised)
#Zoe Rand
#Catch data and raw Discovery mark data are available from statistics@iwc.int

# Libraries ---------------------------------------------------------------
library(tidyverse)
library(maptools)
library(sf)
library(PNWColors) #for color palette
library(geomtextpath) #basin labels on maps
#devtools::install_github("LKremer/ggpointdensity") #newer version is faster
library(ggpointdensity)
library(marmap) #for bathymetry
library(patchwork)
library(posterior)
library(grid) #for figure 5
library(ggplotify) #for figure 5
library(scales) #for nicer labels
library(cowplot) #for figure 5
library(bayesplot) #for posterior plots
library(ggridges) #for posterior plots
library(KScorrect) #for log uniform distribution
# Functions ---------------------------------------------------------------

source("Functions/plotcredibleintervalfunc.R") #function from https://medium.com/@jireh/a-clever-use-of-ggplot-internals-bbb168133909

# Data --------------------------------------------------------------------
ANT_BWCatch<-read_csv("Data/ABW_catch_clean.csv") #includes exact locations
Catch<-read_csv("Data/catchcorrected.csv")
ABW_norecaps<-read_csv("Data/ABW_marks_norecaps_clean.csv")
ABW_recaps<-read_csv("Data/ABW_marks_recaps_diffseasons.csv")
ABW_same_season<-read_csv("Data/ABW_marks_recaps_sameseasons.csv")
#Abund_dat<-read_csv("Custom_Model/Data/basin_abundances_withMat.csv")
Abund_dat<-read_csv("Data/basin_abundances_withHamabe.csv")
Rec_dat<-read_csv("Data/mark_recoveries_3groups.csv")

ABW_norecaps<-ABW_norecaps %>%
  mutate(MARK_basin = recode(MARK_basin, A = "Atlantic", B = "Indian", C = "Pacific"))

Years<-seq(min(Catch$Year), max(Abund_dat$Year), by = 1)

Catch_long<-Catch %>% pivot_longer(c("Atlantic", "Indian", "Pacific"), names_to = "basin", values_to = "catch")
tag_groups<-c("AA", "AB", "AC", "BA", "BB","BC", "CA", "CB","CC")
Rec_long<-Rec_dat %>% pivot_longer(all_of(tag_groups), names_to = "group", values_to = "recoveries") 

#Colors
#for basins
#from ColorBrewer
clrs<-c('#1b9e77','#d95f02','#7570b3')
#for other parameters (not basin specific)
#constant tag loss:
pal<-pnw_palette("Starfish", 5)
bckgrd<-"transparent"
#Basins
#Atlantic is -67.26 to 20 Longitudes
#Indian is 20 to 146.9167 Longitude
#Pacific is everything else (146.9 to -67) 
cutpoints<-c(-67.267, 20, 146.9167)
cutpoints_IWC<-c(-170, -120, -60, 0, 70, 130)
# Model Fits --------------------------------------------------------------

fit<-readRDS("FINAL Code/Results/FitNB_Hamabe_ub05_112323.RDS")
draws_array<-fit$draws()
draws_mat<-as_draws_matrix(draws_array)
#thinned
nthin <- 20 #takes every nthin draw
StanThinned<-thin_draws(draws_array, nthin)

#Simulations
HighMove1<-readRDS("Custom_Model/Simulation Results/simulationHighMove_51523.RDS")
HighMove2<-readRDS("Custom_Model/Simulation Results/simulationHighMove_51823.RDS")
HighMove3<-readRDS("Custom_Model/Simulation Results/simulationHighMove_71723.RDS")
LowMove1<-readRDS("Custom_Model/Simulation Results/simulationLowMove_11323.RDS")


# Figure 2 (Map of catches and mark recoveries) ---------------------------

Figure2<-function(){
  ##World Map
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
  
  #catch densities
  #adjust bigger makes it smoother--1 uses default
  #after_stat(ndensity)means that 1 is the max and 0 is the least dense
  g1<-g0 + ggpointdensity::geom_pointdensity(data = ANT_BWCatch, aes (x = Longitude, y = Latitude, group = NULL, color = after_stat(ndensity)), adjust = 1, method = "kde2d", alpha = 0.8) + 
    scale_color_viridis_c(trans = "log", labels = scales::number_format(accuracy = NULL)) + 
    guides(color = guide_colourbar(barwidth = 2, barheight = 10, nbin = 300, draw.llim = TRUE, draw.ulim = TRUE, label.theme = element_text(size = 10), title.theme = element_text(size = 12))) + 
    labs(color = "Relative density")
  g1
  
  #lats and longs for lines
  IWC_dat<-tibble(LatFrom = rep(-90, length(cutpoints_IWC)), LatTo = rep(-70,length(cutpoints_IWC)),Lon = cutpoints_IWC)
  Basin_dat<-tibble(LatFrom = rep(-90, length(cutpoints)), LatTo = rep(-50,length(cutpoints)),Lon = cutpoints)
  
  #also adding basin labels
  
  basin_labs<-tibble(x = c(-60, 60,-120), xend = c(5,120,-170), y = c(-42.5, -42.5, -42.5), yend = c(-43, -43, -43), label = c("Atlantic", "Indian", "Pacific"))
  basin_labs$Basin<-ifelse(basin_labs$label == "Atlantic", "A", ifelse(basin_labs$label == "Indian", "B", "C"))
  
  #adding basin labels and IWC areas to catch density
  g2 <- g1+ 
    geom_segment(data = Basin_dat, aes(x = Lon, xend = Lon, y = LatFrom, yend = LatTo), inherit.aes = F) + 
    geom_textsegment(data = basin_labs, aes(x = x, y = y, xend = xend, yend = yend, label = label, group = NULL), 
                     linecolor = "transparent", linetype = "dotted",vjust = c(0.9, 0.9, 0.1), family = "Arial", size = rel(5)) + 
    geom_segment(data = IWC_dat, aes(x = Lon, xend = Lon, y = LatFrom, yend = LatTo), 
                 color = "gray40", linetype = "dashed", inherit.aes = F) + 
    annotate("text", x = c(0, 130, -120, -60, 70, -170, 173,173), 
             y = c(-50,-48,-48, -48, -48, -48, -50,-60), 
             hjust=c(0.5,0.1,0.9,0.7, 0.1, 0.5, 0.5,0.3),
             vjust=c(-0.4,0.9,0.9,0.7, 0.9,0.9,0.5,0.5),
             label=c("0°","130°E", "120°W","60°W", "70°E", "170°W","50°S", "60°S"), 
             color='gray40', size=rel(2), family = "Arial") + 
    annotate("text", x = c(-35, 45, 100, 170, -145, -95), 
             y = c(-80, -80, -80, -80, -80, -80),
             hjust = c(0.5,0.5, 0.5, 0.5, 0.5, 0.5),
             vjust = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
             label = c("II", "III", "IV", "V", "VI", "I"),
             color = "black", size = rel(3), family = "Arial")+
    labs(tag = "a) Catches") + 
    theme(plot.tag = element_text(size = rel(8)))
  
  g2
  
  
 
  # catch locations--using jitter to represent uncertainty
  # and mark locations
  g3<-g0 + 
    geom_jitter(data=ANT_BWCatch, aes(x=Longitude, y=Latitude, group=NULL), pch = 21, color = "gray70", fill = "gray90",
                       size=0.8, alpha = 0.3, width = 0.5, height = 0.5) + 
    geom_point(data=ABW_norecaps, aes(x=Longitude, y=Latitude, group=NULL, color = MARK_basin, fill = MARK_basin),
                      size=0.8, shape = 24, alpha = 0.6) + 
    scale_color_manual(values = clrs) + scale_fill_manual(values = clrs) +
    geom_segment(data = Basin_dat, aes(x = Lon, xend = Lon, y = LatFrom, yend = LatTo), inherit.aes = F) + 
    annotate("text", x = c(173,173, -67.26, 20, 147), 
             y = c(-50,-60, -48, -48, -48), 
             hjust=c(0.5,0.3, 0.9, 0.5,0.1),
             vjust=c(0.5,0.5, 0.9, 0.7, 0.9),
             label=c("50°S", "60°S", "67.26°W", "20°E", "146.92°E"), 
             color='gray40', size=rel(2), family = "Arial") +
    geom_textsegment(data = basin_labs, aes(x = x, y = y, xend = xend, yend = yend, label = label, color = as.factor(label), group = NULL), 
                     linecolor = "transparent", linetype = "dotted",vjust = c(0.9, 0.9, 0.1), family = "Arial", size = rel(5)) + theme(legend.position = "None",
                                                                                                                                       plot.tag = element_text(size = rel(8))) +
    labs(tag = "b) Mark releases")
    
    
    
  g3
  
  #marking and recovery locations for inter-season recoveries
  BasinIDs<-ABW_recaps %>% dplyr::select(`MarkNo.`, MARK_basin)
  ABW_recaps_plot<-ABW_recaps %>% dplyr::select(c(`MarkNo.`, lat_mark, long_mark, Lat_rec, Long_rec, MARK_basin, recMARK_basin)) %>%
    rename(lat_rec = Lat_rec, long_rec = Long_rec) %>% pivot_longer(starts_with("lat"), names_to = "dir", names_prefix = "lat_", values_to = "Lat") %>% 
      pivot_longer(starts_with("long"), names_to = "dir2", names_prefix = "long_", values_to = "Long") %>% filter(dir == dir2) %>% 
    pivot_longer(c(MARK_basin, recMARK_basin), names_to = "dir3",  values_to = "Basin") %>%
    filter(!(dir == "rec" & dir3 !="recMARK_basin")) %>% filter(!(dir == "mark" & dir3 !="MARK_basin")) %>% arrange(Long, dir) %>%
    group_by(MarkNo.) %>% left_join(BasinIDs) %>%
    ungroup() %>% arrange(dir) %>% mutate(LongNew = Long)
     
  
  MarkIDs<-ABW_recaps_plot %>% dplyr::select(`MarkNo.`) %>% arrange(`MarkNo.`) %>% unlist() %>% unique()
  
  ABW_recaps_sub1<-ABW_recaps_plot %>% filter(`MarkNo.` == MarkIDs[12]) %>% arrange(dir) %>%mutate(LongNew = ifelse(Long < 0, 180+abs(-180-Long), Long))
  ABW_recaps_sub2<-ABW_recaps_plot %>% filter(`MarkNo.` == MarkIDs[15]) %>% arrange(dir) %>%mutate(LongNew = ifelse(Long < 0, 180+abs(-180-Long), Long))
  ABW_recaps_sub3<-ABW_recaps_plot %>% filter(`MarkNo.` == MarkIDs[16]) %>% arrange(dir) %>%mutate(LongNew = ifelse(Long < 0, 180+abs(-180-Long), Long))
  ABW_recaps_sub4<-ABW_recaps_plot %>% filter(`MarkNo.` == MarkIDs[17]) %>% arrange(dir) %>%mutate(LongNew = ifelse(Long < 0, 180+abs(-180-Long), Long))

  ABW_recaps_plot[ABW_recaps_plot$MarkNo. == MarkIDs[12], "LongNew"] <-ABW_recaps_sub1$LongNew
  ABW_recaps_plot[ABW_recaps_plot$MarkNo. == MarkIDs[15], "LongNew"] <-ABW_recaps_sub2$LongNew
  ABW_recaps_plot[ABW_recaps_plot$MarkNo. == MarkIDs[16], "LongNew"] <-ABW_recaps_sub3$LongNew
  ABW_recaps_plot[ABW_recaps_plot$MarkNo. == MarkIDs[17], "LongNew"] <-ABW_recaps_sub4$LongNew


  g4<-g0 + geom_point(data = ABW_recaps_plot, aes(x = LongNew, y = Lat, color = MARK_basin, shape = dir, group = `MarkNo.`)) + 
      geom_path(data = ABW_recaps_plot, aes(x = LongNew, y = Lat, group = `MarkNo.`, color = MARK_basin)) + theme(legend.position = "None") + 
      scale_color_manual(values = clrs) + scale_shape_manual(values = c(17, 19))
  
  
  
  g5<-g4 + 
    annotate("text", x = c(173,173), 
             y = c(-50,-60), 
             hjust=c( 0.5,0.3),
             vjust=c(0.5,0.5),
             label=c("50°S", "60°S"), 
             color='gray40', size=rel(3), family = "Arial") + 
    geom_segment(data = Basin_dat, aes(x = Lon, xend = Lon, y = LatFrom, yend = LatTo), inherit.aes = F) + 
    geom_textsegment(data = basin_labs, aes(x = x, y = y, xend = xend, yend = yend, label = label, color = as.factor(Basin), group = NULL), 
                     linecolor = "transparent", linetype = "dotted",vjust = c(0.9, 0.9, 0.1), family = "Arial", size = rel(5)) + 
    labs(tag = "d) Between seasons", caption = "2/46 < 250 km") +
    theme(legend.position = "None", plot.tag = element_text(size = rel(8)),
          plot.caption = element_text(size = rel(8)))
  
  g5
  
  #marking and recovery locations for same-season recoveries
  BasinIDs<-ABW_same_season %>% dplyr::select(`MarkNo.`, MARK_basin)
  ABW_ss_plot<-ABW_same_season %>% dplyr::select(c(`MarkNo.`, lat_mark, long_mark, Lat_rec, Long_rec, MARK_basin, recMARK_basin)) %>%
    rename(lat_rec = Lat_rec, long_rec = Long_rec) %>% pivot_longer(starts_with("lat"), names_to = "dir", names_prefix = "lat_", values_to = "Lat") %>% 
    pivot_longer(starts_with("long"), names_to = "dir2", names_prefix = "long_", values_to = "Long") %>% filter(dir == dir2) %>% 
    pivot_longer(c(MARK_basin, recMARK_basin), names_to = "dir3",  values_to = "Basin") %>%
    filter(!(dir == "rec" & dir3 !="recMARK_basin")) %>% filter(!(dir == "mark" & dir3 !="MARK_basin")) %>% arrange(Long, dir) %>%
    group_by(MarkNo.) %>% left_join(BasinIDs) %>%
    ungroup() %>% arrange(dir) %>% mutate(LongNew = Long)
  
  
  g6<-g0 + geom_point(data = ABW_ss_plot, aes(x = LongNew, y = Lat, color = MARK_basin, shape = dir, group = `MarkNo.`)) + 
    geom_path(data = ABW_ss_plot, aes(x = LongNew, y = Lat, group = `MarkNo.`, color = MARK_basin)) + theme(legend.position = "None") + 
    scale_color_manual(values = clrs) + scale_shape_manual(values = c(17,19))+ annotate("text", x = c(173,173), 
                                                 y = c(-50,-60), 
                                                 hjust=c( 0.5,0.3),
                                                 vjust=c(0.5,0.5),
                                                 label=c("50°S", "60°S"), 
                                                 color='gray40', size=rel(3), family = "Arial") + 
    geom_segment(data = Basin_dat, aes(x = Lon, xend = Lon, y = LatFrom, yend = LatTo), inherit.aes = F) + 
    geom_textsegment(data = basin_labs, aes(x = x, y = y, xend = xend, yend = yend, label = label, color = as.factor(Basin), group = NULL), 
                     linecolor = "transparent", linetype = "dotted",vjust = c(0.9, 0.9, 0.1), family = "Arial", size = rel(5)) + 
    labs(tag = "c) Within seasons", caption = "11/49 < 250 km") +
    theme(legend.position = "None", 
          plot.tag = element_text(size = rel(8)),
          plot.caption = element_text(size = rel(8)))
  
  g6
  
  lay<- "
  AAABBB
  AAABBB
  CCCDDD
  CCCDDD
  "
  g7<-g2 + g3 + g6 + g5 + plot_layout(design = lay) & 
    theme(plot.tag.position = c(0, 1),
          plot.tag = element_text(size = 12, hjust = 0, vjust = 0))
  g7
  
  print(g7)
 return(g7)     
      
}

Figure2()

fig2<-Figure2()
fig2
ggsave("Figures/Figure2.png", dpi = 900)
# Figure 3 ----------------------------------------------------------------

Figure3<-function(){
  catchplot<-ggplot(data = Catch_long) + 
    geom_bar(aes(x = Year, y = catch, fill = basin),stat = "identity") + 
    facet_grid(rows = vars(basin), space = "free_y", scales = "free_y") +
    scale_x_continuous(breaks = seq(1904, 1973, by = 4), 
                       labels = as.character(seq(1904, 1973, by = 4)), 
                       expand = c(0,0)) + scale_fill_manual(values = clrs)+
    labs(y = "Catches", x = "Season") + 
    theme_minimal() + 
    theme(text = element_text(size = 12),
          legend.position = "none", panel.background = element_blank(), 
          plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(), 
          axis.text = element_text(size = rel(1)),
          #axis.ticks = element_line(colour = "#0f3040"), 
          axis.line = element_line(color = "black"),
          axis.title = element_text(size = rel(1.1)), 
          strip.text = element_text(size = rel(1.05))
          #title = element_text(colour = "#0f3040")
    )
  print(catchplot)
  return(catchplot)
  
}

Figure3()

fig3<-Figure3()
fig3
ggsave("Figures/Figure3.png", dpi = 600)
# Figure 4 ----------------------------------------------------------------
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
  #hist(priorK$K)
  
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
    coord_cartesian(xlim = c(130000, 220000)) +
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
  #hist(priorr$r)
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
  print(combKtlandr)
  return(combKtlandr)
}


Figure4()
fig4<-Figure4()
ggsave("Figures/Figure4.png", dpi = 600)

# Figure 5 ----------------------------------------------------------------

Figure5<-function(){
  #turn CVS into sigmas
  Abund_dat$C<-exp(1.96 * sqrt(log(1+Abund_dat$CV^2)))
  Abund_dat$Lower<-Abund_dat$PopSize/Abund_dat$C
  Abund_dat$Upper<-Abund_dat$PopSize*Abund_dat$C
  
  #thin draws so easier to plot
  ThinnedPPBr<-subset_draws(draws_array, variable = "PostPredBr")
  #3/3/23 note these are actually hamabe but didn't want to change variable names
  ThinnedPPMatInd<-subset_draws(draws_array, variable = "PostPredMatInd")
  ThinnedPPMatPac<-subset_draws(draws_array, variable = "PostPredMatPac")
  
  Quants_PPBr<-summarise_draws(ThinnedPPBr, "median", ~quantile(.x, probs = c(0.025, 0.975)))[c("median", "2.5%", "97.5%")]
  QuantsMIQ<-summarise_draws(ThinnedPPMatInd, "median", ~quantile(.x, probs = c(0.025, 0.975)))[c("median", "2.5%", "97.5%")]
  QuantsMPQ<-summarise_draws(ThinnedPPMatPac, "median", ~quantile(.x, probs = c(0.025, 0.975)))[c("median", "2.5%", "97.5%")]
  
  
  Quants_PPBr_violin<-ThinnedPPBr %>% as_draws_matrix %>% as_tibble %>% pivot_longer(everything(), names_to = "Par", values_to = "PopSize") %>% add_column(Year = rep(Abund_dat$Year[1:9], 4000), Basin = rep(Abund_dat$Basin[1:9], 4000))
  Quants_MIQ_violin<-ThinnedPPMatInd %>% as_draws_matrix %>% as_tibble %>% pivot_longer(everything(), names_to = "Par", values_to = "PopSize") %>% add_column(Year = rep(Abund_dat$Year[10:19], 4000), Basin = rep(Abund_dat$Basin[10:19], 4000))
  Quants_MPQ_violin<-ThinnedPPMatPac %>% as_draws_matrix %>% as_tibble %>% pivot_longer(everything(), names_to = "Par", values_to = "PopSize") %>% add_column(Year = rep(Abund_dat$Year[20:29], 4000), Basin = rep(Abund_dat$Basin[20:29], 4000))
  
  PopQuantsViolin<-bind_rows(Quants_PPBr_violin, Quants_MIQ_violin, Quants_MPQ_violin) %>% add_column(Source = c(rep(Abund_dat$Source[1:9], 4000), rep(Abund_dat$Source[10:19], 8000)))
  
  source.labs<-c("Branch 2007", "Hamabe et al. 2023")
  
  names(source.labs)<-c("Branch", "Hamabe")
  
  bckgrd<-"transparent"
  PPAbund<-ggplot() + geom_violin(data = PopQuantsViolin, aes(x = as.factor(Year), y = PopSize), fill = "gray60", color = "gray60", alpha = 0.5) + 
    geom_point(data = Abund_dat, aes(x = as.factor(Year), y = PopSize, color = Basin, shape = as.factor(Source))) + 
    geom_linerange(data = Abund_dat, aes(x = as.factor(Year), ymin = Lower, ymax = Upper, color = Basin), size = 0.5) +
    facet_grid(Source ~ Basin, scales = "free_x", labeller = labeller(Source = source.labs)) + 
    scale_color_manual(values = clrs, guide = "none") + #scale_fill_manual(values = clrs)+
    scale_shape_manual(values = c(19, 8), guide = "none") + 
    scale_x_discrete(expand = c(0, 1)) + 
    scale_y_continuous(limits = c(0, 2600), breaks = seq(0, 2400, by  = 400)) + 
    labs(y = "Population Size", x = "Year", shape = "Source")  + 
    theme_bw() + theme(
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
  #getting rid of extra panel
  PPAbundgrob<-ggplotGrob(PPAbund)
  
  #gtable_show_layout(PPAbundgrob) 
  #PPAbundgrob$layout$name
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
  #PPAbund_violin<-as.ggplot(~plot(PPAbundgrob))
  #PPAbund_violin

  #population trajectory
  bckgrd<-"white"
  txt<-"black"
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
    labs(y = "Population size")+ 
    scale_y_continuous(labels = comma, breaks = c(0,5000, 10000,50000,75000,100000, 120000), expand = c(0,0))  + 
    scale_x_continuous(expand = c(0,0), breaks = c(seq(1910, 2000, by =10 ))) + expand_limits(y = c(0,120000), x = c(1913, 2001)) +scale_color_manual(values = clrs) + scale_fill_manual(values = clrs) +  
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
  
  #add title
  popplot3<-popplot + ggtitle("Posterior distribution")
  
  #combine
  abundfig2<-ggdraw(popplot3) +
    draw_plot(PPAbundgrob, 0.33, 0.23, 0.65, 0.75) + 
    draw_label(source.labs[1], x = c(0.98), y = c(0.75), angle = c(270), size = 14, fontfamily = "Arial") + 
    draw_label(source.labs[2], x = c(0.98), y = c(0.40), angle = c(270), size = 14, fontfamily = "Arial") +
    draw_label("Survey coverage: 100%", x = c(0.67), y = c(0.86), size = 14) +
    draw_label("Survey coverage: \n 60%", x = c(0.675), y = c(0.58), size = 14) + 
    draw_label("Survey coverage: \n 29%", x = c(0.89), y = c(0.58), size = 14)
  print(abundfig2)
  return(abundfig2)
}

Figure5()
fig5<-Figure5()
ggsave("Figures/Figure5.png", dpi = 900)

# Figure 6 ----------------------------------------------------------------
Figure6<-function(){
  #excluded areas with no recoveries
  Rec_long<-Rec_long %>%
    filter(group %in% tag_groups[c(1:2,4:6,8:9)])
  #adding marking basin
  Rec_long<-Rec_long %>% mutate(mark_basin = if_else("A" == str_sub(group, 1,1), "Atlantic", if_else("B" == str_sub(group, 1,1), "Indian", "Pacific")))
  
  #plot of recovery data
  group_names<-as_labeller(c("AA" = "Atl_Atl", "AB" = "Atl_Ind", "AC" = "Atl_Pac", "BA" = "Ind_Atl", "BB" = "Ind_Ind", "BC" = "Ind_Pac", "CA" = "Pac_Atl","CB" = "Pac_Ind", "CC" = "Pac_Pac"))
  
  
  f1<-ggplot() + geom_bar(data = Rec_long, aes(x = Yr_rec, y = recoveries, fill = mark_basin), stat = "identity" , alpha = 0.8) + 
    scale_x_continuous(expand = c(0,0), limits = c(1926, 1968.05)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
    labs(x = "Season recovered", y = "Number of marks recovered") + scale_fill_manual(values = clrs) + theme_light() + theme(legend.position = "none",                                                                                                                        axis.title = element_text(colour = "#256482"), 
                                                                                                                             panel.border = element_blank(), 
                                                                                                                             axis.title.x = element_text(color = "black"), 
                                                                                                                             axis.title.y  = element_text(color = "black"))
  f1
  
  #adding posterior predictive fit
  nyears<-46
  ThinnedPPRec<-subset_draws(StanThinned, variable = "PostPredRec")
  PPRecAt_At<-summarise_draws(ThinnedPPRec, "median", ~quantile(.x, probs = c(0.025, 0.975)))[1:nyears, c("median", "2.5%", "97.5%")]
  PPRecAt_In<-summarise_draws(ThinnedPPRec, "median", ~quantile(.x, probs = c(0.025, 0.975)))[(nyears+1):(2*nyears), c("median", "2.5%", "97.5%")]
  PPRecAt_Pac<-summarise_draws(ThinnedPPRec, "median", ~quantile(.x, probs = c(0.025, 0.975)))[(2*nyears + 1):(3*nyears), c("median", "2.5%", "97.5%")]
  PPRecIn_At<-summarise_draws(ThinnedPPRec, "median", ~quantile(.x, probs = c(0.025, 0.975)))[(3*nyears + 1):(4*nyears), c("median", "2.5%", "97.5%")]
  PPRecIn_In<-summarise_draws(ThinnedPPRec, "median", ~quantile(.x, probs = c(0.025, 0.975)))[(4*nyears + 1):(5*nyears), c("median", "2.5%", "97.5%")]
  PPRecIn_Pac<-summarise_draws(ThinnedPPRec, "median", ~quantile(.x, probs = c(0.025, 0.975)))[(5*nyears + 1):(6*nyears), c("median", "2.5%", "97.5%")]
  PPRecPac_At<-summarise_draws(ThinnedPPRec, "median", ~quantile(.x, probs = c(0.025, 0.975)))[(6*nyears + 1):(7*nyears), c("median", "2.5%", "97.5%")]
  PPRecPac_In<-summarise_draws(ThinnedPPRec, "median", ~quantile(.x, probs = c(0.025, 0.975)))[(7*nyears + 1):(8*nyears), c("median", "2.5%", "97.5%")]
  PPRecPac_Pac<-summarise_draws(ThinnedPPRec, "median", ~quantile(.x, probs = c(0.025, 0.975)))[(8*nyears + 1):(9*nyears), c("median", "2.5%", "97.5%")]
  
  PPRecAt_At$group<-rep("AA", nrow(PPRecAt_At))
  PPRecAt_In$group<-rep("AB", nrow(PPRecAt_In))
  PPRecAt_Pac$group<-rep("AC", nrow(PPRecAt_Pac))
  PPRecIn_At$group<-rep("BA", nrow(PPRecIn_At))
  PPRecIn_In$group<-rep("BB", nrow(PPRecIn_In))
  PPRecIn_Pac$group<-rep("BC", nrow(PPRecIn_Pac))
  PPRecPac_At$group<-rep("CA", nrow(PPRecPac_At))
  PPRecPac_In$group<-rep("CB", nrow(PPRecPac_In))
  PPRecPac_Pac$group<-rep("CC", nrow(PPRecPac_Pac))
  #combine:
  PPRecStan<-bind_rows(PPRecAt_At, PPRecAt_In, PPRecAt_Pac,PPRecIn_At, PPRecIn_In, PPRecIn_Pac,PPRecPac_At, PPRecPac_In, PPRecPac_Pac)
  PPRecStan$Year<-rep(Years[23:68], 9)
  head(PPRecStan)
  
  PPRecStan<-PPRecStan %>% rename(Lower = "2.5%", Upper = "97.5%")
  #not including 0s and adjusting upper limits so the credible intervals fit in the plot
  NTRStan_filt<-PPRecStan %>% mutate(Upper2 = ifelse(Upper > 6, 6, Upper))
  
  
  
  #facet 
  f2<-f1 + geom_point(data = NTRStan_filt, aes(x = Year, y = median), size = 0.8) + 
    geom_linerange(data = NTRStan_filt, aes(x = Year, ymin = Lower, ymax = Upper2), size = 0.5) +
    facet_grid(rows = vars(group), scales = "free_y", space = "free_y", switch = "y", labeller = group_names) + 
    scale_y_continuous(limits = c(0, 6), expand = expansion(mult = c(0, 0.1))) + 
    theme(text = element_text(size = 18, family = "Arial"),
          strip.background = element_blank(), 
          strip.text = element_text(color = "black", size = rel(0.65)),
          strip.placement = "outside",
          panel.spacing = unit(0.5, "lines"), 
          panel.grid = element_blank())
  print(f2)
  return(f2)
}

Figure6()

fig6<-Figure6()
ggsave("Figures/Figure6.png", dpi = 600)

# Figure 7 ----------------------------------------------------------------
Figure7<- function(){#creating priors
  priors_movemat<-matrix(data = NA, nrow = 1000, ncol = 6)
  for(i in 1:ncol(priors_movemat)){
    priors_movemat[,i]<-runif(1000, min = 0, max = 0.5)
    }
  priors_movemat_y<-matrix(data = NA, nrow = 1000, ncol = 6)


  colnames(priors_movemat)<-c("MAB", "MAC", "MBA", "MBC", "MCA", "MCB")
  
  priors<-as_tibble(priors_movemat) %>% pivot_longer(everything()) %>% add_column(yval = 2)
  
  cols<-c("#65BCB6", "#65BCB6",
          "#65BCB6", "#65BCB6", 
          pal[1], pal[1])
  color_scheme_set(cols)
  m1<-mcmc_areas(
    draws_array, 
    pars = c("MAA","MAB", "MAC", "MBA", "MBB", "MBC", "MCA", "MCB", "MCC"),
    prob = 0.95, 
    point_est = "none",
    area_method = "equal height",
    ) + 
    geom_density_ridges(data = priors, 
                        aes(x = value, y = name, height = yval, group = name), 
                        scale = 0.5,fill = NA, color = "gray40", 
                        linetype = "dashed", rel_min_height = 0.01, stat = "identity") + 
  scale_x_continuous(limits = c(0, 1.02),expand = c(0,0)) + 
  scale_y_discrete(limits = c("MCC", "MCB", "MCA", "MBC", "MBB", "MBA", "MAC","MAB", "MAA"), 
                     labels = c("MAA" = "Atl_Atl", "MAB" = "Atl_Ind", "MAC" = "Atl_Pac", "MBA" = "Ind_Atl", "MBB" = "Ind_Ind", 
                                "MBC" = "Ind_Pac", "MCA" = "Pac_Atl", "MCB" = "Pac_Ind", "MCC" = "Pac_Pac")) + 
  #geom_vline(aes(xintercept = 0.5)) +
    labs(y = "Parameter (from_to)", x = "Movement probability") + 
    theme_classic() + 
    theme(text = element_text(size = 20, family = "Arial"))

  print(m1)
  return(m1)
}

Figure7()

fig7<-Figure7()
ggsave("Figures/Figure7.png", dpi = 900)


# Figure 8 ----------------------------------------------------------------
#Simulation results

Figure8<-function(){
  ntrials<-15 +15 + 30
  
  #movement rates
  get.movement.rates<-function(post_fit){
    y<-as.vector(as_draws_matrix(post_fit)[,1])
    for(a in 2:9){
      y<-c(y,as.vector(as_draws_matrix(post_fit)[,a]))
    }
    return(y)
  }
  
  y1<-lapply(HighMove1$part_results$fits, get.movement.rates)
  y2<-lapply(HighMove2$part_results$fits, get.movement.rates)
  y3<-lapply(HighMove3$part_results$fits, get.movement.rates)
  
  #putting all together 
  y1tab<-bind_cols(y1)
  y2tab<-bind_cols(y2)
  y3tab<-bind_cols(y3)
  clnames<-paste0("sim", 1:ntrials)
  tibble_highmove<-bind_cols(c(y1tab, y2tab, y3tab))
  colnames(tibble_highmove)<-clnames
  grp<-rep(c("MAA", "MAB", "MAC", "MBA", "MBB", "MBC", "MCA", "MCB", "MCC"), each = (1600)) #1600 is number of MCMC draws
  tibble_highmove$group<-grp
  
  
  truth_tab<-tibble(grp = as.factor(c("MCC", "MCB", "MCA", "MBC", "MBB", "MBA", "MAC","MAB", "MAA")), 
                    truth = c(0.5, 0.25, 0.25, 0.25, 0.5, 0.25, 0.25, 0.25, 0.5), 
                    grpend = c("MCB", "MCA", "MBC", "MBB", "MBA", "MAC","MAB", "MAA", "MCd"))
  
  plot_sim_high_move<-ggplot(tibble_highmove) + geom_density_ridges(aes(x = sim1, y = group), rel_min_height = 0.02, alpha = 0.05, scale = 0.8, fill = pal[1], color = pal[1], size = 0.1) +  
    scale_y_discrete(limits = c("MCC", "MCB", "MCA", "MBC", "MBB", "MBA", "MAC","MAB", "MAA", "MCd"), breaks = c("MAA", "MAB", "MAC", "MBA", "MBB", "MBC", "MCA", "MCB", "MCC", "MCd"), labels = c("MAA" = "Atl_Atl", "MAB" = "Atl_Ind", "MAC" = "Atl_Pac", "MBA" = "Ind_Atl", "MBB" = "Ind_Ind", "MBC" = "Ind_Pac", "MCA" = "Pac_Atl", "MCB" = "Pac_Ind", "MCC" = "Pac_Pac", "MCd" = "MCd"), expand = c(0,0.01)) +
    scale_x_continuous(expand = c(0,0), limits = c(0, 1.05)) +
    labs(x = "Movement (m)", y = "density") + 
    theme_classic() +  
    theme(text = element_text(family = "Arial", size = 14), 
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank()) +
    ggtitle("(a)  High movement")
  #plot_sim_high_move
  
  for(i in 2:ntrials){
    plot_sim_high_move<-plot_sim_high_move + geom_density_ridges(data = tibble_highmove, aes_string(x = paste0("sim", i), y = "group"), rel_min_height = 0.02, alpha = 0.05,scale = 0.8, fill = pal[1],color = pal[1], size = 0.1) 
  }
  
  plot_sim_high_move <- plot_sim_high_move +
    annotate("text", x = c(0.83, 0.52, 0.52, 0.52, 0.83, 0.52, 0.52, 0.52, 0.83), 
             y = c("MAA", "MAB", "MAC", "MBA", "MBB", "MBC", "MCA", "MCB", "MCC"), 
             label = c("Atlantic \u2192 Atlantic", "Atlantic \u2192 Indian", "Atlantic \u2192 Pacific", 
                       "Indian \u2192 Atlantic", "Indian \u2192 Indian", "Indian\u2192 Pacific", 
                       "Pacific \u2192 Atlantic", "Pacific \u2192 Indian", "Pacific \u2192 Pacific"), size = 3.5, vjust = c(-2, -1, -1, -1, -3, -1, -1, -1, -3)) +
    geom_segment(data = truth_tab, aes(x = truth, xend = truth, y = grp, yend = grpend))
  
  
  #plot_sim_high_move
  
  #other parameters
  #columns 10-16
  get.other.pars<-function(post_fit){
    y<-as.vector(as_draws_matrix(post_fit)[,10])
    for(a in 11:16){
      y<-c(y,as.vector(as_draws_matrix(post_fit)[,a]))
    }
    return(y)
  }
  
  z1<-lapply(HighMove1$part_results$fits, get.other.pars)
  z2<-lapply(HighMove2$part_results$fits, get.other.pars)
  z3<-lapply(HighMove3$part_results$fits, get.other.pars)
  
  #putting all together 
  z1tab<-bind_cols(z1)
  z2tab<-bind_cols(z2)
  z3tab<-bind_cols(z3)
  
  tibble_otherpars_high<-bind_cols(c(z1tab, z2tab, z3tab))
  colnames(tibble_otherpars_high)<-clnames
  
  grp2<-rep(c("tl", "r", "lnK", "K1", "K2", "K3", "theta"), each = (1600))
  tibble_otherpars_high$group<-grp2
  
  #plot tag loss
  tlplot_high<-tibble_otherpars_high %>% filter(group == "tl") %>% ggplot() + 
    geom_density(aes(x = sim1, y = after_stat(density)),alpha = 0.05, fill = pal[2], color = pal[2], size = 0.1) + 
    theme_classic() + labs(x = "Mark loss (l)", y = "density") + 
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0), limits = c(0, 1.05)) + theme(text = element_text(size = 14, family = "Arial"),
                                                                     axis.text.y = element_blank(), 
                                                                     axis.title.y = element_blank(), 
                                                                     axis.ticks.y = element_blank())
  #tlplot_high
  
  for(i in 2:ntrials){
    tlplot_high<-tlplot_high + geom_density(aes(x = .data[[paste0("sim", i)]], y = after_stat(density)),alpha = 0.05, fill = pal[2], color = pal[2], size = 0.1)
  }
  
  tlplot_high<-tlplot_high + geom_vline(aes(xintercept = 0.96), size = 0.8) 
  
  #tlplot_high
  
  #plot r
  
  rplot_high<-tibble_otherpars_high %>% filter(group == "r") %>% ggplot() + 
    geom_density(aes(x = sim1, y = after_stat(density)),alpha = 0.05, fill = pal[4], color = pal[4], size = 0.1) + 
    theme_classic() + labs(x = "Intrinsic growth (r)", y = "density") + 
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0), limits = c(0, 0.115)) + theme(text = element_text(size = 14, family = "Arial"),
                                                                      axis.ticks.y = element_blank(),
                                                                      axis.title.y = element_blank(), 
                                                                      axis.text.y = element_blank())
  
  for(i in 2:ntrials){
    rplot_high<-rplot_high + geom_density(aes(x = .data[[paste0("sim", i)]], y = after_stat(density)),alpha = 0.05, fill = pal[4], color = pal[4], size = 0.1)
  }
  
  rplot_high<-rplot_high + geom_vline(aes(xintercept = 0.073), size = 0.8)
  
  #rplot_high
  
  
  #plot theta
  thetaplot_high<-tibble_otherpars_high %>% filter(group == "theta") %>% ggplot() + 
    geom_density(aes(x = sim1, y = after_stat(density)),alpha = 0.05, fill = pal[5], color = pal[5], size = 0.1) + 
    theme_classic() + labs(x = "Overdispersion (\u03B8)", y = "density") + 
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0), limits = c(0, 4)) + theme(text = element_text(size = 14, family = "Arial"),
                                                                  axis.text.y = element_blank(),
                                                                  axis.ticks.y = element_blank(),
                                                                  axis.title.y = element_blank())
  
  for(i in 2:ntrials){
    thetaplot_high<-thetaplot_high + geom_density(aes(x = .data[[paste0("sim", i)]], y = after_stat(density)),alpha = 0.05, fill = pal[5], color = pal[5], size = 0.1)
  }
  
  thetaplot_high<-thetaplot_high +  geom_vline(aes(xintercept = 0.8), size = 0.8)
  
  #thetaplot_high
  
  Basinplot_high<-tibble_otherpars_high %>% filter(group %in% c("K1", "K2", "K3")) %>% ggplot() + 
    geom_density_ridges(aes(x = sim1, y = group, fill = group, color = group), rel_min_height = 0.02, alpha = 0.05,scale = 1, size = 0.1) + 
    scale_y_discrete(breaks = c("K1", "K2", "K3"), labels = c("K1" = "Atl", "K2" = "Ind", "K3" = "Pac"), expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0),limits = c(30000, 120000), breaks = c(40000, 60000, 80000,  100000), labels = c("40000","60000", "80000", "100000")) +
    scale_fill_manual(values = clrs) +
    scale_color_manual(values = clrs)+
    labs(x = "Carrying capacity (K)") + 
    theme_classic() +  
    theme(text = element_text(family = "Arial", size = 14),
          legend.position = "none", 
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), 
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = rel(0.8))) 
  
  
  for(i in 2:ntrials){
    Basinplot_high<-Basinplot_high+ geom_density_ridges(aes(x = .data[[paste0("sim", i)]], y = group, color = group, fill = group),alpha = 0.05, rel_min_height = 0.02, scale = 1, size = 0.1)
  }
  
  Basinplot_high <- Basinplot_high + annotate("text", x = c(95000, 95000, 95000), y = c("K3", "K2", "K1"), label = c("Pacific", "Indian", "Atlantic"), vjust = c(-1,-1,-1), family = "Arial", size = 4, color = c(clrs[3], clrs[2], clrs[1]), fontface = "bold") +
    geom_vline(aes(xintercept = 66667), size = 0.8)
  
  #Basinplot_high
  
  
  highplot<-plot_sim_high_move |((rplot_high | tlplot_high)/(thetaplot_high | Basinplot_high)) +
    plot_layout(widths = c(1, 2))
  
  #highplot
  
  #low movement
  
  y1<-lapply(LowMove1$part_results$fits, get.movement.rates)
  
  #remove simulations that did not converge
  y1<-y1[-c(9, 57)]
  #putting all together 
  
  y1tab<-bind_cols(y1)
  clnames<-paste0("sim", 1:ntrials)
  tibble_lowmove<-y1tab
  colnames(tibble_lowmove)<-clnames
  
  grp<-rep(c("MAA", "MAB", "MAC", "MBA", "MBB", "MBC", "MCA", "MCB", "MCC"), each = (2500)) #1600 is number of MCMC draws
  tibble_lowmove$group<-grp
  
  
  truth_tab_low<-tibble(grp = as.factor(c("MCC", "MCB", "MCA", "MBC", "MBB", "MBA", "MAC","MAB", "MAA")), 
                        truth = c(0.98, 0.01, 0.01, 0.01, 0.98, 0.01, 0.01, 0.01, 0.98), 
                        grpend = c("MCB", "MCA", "MBC", "MBB", "MBA", "MAC","MAB", "MAA","MCd"))
  
  plot_sim_low_move<-ggplot(tibble_lowmove) + geom_density_ridges(aes(x = sim1, y = group), rel_min_height = 0.02, alpha = 0.05, scale = 1, fill = pal[1],color = pal[1], size = 0.1) +  
    scale_y_discrete(limits = c("MCC", "MCB", "MCA", "MBC", "MBB", "MBA", "MAC","MAB", "MAA", "MCd"), breaks = c("MAA", "MAB", "MAC", "MBA", "MBB", "MBC", "MCA", "MCB", "MCC"), labels = c("MAA" = "Atl_Atl", "MAB" = "Atl_Ind", "MAC" = "Atl_Pac", "MBA" = "Ind_Atl", "MBB" = "Ind_Ind", "MBC" = "Ind_Pac", "MCA" = "Pac_Atl", "MCB" = "Pac_Ind", "MCC" = "Pac_Pac"), expand = c(0,0.05)) +
    scale_x_continuous(expand = c(0,0), limits = c(0, 1)) +
    labs(x = "Movement (m)", y = "density") + 
    theme_classic() +  
    theme(text = element_text(family = "Arial", size = 14), 
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank()) +
    ggtitle("(b)  Low movement")
  plot_sim_low_move
  
  for(i in 2:(ntrials-2)){
    plot_sim_low_move<-plot_sim_low_move + geom_density_ridges(aes_string(x = paste0("sim", i), y = "group"), rel_min_height = 0.02, alpha = 0.05,scale = 0.8, fill = pal[1],color = pal[1], size = 0.1) 
  }
  
  plot_sim_low_move <- plot_sim_low_move +
    annotate("text", x = c(0.7, 0.25, 0.25, 0.25, 0.7, 0.25, 0.25, 0.25, 0.7), 
             y = c("MAA", "MAB", "MAC", "MBA", "MBB", "MBC", "MCA", "MCB", "MCC"), 
             label = c("Atlantic \u2192 Atlantic", "Atlantic \u2192 Indian", "Atlantic \u2192 Pacific", 
                       "Indian \u2192 Atlantic", "Indian \u2192 Indian", "Indian\u2192 Pacific", 
                       "Pacific \u2192 Atlantic", "Pacific \u2192 Indian", "Pacific \u2192 Pacific"), size = 3.5, vjust = rep(-1.7, 9)) +
    geom_segment(data = truth_tab_low, aes(x = truth, xend = truth, y = grp, yend = grpend)) 
  
  
  #plot_sim_low_move
  
  #other parameters
  #columns 10-16
  
  z1<-lapply(LowMove1$part_results$fits, get.other.pars)
  
  #remove simulations that did not converge
  z1<-z1[-c(9, 57)]
  
  #putting all together 
  z1tab<-bind_cols(z1)
  
  
  clnames<-paste0("sim", 1:(ntrials-2))
  tibble_otherpars_low<-z1tab
  colnames(tibble_otherpars_low)<-clnames
  
  grp2<-rep(c("tl", "r", "lnK", "K1", "K2", "K3", "theta"), each = (2500))
  tibble_otherpars_low$group<-grp2
  
  #plot tag loss
  tlplot_low<-tibble_otherpars_low %>% filter(group == "tl") %>% ggplot() + 
    geom_density(aes(x = sim1, y = after_stat(density)),alpha = 0.05, fill = pal[2], color = pal[2], size = 0.1) + 
    theme_classic() + labs(x = "Mark loss (l)", y = "density") + 
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0), limits = c(0, 1.05)) + theme(text = element_text(size = 14, family = "Arial"),
                                                                     axis.text.y = element_blank(), 
                                                                     axis.title.y = element_blank(), 
                                                                     axis.ticks.y = element_blank())
  #tlplot_low

  
  for(i in 2:(ntrials-2)){
    tlplot_low<-tlplot_low + geom_density(aes(x = .data[[paste0("sim", i)]], y = after_stat(density)),alpha = 0.05, fill = pal[2], color = pal[2], size = 0.1)
  }
  
  tlplot_low<-tlplot_low + geom_vline(aes(xintercept = 0.96), size = 0.8) 
  
  #tlplot_low
  
  #plot r
  
  rplot_low<-tibble_otherpars_low %>% filter(group == "r") %>% ggplot() + 
    geom_density(aes(x = sim1, y = after_stat(density)),alpha = 0.05, fill = pal[4], color = pal[4], size = 0.1) + 
    theme_classic() + labs(x = "Intrinsic growth (r)", y = "density") + 
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0), limits = c(0, 0.118)) + theme(text = element_text(size = 14, family = "Arial"),
                                                                      axis.ticks.y = element_blank(),
                                                                      axis.title.y = element_blank(), 
                                                                      axis.text.y = element_blank())
  
  for(i in 2:(ntrials-2)){
    rplot_low<-rplot_low + geom_density(aes(x = .data[[paste0("sim", i)]], y = after_stat(density)),alpha = 0.05, fill = pal[4], color = pal[4], size = 0.1)
  }
  
  rplot_low<-rplot_low + geom_vline(aes(xintercept = 0.073), size = 0.8)
  
  #rplot_low
  
  
  #plot theta
  thetaplot_low<-tibble_otherpars_low %>% filter(group == "theta") %>% ggplot() + 
    geom_density(aes(x = sim1, y = after_stat(density)),alpha = 0.05, fill = pal[5], color = pal[5], size = 0.1) + 
    theme_classic() + labs(x = "Overdispersion (\u03B8)", y = "density") + 
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0), limits = c(0, 4)) + theme(text = element_text(size = 14, family = "Arial"),
                                                                  axis.text.y = element_blank(),
                                                                  axis.ticks.y = element_blank(),
                                                                  axis.title.y = element_blank())
  
  for(i in 2:(ntrials-2)){
    thetaplot_low<-thetaplot_low + geom_density(aes(x = .data[[paste0("sim", i)]], y = after_stat(density)),alpha = 0.05, fill = pal[5], color = pal[5], size = 0.1)
  }
  
  thetaplot_low<-thetaplot_low +  geom_vline(aes(xintercept = 0.8), size = 0.8)
  
  #thetaplot_low
  
  Basinplot_low<-tibble_otherpars_low %>% filter(group %in% c("K1", "K2", "K3")) %>% ggplot() + 
    geom_density_ridges(aes(x = sim1, y = group, fill = group, color = group), rel_min_height = 0.02, alpha = 0.05,scale = 1, size = 0.1) + 
    scale_y_discrete(breaks = c("K1", "K2", "K3"), labels = c("K1" = "Atl", "K2" = "Ind", "K3" = "Pac"), expand = c(0,0)) +
    scale_x_continuous(limits = c(25000, 115000),expand = c(0,0), breaks = c(40000, 60000, 80000,  100000), labels = c("40000", "60000","80000", "100000")) +
    scale_fill_manual(values = clrs) +
    scale_color_manual(values = clrs)+
    labs(x = "Carrying capacity (K)") + 
    theme_classic() +  
    theme(text = element_text(family = "Arial", size = 14),
          legend.position = "none", 
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), 
          axis.title.y = element_blank(), 
          axis.text.x = element_text(size = rel(0.8))) 
  
  
  
  for(i in 2:(ntrials-2)){
    Basinplot_low<-Basinplot_low+ geom_density_ridges(aes(x = .data[[paste0("sim", i)]], y = group, color = group, fill = group),alpha = 0.05, rel_min_height = 0.02, scale = 1, size = 0.1)
  }
  
  Basinplot_low <- Basinplot_low + annotate("text", x = c(99000, 98000,93000), y = c("K3", "K2", "K1"), label = c("Pacific", "Indian", "Atlantic"), vjust = c(-1,-1,-1), family = "Arial", size = 4, color = c(clrs[3], clrs[2], clrs[1]), fontface = "bold") +
    geom_vline(aes(xintercept = 66667), size = 0.8)
  
  #Basinplot_low
  
  
  lowplot<-plot_sim_low_move |((rplot_low | tlplot_low)/(thetaplot_low | Basinplot_low)) + plot_layout(widths = c(1, 2))
  #lowplot
  #combining into one plot
  comb<-highplot/lowplot 
  print(comb)
  return(comb)
}

Figure8()

fig8<-Figure8()

ggsave("Figures/Figure8.png", fig8, width = 8, height = 9, units = "in", dpi = 600)
