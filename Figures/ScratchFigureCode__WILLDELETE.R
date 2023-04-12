
library(maptools)
library(maps)

library(ggplot2)
library(mapproj)
library(geomtextpath)
library(posterior)
library(scales)
library(gridExtra)
library(facetscales)
library(bayesplot)
library(ggridges)
library(gginnards)
library(patchwork)
library(grid)
library(ggplotify)
library(cowplot)
library(PNWColors)
library(KScorrect)
# Data --------------------------------------------------------------------
ANT_BWCatch<-read_csv("Custom_model/Data/catches_with_locations.csv")
Catch<-read_csv("Custom_Model/Data/catchcorrected.csv")
ABW_norecaps<-read_csv("Custom_Model/Data/marks_norecaps_clean.csv")
ABW_recaps<-read_csv("Custom_Model/Data/marks_recaps_diffseasons.csv")
#Abund_dat<-read_csv("Custom_Model/Data/basin_abundances_withMat.csv")
Abund_dat<-read_csv("Custom_Model/Data/basin_abundances_withHamabe.csv")
Rec_dat<-read_csv("Custom_Model/Data/mark_recoveries_3groups.csv")

ABW_norecaps<-ABW_norecaps %>%
  mutate(MARK_basin = recode(MARK_basin, A = "Atlantic", B = "Indian", C = "Pacific"))

Years<-seq(min(Catch$Year), max(Abund_dat$Year), by = 1)

Catch_long<-Catch %>% pivot_longer(c("Atlantic", "Indian", "Pacific"), names_to = "basin", values_to = "catch")
tag_groups<-c("AA", "AB", "AC", "BA", "BB","BC", "CA", "CB","CC")
Rec_long<-Rec_dat %>% pivot_longer(all_of(tag_groups), names_to = "group", values_to = "recoveries") 
# Model Fits --------------------------------------------------------------
#CMD stan object
#fit<-readRDS("Custom_Model/Results/FitNB12522.RDS")
fit<-readRDS("Custom_Model/Results/FitNB_Hamabe_22323.RDS")
draws_array<-fit$draws()
draws_mat<-as_draws_matrix(draws_array)
#thinned
nthin <- 20 #takes every nthin draw
StanThinned<-thin_draws(draws_array, nthin)

#Simulations
ConstLowMove<-readRDS("Custom_Model/Results/simulationConstLowMove121322.RDS")
ConstHighMove<-readRDS("Custom_Model/Results/simulationConstHighMove121322.RDS")

#Time Varying
tv_fit<-readRDS("Custom_Model/Results/timevaryingNB12722.RDS")
tv_draws_array<-tv_fit$draws()
apply(tv_draws_array[,,"lp__"], 2, summary) #chain one just got stuck somewhere so removing it for now
#tv_draws_array<-tv_draws_array[,c(2,3,4),]
tv_draws_mat<-as_draws_matrix(tv_draws_array)
#thinned
nthin<-20
tv_StanThinned<-thin_draws(tv_draws_array, nthin)

# Colors ------------------------------------------------------------------
#for basins
#from ColorBrewer
clrs<-c('#1b9e77','#d95f02','#7570b3')
#for other parameters (not basin specific)
#constant tag loss:
pal<-pnw_palette("Starfish", 5)
pal
#non-constant tag loss
#pal2<-pnw_palette("Bay", 5)
# Functions ---------------------------------------------------------------
source("Figures/plotcredibleintervalfunc.R") #function from https://medium.com/@jireh/a-clever-use-of-ggplot-internals-bbb168133909


# Map of Catches+Marks+IWC+Basins -----------------------------------------

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

#adding basin labels
#basin_labs<-tibble(x = c(-60, 60,-120), xend = c(5,120,-170), y = c(-42.5, -42.5, -42.5), yend = c(-43, -43, -43), label = c("Atlantic", "Indian", "Pacific"))
#g3<-g2 + geom_textsegment(data = basin_labs, aes(x = x, y = y, xend = xend, yend = yend, label = label, color = as.factor(label), group = NULL), 
#linecolor = "transparent", vjust = c(0.9, 0.9, 0.9), family = "Times", size = 6)

#g3
g2
ggsave(g2, filename = "Figures/CatchesAndMarks.png", dpi = 600)

# Fitted Population Trajectories + Catches by basin ------------------------------
#abundance data
#turn CVS into sigmas
Abund_dat$C<-exp(1.96 * sqrt(log(1+Abund_dat$CV^2)))
Abund_dat$Lower<-Abund_dat$PopSize/Abund_dat$C
Abund_dat$Upper<-Abund_dat$PopSize*Abund_dat$C

bckgrd<-"transparent"
library(grid)
library(gtable)
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


#PopQuants<-bind_rows(Quants_PPBr, QuantsMIQ, QuantsMPQ)
#PopQuants<-PopQuants %>% add_column(Year = Abund_dat$Year, Basin = Abund_dat$Basin, Data = Abund_dat$PopSize, 
#Source = Abund_dat$Source) %>%
#rename(Lower = "2.5%", Upper = "97.5%")

PopQuantsViolin<-bind_rows(Quants_PPBr_violin, Quants_MIQ_violin, Quants_MPQ_violin) %>% add_column(Source = c(rep(Abund_dat$Source[1:9], 4000), rep(Abund_dat$Source[10:19], 8000)))
source.labs<-c("Branch 2007", "Hamabe et al. 2023")
names(source.labs)<-c("Branch", "Hamabe")

#PPAbund<-ggplot(PopQuants) + geom_line(aes(x = Year, y = median, color = Basin)) + 
#geom_ribbon(aes(x = Year, ymin = Lower, ymax = Upper, color = Basin), alpha = 0, linetype = "dashed") +
#geom_point(aes(x = Year, y = Data, color = Basin)) + facet_grid(Source ~ Basin, scales = "free_x") + 
#scale_color_manual(values = clrs) + scale_x_continuous(breaks = seq(1981, 2008, by = 5)) + 
#labs(y = "Population Size", x = "Year")  + 
#theme_bw() + theme(legend.position = "none",
#legend.background = element_rect(fill = bckgrd),
#panel.background = element_rect(fill = bckgrd), 
#panel.border = element_blank(),
#panel.grid = element_blank(),
#plot.background = element_rect(fill = bckgrd),
#axis.text = element_text(colour = "#256482", size = rel(0.8)), 
#axis.text.x = element_text(angle =90), 
#axis.ticks = element_line(colour = "#256482"), 
#axis.line = element_line(colour = "#256482", size = rel(1)), 
#axis.title = element_text(colour = "#256482"),
#strip.background = element_blank())

#PPAbund


#As Violin Plot: 
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



#bars with errorpoints
#p1<-ggplot() + geom_bar(data = Abund_dat, aes(x = Year, y = PopSize, fill = Basin), stat = "identity" , alpha = 0.8) + 
#scale_x_continuous(expand = c(0,0), breaks = seq(1979, 2008, by = 5)) +
#labs(x = "Year", y = "Population Size") + scale_fill_manual(values = clrs) + theme_light() 


#p1
#p2<- p1 + geom_point(data = PopQuants, aes(x = Year, y = median), size = 0.8) + 
#geom_linerange(data = PopQuants, aes(x = Year, ymin = Lower, ymax = Upper), size = 0.5) +
#facet_grid(Source ~ Basin, scales = "free_x") + 
#theme(strip.background = element_blank(), panel.spacing = unit(0.5, "lines"))

#p2
#PPAbund<-p2
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
abundfig2

ggsave(abundfig2, file = "Figures/populationmodel.png", dpi = 600, width = 9.5, height = 7, units = "in")

#catches by basin
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
catchplot

ggsave(catchplot, filename = "Figures/CatchBasin.png", dpi = 600)



# Mark Recovery Fits ------------------------------------------------------
#posterior predictive
#posterior predictive for movement

bckgrd<-"#FFFFFF"
#excluded areas with no recoveries
Rec_long<-Rec_long %>%
  filter(group %in% tag_groups[c(1:2,4:6,8:9)])
#adding marking basin
Rec_long<-Rec_long %>% mutate(mark_basin = if_else("A" == str_sub(group, 1,1), "Atlantic", if_else("B" == str_sub(group, 1,1), "Indian", "Pacific")))

#plot of recovery data
group_names<-as_labeller(c("AA" = "Atl_Atl", "AB" = "Atl_Ind", "AC" = "Atl_Pac", "BA" = "Ind_Atl", "BB" = "Ind_Ind", "BC" = "Ind_Pac", "CA" = "Pac_Atl","CB" = "Pac_Ind", "CC" = "Pac_Pac"))


f1<-ggplot() + geom_bar(data = Rec_long, aes(x = Yr_rec, y = recoveries, fill = mark_basin), stat = "identity" , alpha = 0.8) + 
  scale_x_continuous(expand = c(0,0)) +
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
PPRecStan$Year<-rep(Years[15:60], 9)
head(PPRecStan)

PPRecStan<-PPRecStan %>% rename(Lower = "2.5%", Upper = "97.5%")
#not including 0s
NTRStan_filt<-PPRecStan #%>%
#filter(q95 >=0.10)

scales_y2<-list(
  AA = scale_y_continuous(limits = c(0,6), breaks = seq(0,6, by = 1), expand = c(0,0), minor_breaks = 0),
  AB = scale_y_continuous(limits = c(0,3.7), breaks = seq(0,3, by = 1), expand = c(0,0),minor_breaks = 0),
  AC = scale_y_continuous(limits = c(0,3.7),breaks = seq(0,3, by = 1), expand = c(0,0),minor_breaks = 0),
  BA = scale_y_continuous(limits = c(0,3.7),breaks = seq(0,3, by = 1), expand = c(0,0),minor_breaks = 0),
  BB = scale_y_continuous(limits = c(0,5.2),breaks = seq(0,5, by = 1), expand = c(0,0),minor_breaks = 0),
  BC = scale_y_continuous(limits = c(0,3.7),breaks = seq(0,3, by = 1), expand = c(0,0),minor_breaks = 0),
  CA = scale_y_continuous(limits = c(0,3.7), breaks = seq(0,3, by = 1), expand = c(0,0),minor_breaks = 0),
  CB = scale_y_continuous(limits = c(0,3.7),breaks = seq(0,3, by = 1), expand = c(0,0),minor_breaks = 0),
  CC = scale_y_continuous(limits = c(0,3.7),breaks = seq(0,3, by = 1), expand = c(0,0),minor_breaks = 0)
)

f2<-f1 + geom_point(data = NTRStan_filt, aes(x = Year, y = median), size = 0.8) + 
  geom_linerange(data = NTRStan_filt, aes(x = Year, ymin = Lower, ymax = Upper), size = 0.5) +
  facet_grid_sc(rows = vars(group), scales = list(y = scales_y2), space = "free_y", switch = "y", labeller = group_names) + 
  theme(text = element_text(size = 18, family = "Arial"),
        strip.background = element_blank(), 
        strip.text = element_text(color = "black", size = rel(0.65)),
        strip.placement = "outside",
        panel.spacing = unit(0.5, "lines"), 
        panel.grid = element_blank())
f2

ggsave(f2, filename = "Figures/TagModelPostPred.png", dpi = 600)
# Movement posteriors and priors ------------------------------------------
#creating priors
priors_movemat<-matrix(data = NA, nrow = 1000, ncol = 6)
for(i in 1:ncol(priors_movemat)){
  priors_movemat[,i]<-runif(1000, min = 0, max = 0.34)
}
priors_movemat_y<-matrix(data = NA, nrow = 1000, ncol = 6)
priors_movemat_y<-dunif(priors_movemat, min = 0, max = 0.34)

colnames(priors_movemat)<-c("MAB", "MAC", "MBA", "MBC", "MCA", "MCB")
colnames(priors_movemat_y)<-c("MAB", "MAC", "MBA", "MBC", "MCA", "MCB")

priors<-as_tibble(priors_movemat) %>% pivot_longer(everything())
priorsy<-as_tibble(priors_movemat_y) %>% pivot_longer(everything()) %>% rename( yval = value)

priors<-priors %>% left_join(priorsy, by = "name")

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
) + geom_density_ridges(data = priors, aes(x = value, y = name, height = yval, group = name), scale = 0.5, fill = NA, color = "gray40", linetype = "dashed", rel_min_height = 0.01) + 
  scale_x_continuous(limits = c(0, 1.02),expand = c(0,0)) + scale_y_discrete(limits = c("MCC", "MCB", "MCA", "MBC", "MBB", "MBA", "MAC","MAB", "MAA"), labels = c("MAA" = "Atl_Atl", "MAB" = "Atl_Ind", "MAC" = "Atl_Pac", "MBA" = "Ind_Atl", "MBB" = "Ind_Ind", "MBC" = "Ind_Pac", "MCA" = "Pac_Atl", "MCB" = "Pac_Ind", "MCC" = "Pac_Pac")) + 
  labs(y = "Parameter (from_to)", x = "Movement probability") + theme_classic() + theme(text = element_text(size = 18, family = "Arial"))

m1

ggsave(m1, filename = "Figures/movementposts.png", dpi = 600)



# Posteriors and priors for other key pars --------------------------------

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

ggsave(combKtlandr, filename = "Figures/Ktlandrpost.png", dpi = 600)


# Simualation Results ------------------------------------------------------

#1) Constant High Movement

y<-vector()
y1<-vector()
y2<-vector()
y3<-vector()
y4<-vector()

for(a in 1:9){
  y<-c(y,as.vector(as_draws_matrix(ConstHighMove$part_results$fits[[1]])[,a]))
  y1<-c(y1,as.vector(as_draws_matrix(ConstHighMove$part_results$fits[[2]])[,a]))
  y2<-c(y2,as.vector(as_draws_matrix(ConstHighMove$part_results$fits[[3]])[,a]))
  y3<-c(y3,as.vector(as_draws_matrix(ConstHighMove$part_results$fits[[4]])[,a]))
  y4<-c(y4,as.vector(as_draws_matrix(ConstHighMove$part_results$fits[[5]])[,a]))
}

#yrep<-matrix(data = c(y1,y2,y3,y4), nrow = 4, ncol = (600*9), byrow = TRUE)

grp<-rep(c("MAA", "MAB", "MAC", "MBA", "MBB", "MBC", "MCA", "MCB", "MCC"), each = (1600))

tibble_highmove<-tibble(sim1 = y, sim2 = y1, sim3 = y2, sim4 = y3, sim5 = y4, group = grp)
tibble_highmove_long<-tibble_highmove %>% pivot_longer(c(sim1,sim2,sim3,sim4, sim5), names_to = "sim", values_to = "vals")



for(i in 1:5){
  print(mcmc_rhat_hist(ConstHighMove$part_results$rhats[[i]]))
} #all converged

truth_tab<-tibble(grp = as.factor(c("MCC", "MCB", "MCA", "MBC", "MBB", "MBA", "MAC","MAB", "MAA")), 
                  truth = c(0.5, 0.25, 0.25, 0.25, 0.5, 0.25, 0.25, 0.25, 0.5), 
                  grpend = c("MCB", "MCA", "MBC", "MBB", "MBA", "MAC","MAB", "MAA", "MCd"))

plot_sim_high_move<-ggplot(tibble_highmove) + geom_density_ridges(aes(x = sim1, y = group), rel_min_height = 0.02, alpha = 0.3, scale = 0.8, fill = pal[1], color = pal[1]) +  
  geom_density_ridges(aes(x = sim2, y = group), rel_min_height = 0.02, alpha = 0.3,scale = 0.8, fill = pal[1],color = pal[1]) + 
  geom_density_ridges(aes(x = sim3, y = group), rel_min_height = 0.02, alpha = 0.3, scale = 0.8, fill = pal[1], color = pal[1]) +  
  geom_density_ridges(aes(x = sim4, y = group), rel_min_height = 0.02, alpha = 0.3, scale = 0.8, fill = pal[1], color = pal[1]) +  
  geom_density_ridges(aes(x = sim5, y = group), rel_min_height = 0.02, alpha = 0.3, scale = 0.8, fill = pal[1], color = pal[1]) + 
  scale_y_discrete(limits = c("MCC", "MCB", "MCA", "MBC", "MBB", "MBA", "MAC","MAB", "MAA", "MCd"), breaks = c("MAA", "MAB", "MAC", "MBA", "MBB", "MBC", "MCA", "MCB", "MCC", "MCd"), labels = c("MAA" = "Atl_Atl", "MAB" = "Atl_Ind", "MAC" = "Atl_Pac", "MBA" = "Ind_Atl", "MBB" = "Ind_Ind", "MBC" = "Ind_Pac", "MCA" = "Pac_Atl", "MCB" = "Pac_Ind", "MCC" = "Pac_Pac", "MCd" = "MCd"), expand = c(0,0.01)) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 1.05)) +
  labs(x = "Movement (m)", y = "density") + 
  annotate("text", x = c(0.83, 0.52, 0.52, 0.52, 0.83, 0.52, 0.52, 0.52, 0.83), 
           y = c("MAA", "MAB", "MAC", "MBA", "MBB", "MBC", "MCA", "MCB", "MCC"), 
           label = c("Atlantic \u2192 Atlantic", "Atlantic \u2192 Indian", "Atlantic \u2192 Pacific", 
                     "Indian \u2192 Atlantic", "Indian \u2192 Indian", "Indian\u2192 Pacific", 
                     "Pacific \u2192 Atlantic", "Pacific \u2192 Indian", "Pacific \u2192 Pacific"), size = 3.5, vjust = c(-3, -1, -1, -1, -3, -1, -1, -1, -3)) +
  geom_segment(data = truth_tab, aes(x = truth, xend = truth, y = grp, yend = grpend)) + 
  theme_classic() +  
  theme(text = element_text(family = "Arial", size = 14), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank()) +
  ggtitle("(a)  High movement")
plot_sim_high_move

#other parameters
z<-vector()
z1<-vector()
z2<-vector()
z3<-vector()
z4<-vector()

for(a in 10:16){
  z<-c(z,as.vector(as_draws_matrix(ConstHighMove$part_results$fits[[1]])[,a]))
  z1<-c(z1,as.vector(as_draws_matrix(ConstHighMove$part_results$fits[[2]])[,a]))
  z2<-c(z2,as.vector(as_draws_matrix(ConstHighMove$part_results$fits[[3]])[,a]))
  z3<-c(z3,as.vector(as_draws_matrix(ConstHighMove$part_results$fits[[4]])[,a]))
  z4<-c(z4,as.vector(as_draws_matrix(ConstHighMove$part_results$fits[[5]])[,a]))
}

grp<-rep(c("tl", "r", "lnK", "K1", "K2", "K3", "theta"), each = (1600))

tibble_otherpars_high<-tibble(sim1 = z, sim2 = z1, sim3 = z2, sim4 = z3, sim5 = z4, group = grp)
tibble_otherpars_high_long<-tibble_otherpars_high %>% pivot_longer(c(sim1,sim2,sim3,sim4, sim5), names_to = "sim", values_to = "vals")

tlplot_high<-tibble_otherpars_high %>% filter(group == "tl") %>% ggplot() + 
  geom_density(aes(x = sim1, y = after_stat(density)),alpha = 0.3, fill = pal[2], color = pal[2]) + 
  geom_density(aes(x = sim2, y = after_stat(density)),alpha = 0.3, fill = pal[2], color = pal[2]) + 
  geom_density(aes(x = sim3, y = after_stat(density)),alpha = 0.3, fill = pal[2], color = pal[2]) + 
  geom_density(aes(x = sim4, y = after_stat(density)),alpha = 0.3, fill = pal[2], color = pal[2]) + 
  geom_density(aes(x = sim5, y = after_stat(density)),alpha = 0.3, fill = pal[2], color = pal[2]) +
  geom_vline(aes(xintercept = 0.96), size = 0.8) +
  theme_classic() + labs(x = "Mark loss (l)", y = "density") + 
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 1.05)) + theme(text = element_text(size = 14, family = "Arial"),
                                                                   axis.text.y = element_blank(), 
                                                                   axis.title.y = element_blank(), 
                                                                   axis.ticks.y = element_blank())
tlplot_high

rplot_high<-tibble_otherpars_high %>% filter(group == "r") %>% ggplot() + 
  geom_density(aes(x = sim1, y = after_stat(density)),alpha = 0.3, fill = pal[4], color = pal[4]) + 
  geom_density(aes(x = sim2, y = after_stat(density)),alpha = 0.3, fill = pal[4], color = pal[4]) + 
  geom_density(aes(x = sim3, y = after_stat(density)),alpha = 0.3, fill = pal[4], color = pal[4]) + 
  geom_density(aes(x = sim4, y = after_stat(density)),alpha = 0.3, fill = pal[4], color = pal[4]) + 
  geom_density(aes(x = sim5, y = after_stat(density)),alpha = 0.3, fill = pal[4], color = pal[4]) +
  geom_vline(aes(xintercept = 0.073), size = 0.8) +
  theme_classic() + labs(x = "Intrinsic growth (r)", y = "density") + 
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 0.115)) + theme(text = element_text(size = 14, family = "Arial"),
                                                                    axis.ticks.y = element_blank(),
                                                                    axis.title.y = element_blank(), 
                                                                    axis.text.y = element_blank())
rplot_high


thetaplot_high<-tibble_otherpars_high %>% filter(group == "theta") %>% ggplot() + 
  geom_density(aes(x = sim1, y = after_stat(density)),alpha = 0.3, fill = pal[5], color = pal[5]) + 
  geom_density(aes(x = sim2, y = after_stat(density)),alpha = 0.3, fill = pal[5], color = pal[5]) + 
  geom_density(aes(x = sim3, y = after_stat(density)),alpha = 0.3, fill = pal[5], color = pal[5]) + 
  geom_density(aes(x = sim4, y = after_stat(density)),alpha = 0.3, fill = pal[5], color = pal[5]) + 
  geom_density(aes(x = sim5, y = after_stat(density)),alpha = 0.3, fill = pal[5], color = pal[5]) +
  geom_vline(aes(xintercept = 0.8), size = 0.8) +
  theme_classic() + labs(x = "Overdispersion (\u03B8)", y = "density") + 
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 4)) + theme(text = element_text(size = 14, family = "Arial"),
                                                                axis.text.y = element_blank(),
                                                                axis.ticks.y = element_blank(),
                                                                axis.title.y = element_blank())

thetaplot_high

Basinplot_high<-tibble_otherpars_high %>% filter(group %in% c("K1", "K2", "K3")) %>% ggplot() + 
  geom_density_ridges(aes(x = sim1, y = group, fill = group, color = group), rel_min_height = 0.02, alpha = 0.4,scale = 1) + 
  geom_density_ridges(aes(x = sim2, y = group, fill = group, color = group), rel_min_height = 0.02, alpha = 0.4,scale = 1) + 
  geom_density_ridges(aes(x = sim3, y = group,fill = group, color = group), rel_min_height = 0.02, alpha = 0.4, scale = 1) +  
  geom_density_ridges(aes(x = sim4, y = group, fill = group, color = group), rel_min_height = 0.02, alpha = 0.4, scale = 1) +  
  geom_density_ridges(aes(x = sim5, y = group, fill = group, color = group), rel_min_height = 0.02, alpha = 0.4, scale = 1) + 
  scale_y_discrete(breaks = c("K1", "K2", "K3"), labels = c("K1" = "Atl", "K2" = "Ind", "K3" = "Pac"), expand = c(0,0)) +
  #xlim() +
  scale_x_continuous(expand = c(0,0),limits = c(30000, 120000), breaks = c(40000, 60000, 80000,  100000), labels = c("40000","60000", "80000", "100000")) +
  scale_fill_manual(values = clrs) +
  scale_color_manual(values = clrs)+
  labs(x = "Carrying capacity (K)") + 
  annotate("text", x = c(95000, 95000, 95000), y = c("K3", "K2", "K1"), label = c("Pacific", "Indian", "Atlantic"), vjust = c(-1,-1,-1), family = "Arial", size = 4, color = c(clrs[3], clrs[2], clrs[1]), fontface = "bold") +
  geom_vline(aes(xintercept = 66667), size = 0.8) +
  theme_classic() +  
  theme(text = element_text(family = "Arial", size = 14),
        legend.position = "none", 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = rel(0.8))) 
Basinplot_high

highplot<-plot_sim_high_move |((rplot_high | tlplot_high)/(thetaplot_high | Basinplot_high)) +
  plot_layout(widths = c(1, 2))

highplot

#2) Constant Low Movement
y<-vector()
y1<-vector()
y2<-vector()
y3<-vector()
y4<-vector()

for(a in 1:9){
  y<-c(y,as.vector(as_draws_matrix(ConstLowMove$part_results$fits[[1]])[,a]))
  y1<-c(y1,as.vector(as_draws_matrix(ConstLowMove$part_results$fits[[2]])[,a]))
  y2<-c(y2,as.vector(as_draws_matrix(ConstLowMove$part_results$fits[[3]])[,a]))
  y3<-c(y3,as.vector(as_draws_matrix(ConstLowMove$part_results$fits[[4]])[,a]))
  y4<-c(y4,as.vector(as_draws_matrix(ConstLowMove$part_results$fits[[5]])[,a]))
}

#yrep<-matrix(data = c(y1,y2,y3,y4), nrow = 4, ncol = (600*9), byrow = TRUE)

grp<-rep(c("MAA", "MAB", "MAC", "MBA", "MBB", "MBC", "MCA", "MCB", "MCC"), each = (2000))

tibble_lowmove<-tibble(sim1 = y, sim2 = y1, sim3 = y2, sim4 = y3, sim5 = y4, group = grp)
tibble_lowmove_long<-tibble_lowmove %>% pivot_longer(c(sim1,sim2,sim3,sim4, sim5), names_to = "sim", values_to = "vals")

for(i in 1:5){
  print(mcmc_rhat_hist(ConstLowMove$part_results$rhats[[5]]))
} #all converged

for(i in 1:5){
  print(mcmc_trace(ConstLowMove$part_results$fits[[5]], facet_args = list(ncol = 3)))
} 


truth_tab_low<-tibble(grp = as.factor(c("MCC", "MCB", "MCA", "MBC", "MBB", "MBA", "MAC","MAB", "MAA")), 
                      truth = c(0.98, 0.01, 0.01, 0.01, 0.98, 0.01, 0.01, 0.01, 0.98), 
                      grpend = c("MCB", "MCA", "MBC", "MBB", "MBA", "MAC","MAB", "MAA","MCd"))

plot_sim_low_move<-ggplot(tibble_lowmove) + geom_density_ridges(aes(x = sim1, y = group), rel_min_height = 0.02, alpha = 0.3, scale = 1, fill = pal[1],color = pal[1]) +  
  geom_density_ridges(aes(x = sim2, y = group), rel_min_height = 0.02, alpha = 0.3,scale = 1, fill = pal[1],color = pal[1]) + 
  geom_density_ridges(aes(x = sim3, y = group), rel_min_height = 0.02, alpha = 0.3, scale = 1, fill = pal[1],color = pal[1]) +  
  geom_density_ridges(aes(x = sim4, y = group), rel_min_height = 0.02, alpha = 0.3, scale = 1, fill = pal[1],color = pal[1]) +  
  geom_density_ridges(aes(x = sim5, y = group), rel_min_height = 0.02, alpha = 0.3, scale = 1, fill = pal[1],color = pal[1]) + 
  scale_y_discrete(limits = c("MCC", "MCB", "MCA", "MBC", "MBB", "MBA", "MAC","MAB", "MAA", "MCd"), breaks = c("MAA", "MAB", "MAC", "MBA", "MBB", "MBC", "MCA", "MCB", "MCC"), labels = c("MAA" = "Atl_Atl", "MAB" = "Atl_Ind", "MAC" = "Atl_Pac", "MBA" = "Ind_Atl", "MBB" = "Ind_Ind", "MBC" = "Ind_Pac", "MCA" = "Pac_Atl", "MCB" = "Pac_Ind", "MCC" = "Pac_Pac"), expand = c(0,0.05)) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 1)) +
  labs(x = "Movement (m)", y = "density") + 
  annotate("text", x = c(0.7, 0.25, 0.25, 0.25, 0.7, 0.25, 0.25, 0.25, 0.7), 
           y = c("MAA", "MAB", "MAC", "MBA", "MBB", "MBC", "MCA", "MCB", "MCC"), 
           label = c("Atlantic \u2192 Atlantic", "Atlantic \u2192 Indian", "Atlantic \u2192 Pacific", 
                     "Indian \u2192 Atlantic", "Indian \u2192 Indian", "Indian\u2192 Pacific", 
                     "Pacific \u2192 Atlantic", "Pacific \u2192 Indian", "Pacific \u2192 Pacific"), size = 3.5, vjust = rep(-1.7, 9)) +
  geom_segment(data = truth_tab_low, aes(x = truth, xend = truth, y = grp, yend = grpend)) + 
  theme_classic() +  
  theme(text = element_text(family = "Arial", size = 14), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank()) +
  ggtitle("(b)  Low movement")
plot_sim_low_move

#other parameters
z<-vector()
z1<-vector()
z2<-vector()
z3<-vector()
z4<-vector()

for(a in 10:16){
  z<-c(z,as.vector(as_draws_matrix(ConstLowMove$part_results$fits[[1]])[,a]))
  z1<-c(z1,as.vector(as_draws_matrix(ConstLowMove$part_results$fits[[2]])[,a]))
  z2<-c(z2,as.vector(as_draws_matrix(ConstLowMove$part_results$fits[[3]])[,a]))
  z3<-c(z3,as.vector(as_draws_matrix(ConstLowMove$part_results$fits[[4]])[,a]))
  z4<-c(z4,as.vector(as_draws_matrix(ConstLowMove$part_results$fits[[5]])[,a]))
}

grp<-rep(c("tl", "r", "lnK", "K1", "K2", "K3", "theta"), each = (2000))

tibble_otherpars_low<-tibble(sim1 = z, sim2 = z1, sim3 = z2, sim4 = z3, sim5 = z4, group = grp)
tibble_otherpars_low_long<-tibble_otherpars_low %>% pivot_longer(c(sim1,sim2,sim3,sim4, sim5), names_to = "sim", values_to = "vals")

tlplot_low<-tibble_otherpars_low %>% filter(group == "tl") %>% ggplot() + 
  geom_density(aes(x = sim1, y = after_stat(density)),alpha = 0.3, fill = pal[2], color = pal[2]) + 
  geom_density(aes(x = sim2, y = after_stat(density)),alpha = 0.3, fill = pal[2], color = pal[2]) + 
  geom_density(aes(x = sim3, y = after_stat(density)),alpha = 0.3, fill = pal[2], color = pal[2]) + 
  geom_density(aes(x = sim4, y = after_stat(density)),alpha = 0.3, fill = pal[2], color = pal[2]) + 
  geom_density(aes(x = sim5, y = after_stat(density)),alpha = 0.3, fill = pal[2], color = pal[2]) +
  geom_vline(aes(xintercept = 0.96), size = 0.8) +
  theme_classic() + labs(x = "Mark loss (l)", y = "density") + 
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 1.05)) + theme(text = element_text(size = 14, family = "Arial"),
                                                                   axis.text.y = element_blank(), 
                                                                   axis.title.y = element_blank(), 
                                                                   axis.ticks.y = element_blank())
tlplot_low
rplot_low<-tibble_otherpars_low %>% filter(group == "r") %>% ggplot() + 
  geom_density(aes(x = sim1, y = after_stat(density)),alpha = 0.3, fill = pal[4], color = pal[4]) + 
  geom_density(aes(x = sim2, y = after_stat(density)),alpha = 0.3, fill = pal[4], color = pal[4]) + 
  geom_density(aes(x = sim3, y = after_stat(density)),alpha = 0.3, fill = pal[4], color = pal[4]) + 
  geom_density(aes(x = sim4, y = after_stat(density)),alpha = 0.3, fill = pal[4], color = pal[4]) + 
  geom_density(aes(x = sim5, y = after_stat(density)),alpha = 0.3, fill = pal[4], color = pal[4]) +
  geom_vline(aes(xintercept = 0.073), size = 0.8) +
  theme_classic() + labs(x = "Intrinsic growth (r)", y = "density") + 
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 0.118)) + theme(text = element_text(size = 14, family = "Arial"),
                                                                    axis.ticks.y = element_blank(),
                                                                    axis.title.y = element_blank(), 
                                                                    axis.text.y = element_blank())
rplot_low


thetaplot_low<-tibble_otherpars_low %>% filter(group == "theta") %>% ggplot() + 
  geom_density(aes(x = sim1, y = after_stat(density)),alpha = 0.3, fill = pal[5], color = pal[5]) + 
  geom_density(aes(x = sim2, y = after_stat(density)),alpha = 0.3, fill = pal[5], color = pal[5]) + 
  geom_density(aes(x = sim3, y = after_stat(density)),alpha = 0.3, fill = pal[5], color = pal[5]) + 
  geom_density(aes(x = sim4, y = after_stat(density)),alpha = 0.3, fill = pal[5], color = pal[5]) + 
  geom_density(aes(x = sim5, y = after_stat(density)),alpha = 0.3, fill = pal[5], color = pal[5]) +
  geom_vline(aes(xintercept = 0.8), size = 0.8) +
  theme_classic() + labs(x = "Overdispersion (\u03B8)", y = "density") + 
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 4)) + theme(text = element_text(size = 14, family = "Arial"),
                                                                axis.text.y = element_blank(),
                                                                axis.ticks.y = element_blank(),
                                                                axis.title.y = element_blank())

thetaplot_low

Basinplot_low<-tibble_otherpars_low %>% filter(group %in% c("K1", "K2", "K3")) %>% ggplot() + 
  geom_density_ridges(aes(x = sim1, y = group, fill = group, color = group), rel_min_height = 0.02, alpha = 0.4,scale = 1) + 
  geom_density_ridges(aes(x = sim2, y = group, fill = group, color = group), rel_min_height = 0.02, alpha = 0.4,scale = 1) + 
  geom_density_ridges(aes(x = sim3, y = group,fill = group, color = group), rel_min_height = 0.02, alpha = 0.4, scale = 1) +  
  geom_density_ridges(aes(x = sim4, y = group, fill = group, color = group), rel_min_height = 0.02, alpha = 0.4, scale = 1) +  
  geom_density_ridges(aes(x = sim5, y = group, fill = group, color = group), rel_min_height = 0.02, alpha = 0.4, scale = 1) + 
  scale_y_discrete(breaks = c("K1", "K2", "K3"), labels = c("K1" = "Atl", "K2" = "Ind", "K3" = "Pac"), expand = c(0,0)) +
  scale_x_continuous(limits = c(30000, 120000),expand = c(0,0), breaks = c(40000, 60000, 80000,  100000), labels = c("40000", "60000","80000", "100000")) +
  scale_fill_manual(values = clrs) +
  scale_color_manual(values = clrs)+
  labs(x = "Carrying capacity (K)") + 
  geom_text( data = tibble(rep(NA, 3)), aes(x = c(100000, 96000,92000), y = c("K3", "K2", "K1"), label = c("Pacific", "Indian", "Atlantic")), vjust = c(-1, -1,-1), family = "Arial", size = 4, color = c(clrs[3], clrs[2], clrs[1]), fontface = "bold") +
  geom_vline(aes(xintercept = 66667), size = 0.8) +
  theme_classic() +  
  theme(text = element_text(family = "Arial", size = 14),
        legend.position = "none", 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.x = element_text(size = rel(0.8))) 

Basinplot_low

lowplot<-plot_sim_low_move |((rplot_low | tlplot_low)/(thetaplot_low | Basinplot_low)) + plot_layout(widths = c(1, 2))
lowplot
#combining into one plot
comb<-highplot/lowplot 
comb
#comb2<-plot_grid(highplot, lowplot, nrow = 2)
#comb2

comb3<-(plot_sim_low_move |((rplot_low | tlplot_low)/(thetaplot_low | Basinplot_low)))/(plot_sim_high_move |((rplot_high | tlplot_high)/(thetaplot_high | Basinplot_high))) + plot_layout(widths = c(1, 3,3))
comb3

plot_sim_low_move + rplot_low + tlplot_low + thetaplot_low + Basinplot_low + plot_sim_high_move + rplot_high + tlplot_high + thetaplot_high + Basinplot_high + plot_layout(ncol = 3, nrow = 4, widths = c(1, 2, 2))

design<-"
AAAAABBBCCC
AAAAADDDEEE
FFFFFGGGHHH
FFFFFKKKLLL

"
comb2<-plot_sim_high_move + rplot_high + tlplot_high + thetaplot_high + Basinplot_high + plot_sim_low_move + rplot_low + tlplot_low + thetaplot_low + Basinplot_low +  plot_layout(design = design)
comb2
ggsave(comb2, file = "Figures/SimComb.png", dpi = 600, height = 10, width = 8, units = "in")


# Time Varying Results ----------------------------------------------------
#posterior predictive for movement

bckgrd<-"#FFFFFF"
#excluded areas with no recoveries
Rec_long<-Rec_long %>%
  filter(group %in% tag_groups[c(1:2,4:6,8:9)])
#adding marking basin
Rec_long<-Rec_long %>% mutate(mark_basin = if_else("A" == str_sub(group, 1,1), "Atlantic", if_else("B" == str_sub(group, 1,1), "Indian", "Pacific")))

#plot of recovery data
group_names<-as_labeller(c("AA" = "Atl_Atl", "AB" = "Atl_Ind", "AC" = "Atl_Pac", "BA" = "Ind_Atl", "BB" = "Ind_Ind", "BC" = "Ind_Pac", "CA" = "Pac_Atl","CB" = "Pac_Ind", "CC" = "Pac_Pac"))


f1<-ggplot() + geom_bar(data = Rec_long, aes(x = Yr_rec, y = recoveries, fill = mark_basin), stat = "identity" , alpha = 0.8) + 
  scale_x_continuous(expand = c(0,0)) +
  labs(x = "Season Recovered", y = "Number of Marks Recovered") + scale_fill_manual(values = clrs) + theme_light() + theme(legend.position = "none",                                                                                                                        axis.title = element_text(colour = "#256482"), 
                                                                                                                           panel.border = element_blank())
f1

#adding posterior predictive fit
nyears<-46
ThinnedPPRec<-subset_draws(tv_StanThinned, variable = "PostPredRec")
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
PPRecStan$Year<-rep(Years[15:60], 9)
head(PPRecStan)

PPRecStan<-PPRecStan %>% rename(Lower = "2.5%", Upper = "97.5%")
#not including 0s
NTRStan_filt<-PPRecStan #%>%
#filter(q95 >=0.10)

scales_y2<-list(
  AA = scale_y_continuous(limits = c(0,6), breaks = seq(0,6, by = 1), expand = c(0,0), minor_breaks = 0),
  AB = scale_y_continuous(limits = c(0,3.7), breaks = seq(0,3, by = 1), expand = c(0,0),minor_breaks = 0),
  AC = scale_y_continuous(limits = c(0,3.7),breaks = seq(0,3, by = 1), expand = c(0,0),minor_breaks = 0),
  BA = scale_y_continuous(limits = c(0,3.7),breaks = seq(0,3, by = 1), expand = c(0,0),minor_breaks = 0),
  BB = scale_y_continuous(limits = c(0,5.2),breaks = seq(0,5, by = 1), expand = c(0,0),minor_breaks = 0),
  BC = scale_y_continuous(limits = c(0,3.7),breaks = seq(0,3, by = 1), expand = c(0,0),minor_breaks = 0),
  CA = scale_y_continuous(limits = c(0,3.7), breaks = seq(0,3, by = 1), expand = c(0,0),minor_breaks = 0),
  CB = scale_y_continuous(limits = c(0,3.7),breaks = seq(0,3, by = 1), expand = c(0,0),minor_breaks = 0),
  CC = scale_y_continuous(limits = c(0,3.7),breaks = seq(0,3, by = 1), expand = c(0,0),minor_breaks = 0)
)

f2<-f1 + geom_point(data = NTRStan_filt, aes(x = Year, y = median), size = 0.8) + 
  geom_linerange(data = NTRStan_filt, aes(x = Year, ymin = Lower, ymax = Upper), size = 0.5) +
  facet_grid_sc(rows = vars(group), scales = list(y = scales_y2), space = "free_y", switch = "y", labeller = group_names) + 
  theme(strip.background = element_blank(), strip.placement = "outside", panel.spacing = unit(0.5, "lines"))
f2

ggsave(f2, "Custom_Model/Figures/timevaryingMovementFit.png", dpi = 600)

#posterior predictive for population 

bckgrd<-"#FFFFFF"
txt<-"#000000"
library(grid)
library(gtable)
ThinnedQ<-subset_draws(tv_draws_array, variable = "q")
ThinnedPPBr<-subset_draws(tv_draws_array, variable = "PostPredBr")
ThinnedPPMatInd<-subset_draws(tv_draws_array, variable = "PostPredMatInd")
ThinnedPPMatPac<-subset_draws(tv_draws_array, variable = "PostPredMatPac")
#PPMatIndQ<-rep(ThinnedQ[,, variable = "q[1]"], each = 10) * ThinnedPPMatInd #don't need this because automatically accounted for in rng
#PPMatPacQ<-rep(ThinnedQ[,, variable = "q[2]"], each = 10) * ThinnedPPMatPac #don't need this because automatically accounted for in rng
Quants_PPBr<-summarise_draws(ThinnedPPBr, "median", ~quantile(.x, probs = c(0.025, 0.975)))[c("median", "2.5%", "97.5%")]
QuantsMIQ<-summarise_draws(ThinnedPPMatInd, "median", ~quantile(.x, probs = c(0.025, 0.975)))[c("median", "2.5%", "97.5%")]
QuantsMPQ<-summarise_draws(ThinnedPPMatPac, "median", ~quantile(.x, probs = c(0.025, 0.975)))[c("median", "2.5%", "97.5%")]


Quants_PPBr_violin<-ThinnedPPBr %>% as_draws_matrix %>% as_tibble %>% pivot_longer(everything(), names_to = "Par", values_to = "PopSize") %>% add_column(Year = rep(Abund_dat$Year[1:9], 3000), Basin = rep(Abund_dat$Basin[1:9], 3000))
Quants_MIQ_violin<-ThinnedPPMatInd %>% as_draws_matrix %>% as_tibble %>% pivot_longer(everything(), names_to = "Par", values_to = "PopSize") %>% add_column(Year = rep(Abund_dat$Year[10:19], 3000), Basin = rep(Abund_dat$Basin[10:19], 3000))
Quants_MPQ_violin<-ThinnedPPMatPac %>% as_draws_matrix %>% as_tibble %>% pivot_longer(everything(), names_to = "Par", values_to = "PopSize") %>% add_column(Year = rep(Abund_dat$Year[20:29], 3000), Basin = rep(Abund_dat$Basin[20:29], 3000))


PopQuants<-bind_rows(Quants_PPBr, QuantsMIQ, QuantsMPQ)
PopQuants<-PopQuants %>% add_column(Year = Abund_dat$Year, Basin = Abund_dat$Basin, Data = Abund_dat$PopSize, 
                                    Source = Abund_dat$Source) %>%
  rename(Lower = "2.5%", Upper = "97.5%")

PopQuantsViolin<-bind_rows(Quants_PPBr_violin, Quants_MIQ_violin, Quants_MPQ_violin) %>% add_column(Source = c(rep(Abund_dat$Source[1:9], 3000), rep(Abund_dat$Source[10:19], 6000)))

source.labs<-c("Branch 2007", "Matsuoka & Hakamada 2014")
names(source.labs)<-c("Branch", "Matsouka")

#As Violin Plot: 
PPAbund<-ggplot() + geom_violin(data = PopQuantsViolin, aes(x = as.factor(Year), y = PopSize), alpha = 0.5) + 
  geom_point(data = Abund_dat, aes(x = as.factor(Year), y = PopSize, color = Basin)) + facet_grid(Source ~ Basin, scales = "free_x", labeller = labeller(Source = source.labs)) + 
  scale_color_manual(values = clrs) + #scale_fill_manual(values = clrs)+
  #scale_x_discrete(breaks = seq(1981, 2008, by = 5)) + 
  scale_y_continuous(limits = c(0, 1200), breaks = seq(0, 1200, by  = 200)) + 
  labs(y = "Population Size", x = "Year", tag = "(b)")  + 
  theme_bw() + theme(legend.position = "none",
                     legend.background = element_rect(fill = bckgrd),
                     panel.background = element_rect(fill = bckgrd), 
                     panel.border = element_rect(fill = "transparent"),
                     #panel.grid = element_line(color = txt),
                     plot.background = element_rect(fill = bckgrd),
                     axis.text = element_text(colour = txt, size = rel(0.8)), 
                     axis.text.x = element_text(angle =90), 
                     #axis.ticks = element_line(colour = txt), 
                     #axis.line = element_line(colour = txt, size = rel(1)), 
                     #axis.title = element_text(colour = txt),
                     strip.background = element_blank(), 
                     plot.tag.position = "topleft")


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
#ggsave(PPAbund_violin, "Custom_Model/Figures/PPAbund_violin.png", dpi = 600)


#population trajectory
bckgrd<-"#FFFFFF"
txt<-"#000000"
nyears<-length(Years)


ThinnedNPop<-subset_draws(tv_StanThinned, variable = "Npop")
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
  #geom_point(data = Abund_dat_BR, aes(x = Year, y = PopSize, color = Basin)) + 
  labs(y = "Population size")+ 
  scale_y_continuous(labels = comma, breaks = c(0,5000, 10000,50000,75000,100000), expand = c(0,0))  + 
  scale_x_continuous(expand = c(0,0), breaks = c(seq(1910, 2000, by =10 ))) + expand_limits(y = c(0,100000), x = c(1913, 2001)) +scale_color_manual(values = clrs) + scale_fill_manual(values = clrs) +  
  theme_classic() + theme(legend.background = element_rect(fill = bckgrd),
                          panel.background = element_rect(fill = bckgrd), 
                          plot.background = element_rect(fill = bckgrd),
                          axis.text = element_text(colour = txt, size = rel(0.8)), 
                          #axis.ticks = element_line(colour = "#256482"), 
                          #axis.line = element_line(colour = "#256482", size = rel(1)), 
                          axis.title = element_text(colour = txt), 
                          panel.border = element_blank(), 
                          legend.position = "bottom") 
popplot

popplotlog<-popplot + scale_y_continuous(trans = "log10", breaks = c(50,1000, 10000,50000,100000), labels = comma) +
  theme(legend.position = "bottom")
popplotlog

popplot2<- popplotlog + ggtitle("Posterior distribution of population size")
popplot2 



PPAbund_violin2<-PPAbund_violin+ labs(tag = "(b)") + theme(plot.tag.position = "topright", 
                                                           panel.border = element_blank())
PPAbund_violin2

popplot2 + PPAbund_violin

#library(cowplot)
abundfig<-plot_grid(popplot2, PPAbund_violin, rel_widths = c(1,2))
abundfig

#ggsave(abundfig, "Custom_Model/Figures/abundancemodelfit.png", dpi = 600)


#Carrying Capacity and 3 tag loss parameters posterior distributions
#Total Carrying Capacity (K)
prior_lnK<-matrix(data = runif(1000,9,13), nrow = 1000, ncol = 1)
colnames(prior_lnK)<-"lnK"
priorK<-as_tibble(prior_lnK)
priorK$y<-dunif(priorK$lnK, min = 9, max = 13)

lnkpost<-subset_draws(tv_draws_array, "lnK")
lnKub<-summarise_draws(lnkpost, ~quantile(.x, probs = c(0.025, 0.975)))[,3] %>% rename("q975" = "97.5%")
lnKlb<-summarise_draws(lnkpost, ~quantile(.x, probs = c(0.025, 0.975)))[,2] %>% rename("q025" = "2.5%")
lnkmed<-summarise_draws(lnkpost, "median")
lnkposttab<-lnkpost %>% as_draws_matrix() %>% as_tibble()

kpost<-ggplot() + geom_density(data = lnkposttab, aes(x = lnK, y = after_stat(density)),alpha = 1) + 
  theme_classic() + labs(x = "Carrying capacity (K)", y = "density", tag = "(a)") + 
  geom_density(data = priorK, aes(x = lnK, y = y+5), fill = NA, linetype = "dashed") + 
  scale_x_continuous(expand = c(0,0), limits  = c(12,12.5), breaks = log(c(170000,190000,210000, 230000)), 
                     labels = c("170000","190000","210000", "230000")) + 
  scale_y_continuous(expand = c(0,0)) #+
#coord_fixed(ratio = 0.04,xlim = c(12,12.4), ylim = c(0,11))

kpost

kpost1<-kpost %>% plot_credible_interval(lnKlb$q025, lnKub$q975, fillcol = pal2[1])
kpost1

#r
prior_r<-matrix(data = runif(1000,0,0.114), nrow = 1000, ncol = 1)
colnames(prior_r)<-"r"
priorr<-as_tibble(prior_r)
priorr$y<-dunif(priorr$r, min = 0, max = 0.114)

r_post<-subset_draws(tv_draws_array, "r")
rub<-summarise_draws(r_post,~quantile(.x, probs = c(0.025, 0.975)))[,3] %>% rename("q975" = "97.5%")
rlb<-summarise_draws(r_post,~quantile(.x, probs = c(0.025, 0.975)))[,2] %>% rename("q025" = "2.5%")
rmed<-summarise_draws(r_post, "median")
rposttab<-r_post %>% as_draws_matrix() %>% as_tibble()
rpost<-ggplot() + geom_density(data = rposttab, aes(x = r, y = after_stat(density)),alpha = 1) + 
  theme_classic() + labs(x = "r", y = "density", tag = "(b)") + 
  geom_density(data = priorr, aes(x = r, y = y+5), fill = NA, linetype = "dashed") + 
  #geom_vline(aes(xintercept = rmed$median), color = "#3cadbf", size = 2)+
  scale_x_continuous(expand = c(0,0), limits  = c(0,0.114), breaks = c(0.01, 0.03, 0.05, 0.07, 0.09, 0.11)) + 
  scale_y_continuous(expand = c(0,0)) + theme(panel.background = element_rect(fill = bckgrd), 
                                              plot.background = element_rect(fill = bckgrd),
                                              #axis.text = element_text(colour = "#0f3040", size = rel(0.8)), 
                                              #axis.ticks = element_line(colour = "#0f3040"), 
                                              #axis.line = element_line(colour = "#0f3040", size = rel(1)), 
                                              #axis.title = element_text(colour = "#0f3040"), 
                                              panel.border = element_blank())


rpost

rpost1<-rpost %>% plot_credible_interval(rlb$q025, rub$q975, fillcol = pal2[2])
rpost1
#tag loss
prior_tl<-matrix(data = rbeta(1000, 1,1), nrow = 1000, ncol = 1)
colnames(prior_tl)<-"tl"
priortl<-as_tibble(prior_tl)
priortl$yval<-dbeta(priortl$tl, 1, 1)

#tl1930
tl1post<-subset_draws(tv_draws_array, "tl1930")
tl1ub<-summarise_draws(tl1post, ~quantile(.x, probs = c(0.025, 0.975)))[,3] %>% rename("q975" = "97.5%")
tl1lb<-summarise_draws(tl1post, ~quantile(.x, probs = c(0.025, 0.975)))[,2] %>% rename("q025" = "2.5%")
tl1med<-summarise_draws(tl1post, "median")
tl1posttab<-tl1post %>% as_draws_matrix() %>% as_tibble()

tl1postplot<-ggplot() + geom_density(data = tl1posttab, aes(x = tl1930, y = after_stat(density)),alpha = 1) + 
  theme_classic() + labs(x = "Probability of mark loss (l) 1930 - 1943 ", y = "density", tag = "(c)")+ 
  geom_density(data = priortl, aes(x = tl, y = yval +5), fill = NA, linetype = "dashed") + scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0), limits = c(0.9, 1))

tl1postplot

tl1postplot1<-tl1postplot %>% plot_credible_interval(tl1lb$q025, tl1ub$q975, fillcol = pal2[3])
tl1postplot1

#tl1944
tl2post<-subset_draws(tv_draws_array, "tl1944")
tl2ub<-summarise_draws(tl2post, ~quantile(.x, probs = c(0.025, 0.975)))[,3] %>% rename("q975" = "97.5%")
tl2lb<-summarise_draws(tl2post, ~quantile(.x, probs = c(0.025, 0.975)))[,2] %>% rename("q025" = "2.5%")
tl2med<-summarise_draws(tl2post, "median")
tl2posttab<-tl2post %>% as_draws_matrix() %>% as_tibble()

tl2postplot<-ggplot() + geom_density(data = tl2posttab, aes(x = tl1944, y = after_stat(density)),alpha = 1) + 
  theme_classic() + labs(x = "Probability of mark loss (l) 1944 - 1956 ", y = "density", tag = "(d)")+ 
  geom_density(data = priortl, aes(x = tl, y = yval +5), fill = NA, linetype = "dashed") + scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0), limits = c(0.5, 1))

tl2postplot

tl2postplot1<-tl2postplot %>% plot_credible_interval(tl2lb$q025, tl2ub$q975, fillcol = pal2[4])
tl2postplot1

#tl1957
tl3post<-subset_draws(tv_draws_array, "tl1957")
tl3ub<-summarise_draws(tl3post, ~quantile(.x, probs = c(0.025, 0.975)))[,3] %>% rename("q975" = "97.5%")
tl3lb<-summarise_draws(tl3post, ~quantile(.x, probs = c(0.025, 0.975)))[,2] %>% rename("q025" = "2.5%")
tl3med<-summarise_draws(tl3post, "median")
tl3posttab<-tl3post %>% as_draws_matrix() %>% as_tibble()

tl3postplot<-ggplot() + geom_density(data = tl3posttab, aes(x = tl1957, y = after_stat(density)),alpha = 1) + 
  theme_classic() + labs(x = "Probability of mark loss (l) 1957+ ", y = "density", tag = "(e)")+ 
  geom_density(data = priortl, aes(x = tl, y = yval +5), fill = NA, linetype = "dashed") + scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0), limits = c(0.6, 1))

tl3postplot

tl3postplot1<-tl3postplot %>% plot_credible_interval(tl3lb$q025, tl3ub$q975, fillcol = pal2[5])
tl3postplot1


combKandtl<-(kpost1 + rpost1)/(tl1postplot1 + tl2postplot1 + tl3postplot1)
combKandtl


#movement rates
bckgrd<-"#FFFFFF"
#custom_scheme <- rev(c("#0f3040", "#0f3040",
#"#eef2e8", "#eef2e8",
#"#0092a8","#0092a8"))
#color_scheme_set("viridisA")
#str(color_scheme_get())

priors_movemat<-matrix(data = NA, nrow = 1000, ncol = 6)
for(i in 1:ncol(priors_movemat)){
  priors_movemat[,i]<-runif(1000, min = 0, max = 0.34)
}
priors_movemat_y<-matrix(data = NA, nrow = 1000, ncol = 6)
priors_movemat_y<-dunif(priors_movemat, min = 0, max = 0.34)

colnames(priors_movemat)<-c("MAB", "MAC", "MBA", "MBC", "MCA", "MCB")
colnames(priors_movemat_y)<-c("MAB", "MAC", "MBA", "MBC", "MCA", "MCB")

priors<-as_tibble(priors_movemat) %>% pivot_longer(everything())
priorsy<-as_tibble(priors_movemat_y) %>% pivot_longer(everything()) %>% rename( yval = value)

priors<-priors %>% left_join(priorsy, by = "name")

color_scheme_set("teal")
m1<-mcmc_areas(
  tv_draws_array, 
  pars = c("MAA","MAB", "MAC", "MBA", "MBB", "MBC", "MCA", "MCB", "MCC"),
  prob = 0.95, 
  point_est = "median",
  area_method = "scaled height") + geom_density_ridges(data = priors, aes(x = value, y = name, height = yval, group = name), scale = 0.5, fill = NA, color = "gray40", linetype = "dashed", rel_min_height = 0.01) + 
  scale_x_continuous(expand = c(0,0)) + scale_y_discrete(labels = c("MAA" = "Atl_Atl", "MAB" = "Atl_Ind", "MAC" = "Atl_Pac", "MBA" = "Ind_Atl", "MBB" = "Ind_Ind", "MBC" = "Ind_Pac", "MCA" = "Pac_Atl", "MCB" = "Pac_Ind", "MCC" = "Pac_Pac"), expand = c(0,0)) + 
  labs(y = "Movement Parameter (from_to)", x = "Probability") + theme(text = element_text(size = 15))

m1

#comparing to movement without time varying results
color_scheme_set("green")
m2<-mcmc_areas(draws_array, 
               pars = c("MAA","MAB", "MAC", "MBA", "MBB", "MBC", "MCA", "MCB", "MCC"),
               prob = 0.95, 
               point_est = "median",
               area_method = "scaled height") + geom_density_ridges(data = priors, aes(x = value, y = name, height = yval, group = name), scale = 0.5, fill = NA, color = "gray40", linetype = "dashed", rel_min_height = 0.01) + 
  scale_x_continuous(expand = c(0,0)) + scale_y_discrete(labels = c("MAA" = "Atl_Atl", "MAB" = "Atl_Ind", "MAC" = "Atl_Pac", "MBA" = "Ind_Atl", "MBB" = "Ind_Ind", "MBC" = "Ind_Pac", "MCA" = "Pac_Atl", "MCB" = "Pac_Ind", "MCC" = "Pac_Pac"), expand = c(0,0)) + 
  labs(y = "Movement Parameter (from_to)", x = "Probability") + theme(text = element_text(size = 15))
m2

m1+m2



combined <- rbind(mcmc_areas_data(draws_mat, pars = c("MAA","MAB", "MAC", "MBA", "MBB", "MBC", "MCA", "MCB", "MCC"), prob = 0.95), mcmc_areas_data(tv_draws_mat, pars = c("MAA","MAB", "MAC", "MBA", "MBB", "MBC", "MCA", "MCB", "MCC"), prob = 0.95))
combined$model <- rep(c("Constant", "Time Varying"), each = 18745)
combined$name <-paste0(combined$parameter, combined$model)

#pos <- position_nudge(y = ifelse(combined$model == "Time Varying", 0, 0.1))
#combined %>% 
ggplot(data = combined, aes(x = x, height = plotting_density, y = parameter, color = model, fill = model, group = name)) + 
  geom_density_ridges(scale = 1,  rel_min_height = 0.01, stat = "identity", alpha = 0.8) + 
  scale_x_continuous(expand = c(0,0)) + scale_y_discrete(labels = c("MAA" = "Atl_Atl", "MAB" = "Atl_Ind", "MAC" = "Atl_Pac", "MBA" = "Ind_Atl", "MBB" = "Ind_Ind", "MBC" = "Ind_Pac", "MCA" = "Pac_Atl", "MCB" = "Pac_Ind", "MCC" = "Pac_Pac"), expand = c(0,0)) +
  theme_ridges(grid = FALSE, center_axis_labels = TRUE) + labs(y = "Movement Parameter (from_to)", x = "Probability")



#geom_linerange(aes(xmin = l, xmax = h), position = pos, size=2)+
#geom_linerange(aes(xmin = ll, xmax = hh), position = pos)+
#geom_point(position = pos, color="black")
#summary results
View(head(summary(tv_draws_array), 18))
