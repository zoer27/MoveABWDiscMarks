
# Packages ----------------------------------------------------------------
library(tidyverse)
# Functions ---------------------------------------------------------------

source("Functions/Data_Manip_Functions.R")
source("Functions/NewGroupsFunctions.R")

# Mark Recovery Data  --------------------------------------------------------------------
ABW_marks_recaps<-read_csv("Data/ABW_marks_recaps_diffseasons.csv") #just contains ABW marks with recaptures and recapture information, only different seasons 
ABW_marks_norecaps<-read_csv("Data/ABW_marks_norecaps_clean.csv") #contains marks that were placed but never recaptured


# Catch Data --------------------------------------------------------------
ABW_Catches<-read_csv("Data/ABW_catch_clean.csv")
TB_2004<-read_csv("Data/TB2004_Catches.csv")
TB_2004$Year<-TB_2004$Year - 1 #Trevor used end of whaling season in his paper, but I am using start of whaling season


# Abundance Data ----------------------------------------------------------
TB_2007<-read_csv("Data/TB_2007_Abundance.csv")[1:19, 1:5]
#Mat_Abund_2014<-read_csv("Custom_Model/Data/Matsouka_Abundance.csv")
#excluding years where only areas IV and V were surveyed
#Mat_Abund_filtered<-Mat_Abund_2014 %>% filter(Season >= 1995)
#converting areas into numbers
#Mat_area_nums<-c("IV" = '4', "IIIE" = "3", "V" = "5", "VIW" = "6")
#Mat_Abund_filtered$AreaNum<-recode(Mat_Abund_filtered$Area, !!!Mat_area_nums)
HambeAbund<-read_csv("Data/Hambe2023Abundance.csv")
Ham_area_nums<-c("IV" = '4',  "V" = "5")
HambeAbund$AreaNum<-recode(HambeAbund$Area, !!!Ham_area_nums)


# Groups ------------------------------------------------------------------
#I/II vs III/IV vs V/VI
opt1<-c("1" = "A", "2" = "A", "3" = "B", "4" = "B", "5" = "C", "6" = "C")
#A is 120W to 0
#B is 0 to 130E
#C is 130E to 20W
opt1brks<-c(-120, 0, 130)

#II/III, vs IV/V, vs VI/I
opt2<-c("1" = "C", "2" = "A", "3" = "A", "4" = "B", "5" = "B", "6" = "C")
#A is 60W to 70E
#B is 70E to 170W
#C is 170W to 60W
opt2brks<-c(-60, 70, -170)

# Mark Recovery Data Summary and Cleaning ---------------------------------

#summarizing mark data
#last season is 1972
marking_years<-unique(c(unique(ABW_marks_norecaps$Seasn_mark), unique(ABW_marks_recaps$Seasn_mark + 1900)))
all_years<-seq(1926, 1972, by = 1)

#####Option1 
#tibble for releases
releases_opt1<-releases_tab_newgrps(opt1)
head(releases_opt1)
summary(releases_opt1)



#####Option 2 
#tibble for releases
releases_opt2<-releases_tab_newgrps(opt2)
head(releases_opt2)
summary(releases_opt2)



#recoveries, grouped by where they were released
recoveries_opt1<-recoveries_tab_newgrps(opt1)
head(recoveries_opt1)
summary(recoveries_opt1)
#sum(recoveries_opt1[,-1])

#plotting
tag_groups<-c("AA", "AB", "AC", "BA", "BB","BC", "CA", "CB","CC")
Rec_long<-recoveries_opt1 %>% pivot_longer(all_of(tag_groups), names_to = "group", values_to = "recoveries") %>% 
  mutate(mark_basin = if_else("A" == str_sub(group, 1,1), "I&II", if_else("B" == str_sub(group, 1,1), "III&IV", "V&VI")))
ggplot() + geom_bar(data = Rec_long, aes(x = Yr_rec, y = recoveries, fill = mark_basin), stat = "identity" , alpha = 0.8) + 
  facet_grid(rows = vars(group))


#Option 2
recoveries_opt2<-recoveries_tab_newgrps(opt2)
head(recoveries_opt2)
summary(recoveries_opt2)
Rec_long<-recoveries_opt2 %>% pivot_longer(all_of(tag_groups), names_to = "group", values_to = "recoveries") %>% 
  mutate(mark_basin = if_else("A" == str_sub(group, 1,1), "II&III", if_else("B" == str_sub(group, 1,1), "IV&V", "VI&I")))
ggplot() + geom_bar(data = Rec_long, aes(x = Yr_rec, y = recoveries, fill = mark_basin), stat = "identity" , alpha = 0.8) + 
  facet_grid(rows = vars(group))

#save
#write_csv(recoveries_opt1, file = "Custom_Model/Data/mark_recoveries_opt1.csv")
#write_csv(releases_opt1, file = "Custom_Model/Data/mark_releases_opt1.csv")

#write_csv(recoveries_opt2, file = "Custom_Model/Data/mark_recoveries_opt2.csv")
#write_csv(releases_opt2, file = "Custom_Model/Data/mark_releases_opt2.csv")

# Catches ----------------------------------------------------------
catches_opt1<-cor_catches_by_group(opt1brks, 1)
head(catches_opt1)
summary(catches_opt1)

#plotting
Catch_long<-catches_opt1 %>% pivot_longer(c("A", "B", "C"), names_to = "basin", values_to = "catch")
basin.labs1<-c("I&II", "III&IV", "V&VI")
names(basin.labs1)<-c("A", "B", "C")
ggplot(data = Catch_long) + 
  geom_bar(aes(x = Year, y = catch, fill = basin),stat = "identity") + 
  facet_grid(rows = vars(basin), space = "free_y", scales = "free_y", labeller = labeller(basin = basin.labs1)) +
  scale_x_continuous(breaks = seq(1904, 1973, by = 4), 
                     labels = as.character(seq(1904, 1973, by = 4)), 
                     expand = c(0,0)) +
  labs(y = "Catches", x = "Season")



catches_opt2<-cor_catches_by_group(opt2brks, 2)
head(catches_opt2)
summary(catches_opt2)

#plotting
Catch_long<-catches_opt2 %>% pivot_longer(c("A", "B", "C"), names_to = "basin", values_to = "catch")
basin.labs2<-c("II&III", "IV&V", "VI&I")
names(basin.labs2)<-c("A", "B", "C")
ggplot(data = Catch_long) + 
  geom_bar(aes(x = Year, y = catch, fill = basin),stat = "identity") + 
  facet_grid(rows = vars(basin), space = "free_y", scales = "free_y", labeller = labeller(basin = basin.labs2)) +
  scale_x_continuous(breaks = seq(1904, 1973, by = 4), 
                     labels = as.character(seq(1904, 1973, by = 4)), 
                     expand = c(0,0)) +
  labs(y = "Catches", x = "Season")


#write_csv(catches_opt1, file = "Custom_Model/Data/catchcorrected_opt1.csv")

#write_csv(catches_opt2, file = "Custom_Model/Data/catchcorrected_opt2.csv")



# Checking that years for Catches and Recoveries make sense ---------------
zero_catches<-catches_opt1 %>% filter(A == 0 | B == 0 | C == 0) %>% filter(Year >=1926)
yrs_with_recs<-recoveries_opt1 %>% filter(if_any(c(-Yr_rec), ~ .x >0)) 
zero_catches[which(zero_catches$Year %in% yrs_with_recs$Yr_rec ),]
yrs_with_recs[which(yrs_with_recs$Yr_rec %in% zero_catches$Year),]


zero_catches<-catches_opt2 %>% filter(A == 0 | B == 0 | C == 0) %>% filter(Year >=1926)
yrs_with_recs<-recoveries_opt2 %>% filter(if_any(c(-Yr_rec), ~ .x >0)) 
zero_catches[which(zero_catches$Year %in% yrs_with_recs$Yr_rec ),]
yrs_with_recs[which(yrs_with_recs$Yr_rec %in% zero_catches$Year),]



# Abundance estimates -----------------------------------------------------

#IIIE = 35E-70E
#VIW = 170W-145W

#JARPA abundance
#In option 1: Matsuoka has B (3,4) and both were in odd years, and C(5,6) and both were in even years
#In option 2: Matsuoka has A (3) in odd years, B(4,5) in both even and odd and C(6) in even years, assuming even years for B
#For Matsuoka issues with 0s in area C

#For Hamabe: 
#In option 1: B is odd and C is even
#In option 2: only area B so using even years

#using just Hamabe for right now
Abund_newgrps_1<-abundance_dat_by_IWC(opt1, Mat = F, opt = 1)

#plotting
ggplot() + geom_point(data = Abund_newgrps_1, aes(x = as.factor(Year), y = Total, color = Area, shape = as.factor(Source))) + 
  facet_grid(Source ~ Area, labeller = labeller(Area = basin.labs1))


Abund_newgrps_2<-abundance_dat_by_IWC(opt2, Mat = F, opt = 2)

#plotting
ggplot() + geom_point(data = Abund_newgrps_2, aes(x = as.factor(Year), y = Total, color = Area, shape = as.factor(Source))) + 
  facet_grid(Source ~ Area, labeller = labeller(Area = basin.labs2))

write_csv(Abund_newgrps_1, "Custom_Model/Data/Abund_NewGrps_1.csv")
write_csv(Abund_newgrps_2, "Custom_Model/Data/Abund_NewGrps_2.csv")
