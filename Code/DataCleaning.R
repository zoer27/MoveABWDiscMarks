#Preparing data for the model 
#Note that this does not run because the raw discovery mark data is not included in the repository
#Discovery mark data available from statistics@iwc.int
library(tidyverse)
library(rnaturalearth) #for plotting
library(sf) #for plotting
# Functions ---------------------------------------------------------------

source("FINAL Code/Functions/Data_Manip_Functions.R")
source("FINAL Code/Functions/UsefulFunctions.R") #blue whale functions from Trevor Branch

# Mark Recovery Data ------------------------------------------------------
ABW_marks_recaps<-read_csv("FINAL Code/Data/ABW_Marks_Recaps.csv") #just contains ABW marks with recaptures and recapture information 
ABW_marks_norecaps<-read_csv("FINAL Code/Data/ABW_Marks_NoRecaps.csv") #contains marks that were placed but never recaptured

#removing  a bunch of extra columns
ABW_marks_norecaps<-ABW_marks_norecaps[,1:12]
head(ABW_marks_norecaps)

#recaps imported a bunch of extra rows
ABW_marks_recaps<-ABW_marks_recaps[1:95,]


#N/S and E/W are NAs for some reason in the norecaps--getting this info from the Master sheet

#combining two master/main datasets
main<-read_csv("FINAL Code/Data/MarkMaster_5a.csv", skip = 1)
main_2<-read_csv("FINAL Code/Data/MarkMaster_2a.csv", skip = 1)
main<-select(main, -c(Note, Comment))
main_2<-select(main_2, - c(Expt, Comment))
main<-rbind(main, main_2)

#incorrectly inputted ID for one observation (25810 should be 25811)
temp<-which(main$No. == "25810")
main[temp[2],1]<-"25811"

#adding NSEW to no_recaps dataframe
ABW_marks_norecaps<-get_direction(ABW_marks_norecaps, main) #from Data Manipulation Functions
summary(ABW_marks_norecaps)

#double checking that all marks with recaptures are not in the dataframe without recaptures
No_Recap_names<-unique(ABW_marks_norecaps$MarkID)
Recap_names<-unique(ABW_marks_recaps$MarkNo.)

which(Recap_names %in% No_Recap_names)
Recap_names[46] #one duplicate
ABW_marks_recaps[ABW_marks_recaps$MarkNo. == Recap_names[46],] 
#has recapture information so removing from no recaptures
idx<-which(No_Recap_names %in% Recap_names)
ABW_marks_norecaps<-ABW_marks_norecaps[-idx,]

#check that it worked
No_Recap_names<-unique(ABW_marks_norecaps$MarkID)
Recap_names<-unique(ABW_marks_recaps$MarkNo.)
which(Recap_names %in% No_Recap_names)

#adding season
ABW_marks_norecaps$Seasn_mark<-ifelse(ABW_marks_norecaps$Month <= 7, ABW_marks_norecaps$Year-1, ABW_marks_norecaps$Year)
ABW_marks_recaps$Seasn_mark<-ifelse(ABW_marks_recaps$Mon_mark <= 7, ABW_marks_recaps$Yr_mark -1, ABW_marks_recaps$Yr_mark)

#using Trevor's get subspecies to assign marks to subspecies of blue whales
subsp<-get.subspecies.ZR(ABW_marks_norecaps$Latitude, ABW_marks_norecaps$Longitude)
ABW_marks_norecaps$Sp<-subsp
table(ABW_marks_norecaps$Sp) 
#removing non-ABW
ABW_marks_norecaps<-ABW_marks_norecaps %>%
  filter(Sp == "Antarctic")
#removing the points between 20S and 0S (not ABW)
ABW_marks_norecaps<-filter(ABW_marks_norecaps, !(Latitude < 0 & Latitude > -20))
table(ABW_marks_norecaps$Sp)

#looking at marks
subsp_mark<-get.subspecies.ZR(ABW_marks_recaps$lat_mark, ABW_marks_recaps$long_mark)
table(subsp_mark) #all Antarctic
subsp_rec<-get.subspecies.ZR(ABW_marks_recaps$Lat_rec, ABW_marks_recaps$Long_rec)
table(subsp_rec) #all Antarctic


#plotting to check mark basin designations
locs_sf_ABW<-ABW_marks_norecaps %>%
  st_as_sf(coords = c('Longitude', 'Latitude')) %>%
  st_set_crs(4326)
locs_sf_ABW_mark<-ABW_marks_recaps %>%
  st_as_sf(coords = c('long_mark', 'lat_mark')) %>%
  st_set_crs(4326)
locs_sf_ABW_recap<-ABW_marks_recaps %>% 
  st_as_sf(coords = c("Long_rec","Lat_rec")) %>%
  st_set_crs(4326)

world<-ne_countries(scale = "medium", returnclass = "sf")

whale_marks<-ggplot() + 
  theme_bw() + 
  geom_sf(aes(color = MARK_basin),  data = locs_sf_ABW) +
  geom_sf(aes(), data = world) +
  coord_sf(xlim=c(-180,180),ylim=c(-80,20), expand = FALSE) + 
  theme(legend.position = "bottom")
whale_marks

marks_placed<-ggplot() + 
  theme_bw() + 
  geom_sf(aes(color = MARK_basin),  data = locs_sf_ABW_mark) +
  geom_sf(aes(), data = world) +
  coord_sf(xlim=c(-180,180),ylim=c(-80,20), expand = FALSE) + 
  theme(legend.position = "bottom")

marks_rec<-ggplot() + 
  theme_bw() + 
  geom_sf(aes(color = recMARK_basin),  data = locs_sf_ABW_recap) +
  geom_sf(aes(), data = world) +
  coord_sf(xlim=c(-180,180),ylim=c(-80,20), expand = FALSE) + 
  theme(legend.position = "bottom")

marks_rec #some mismatches here

#Basins
#Atlantic is -67.26 to 20 Longitudes
#Indian is 20 to 146.9167 Longitude
#Pacific is everything else (146.9 to -67) 
cutpoints<-c(-67.267, 20, 146.9167)

#2 marks mistakenly indicated as recovered in the Indian when they are in the Pacific
ABW_marks_recaps %>% filter(Long_rec > -140 & Long_rec < -120 & recMARK_basin == "B")

idx<-ABW_marks_recaps %>% filter(Long_rec > -140 & Long_rec < -120 & recMARK_basin == "B") %>%
  dplyr::select(MarkNo.) %>% unlist()

#remarking them as pacific
ABW_marks_recaps[ABW_marks_recaps$MarkNo. %in% idx, "recMARK_basin"] <-c("C", "C")

#replotting
locs_sf_ABW_recap<-ABW_marks_recaps %>% 
  st_as_sf(coords = c("Long_rec","Lat_rec")) %>%
  st_set_crs(4326)

marks_rec<-ggplot() + 
  theme_bw() + 
  geom_sf(aes(color = recMARK_basin),  data = locs_sf_ABW_recap) +
  geom_sf(aes(), data = world) +
  coord_sf(xlim=c(-180,180),ylim=c(-80,20), expand = FALSE) + 
  theme(legend.position = "bottom")

marks_rec

#save
#write_csv(ABW_marks_norecaps, file = "FINAL Code/Data/ABW_marks_norecaps_clean.csv")
#write_csv(ABW_marks_recaps, file = "FINAL Code/Data/ABW_marks_recaps_clean.csv")

#seperating out same season and different season recaptures
same_seasons<-ABW_marks_recaps %>%filter(Seasn_mark == Seasn_rec)
diff_seasons<-ABW_marks_recaps %>% filter(Seasn_mark != Seasn_rec)

#save
write_csv(diff_seasons, "FINAL Code/Data/ABW_marks_recaps_diffseasons.csv")
write_csv(same_seasons, "FINAL Code/Data/ABW_marks_recaps_sameseasons.csv")


# Marking and Recovery Tables for Models ----------------------------------

#summarizing mark data
#last season is 1972
main$seasn<-ifelse(main$Mon <= 7, main$Yr -1, main$Yr)
marking_years<-unique(main$seasn)
marking_years<-marking_years+1900
marking_years<-marking_years[which(marking_years<= 1972 & marking_years>= 1926)]

#tibble for releases
releases<-tibble("Year"= marking_years)

#number of releases in each year in each basin
temp<-ABW_marks_norecaps %>%
  group_by(Seasn_mark, MARK_basin) %>%
  summarise(count = n()) %>%
  pivot_wider(names_from = MARK_basin, values_from = count) %>%
  replace_na(list(A = 0, B = 0, C = 0)) %>%
  filter(Seasn_mark <= 1972) %>%
  rename(Year = Seasn_mark)

summary(temp)

temp1<-diff_seasons %>%
  group_by(Seasn_mark, MARK_basin) %>%
  summarise(count = n()) %>%
  pivot_wider(names_from = MARK_basin, values_from = count) %>%
  replace_na(list(A = 0, B = 0, C = 0)) %>%
  filter(Seasn_mark <= 1972) %>%
  rename(Year = Seasn_mark) %>%
  mutate(Year = Year + 1900)

summary(temp1)

temp_all<-bind_rows(temp, temp1) %>%
  group_by(Year) %>%
  summarise_all(sum, na.rm = TRUE)

summary(temp_all)

#adding to releases tibble
releases<-left_join(releases,temp_all)
releases<-replace_na(releases, list(A = 0, B = 0, C = 0))
head(releases)
summary(releases)

#addingn years with 0 releases
all_years<-seq(1926, 1972, by = 1)
no_releases<-all_years[which(all_years %!in% releases$Year)]
releases<-add_row(releases, Year = no_releases, A = 0,B = 0, C = 0)

#recoveries, grouped by where they were released
recoveries<-tibble(Yr_rec = unique(releases$Year))
recoveries<-filter(recoveries, recoveries$Yr_rec != 1926) #recoveries start in 1927 because first marks were 1926 and can't be recovered in the same season

temp2<-diff_seasons %>%
  group_by(Seasn_rec, MARK_basin, recMARK_basin) %>%
  summarise(count = n()) %>%
  pivot_wider(names_from = c(MARK_basin,recMARK_basin), values_from = count, names_sep = "") %>%
  add_column(AC = 0, CA = 0) %>%
  replace_na(list(AA = 0, AB = 0, AC = 0,BB = 0,BA = 0,BC = 0, CC = 0, CA = 0, CB = 0))%>%
  dplyr::select(sort(colnames(.)))%>%
  relocate(Seasn_rec, .before = AA) %>%
  rename(Yr_rec = Seasn_rec)

temp2$Yr_rec<-temp2$Yr_rec + 1900

summary(temp2)

recoveries<-left_join(recoveries, temp2)
recoveries<-replace_na(recoveries, list(AA = 0, AB = 0, AC = 0,BB = 0,BA = 0,BC = 0, CC = 0, CA = 0, CB = 0))
head(recoveries)
summary(recoveries)

#adding in years where there were no recoveries
all_years<-seq(1927, 1972, by = 1)
no_recoveries<-all_years[which(all_years %!in% recoveries$Yr_rec)]
recoveries<-add_row(recoveries, Yr_rec = no_recoveries, AA = 0, AB = 0, AC = 0,BB = 0,BA = 0,BC = 0, CC = 0, CA = 0, CB = 0)
tail(recoveries)
head(recoveries)

#saving
write_csv(recoveries, file = "FINAL Code/Data/mark_recoveries_3groups.csv")
write_csv(releases, file = "FINAL Code/Data/mark_releases.csv")


# Catch Data --------------------------------------------------------------
files <- paste0("Custom_Model/Data/IWC_Individual_v7.1/", c( "IO.csv", "NA.csv", "NP.csv", "SA.csv", 
                                                             "SHL.csv", "SHP1.csv", "SHP2.csv",
                                                             "SP.csv", "SU.csv"))
BWOrig<-get.blue.data.z(files = files)
names(BWOrig) #X18 is N/S, X22 is E/W 
#renaming some columns to make things easier and to work with Trevor's functions
BWOrig<- BWOrig %>% rename(X18 = ...18, X22 = ...22, Mn = Mn...17, Mn_1 = Mn...21, L.u = `L-u`)
table(BWOrig$Sp) #4 is True Blue and 15 is Pygmy Blue

#getting rid of catches without location data: 
BWOrig_loc<-BWOrig[-which(is.na(BWOrig$X18)),]
which(is.na(BWOrig_loc$X22)) #checking to make sure the above worked
which(is.na(BWOrig_loc$X18))
filter(BWOrig_loc, Lat == 0 & Mn == 0 & Lon == 0 & Mn_1 == 0) #another check

#View(filter(BWOrig_loc, Lon == 0)) #still kept in locations which were actually caught at 0 longitude

#calculating whaling season
BWOrig_loc$Seasn<-ifelse(BWOrig_loc$Mon <= 7, BWOrig_loc$Year -1, BWOrig_loc$Year)
summary(BWOrig_loc)
hist(BWOrig_loc$Mon)

#getting decimal latlon
latlon<-decimal.latlon(lat = BWOrig_loc$Lat, latmin = BWOrig_loc$Mn, latNS = BWOrig_loc$X18, 
                       lon = BWOrig_loc$Lon, lonmin = BWOrig_loc$Mn_1,lonEW = BWOrig_loc$X22)
BWOrig_loc$Latitude<-latlon$Lat
BWOrig_loc$Longitude<-latlon$Lon
#getting subspecies
spp<-get.subspecies(bluedata = BWOrig_loc, declatlon = TRUE)

BWID<-seq(1:nrow(BWOrig_loc))

BWOrig_loc$Subsp<-spp
#saving ABW catches with all information--including Durban
ABW_catch_clean<-BWOrig_loc %>% filter(Subsp == "Antarctic" | Subsp == "Unknown Durban")
write_csv(ABW_catch_clean, file = "FINAL Code/Data/ABW_catch_clean.csv")

#getting just the info we need for summaries
BWCatch<-cbind(BWID, BWOrig_loc$Seasn, latlon, Subsp = as.factor(spp))
head(BWCatch)
summary(BWCatch)

#just Antarctic blue whales--including Durban
ANT_BWCatch<-filter(BWCatch,Subsp == "Antarctic" | Subsp == "Unknown Durban") %>% rename(Seasn = "BWOrig_loc$Seasn")
summary(ANT_BWCatch)


#number of catches per year in each ocean basin
#A = Atlantic, B = Indian, C = Pacific 
#Atlantic is -67.26 to 20 Longitudes
#Indian is 20 to 146.9167 Longitude
#Pacific is everything else (146.9 to -67) 
cutpoints<-c(-67.267, 20, 146.9167)

catches_basins<-ANT_BWCatch %>%
  rename(Year = Seasn) %>%
  group_by(Year) %>%
  count(cut(Lon, breaks = cutpoints, labels = c('Atlantic', 'Indian')))%>%
  rename(Basin = 'cut(Lon, breaks = cutpoints, labels = c("Atlantic", "Indian"))') %>%
  mutate(Basin = fct_explicit_na(Basin, na_level = "Pacific"))

levels(catches_basins$Basin)
head(catches_basins)
#write_csv(catches_basins, "FINAL Code/Data/catch_basins_not_corrected.csv")
#correcting for catches without individual catch information 
#Catch series from Branch et al 2004
TB_2004<-read_csv("FINAL Code/Data/TB2004_catches.csv")
TB_2004$Year<-TB_2004$Year-1 #Trevor uses end of whaling season as his years but I use start
#only need catches before 1973
TB_2004<-filter(TB_2004, Year < 1973)

#getting total number of catches by year
total_year_catch<-catches_basins %>%
  group_by(Year) %>%
  summarise(total = sum(n))
total_year_catch

#getting proportion of catches in each basin from individual database
prop_catches_basins<- catches_basins %>% 
  left_join(total_year_catch, by = "Year") %>%
  summarise(prop = n/total, .groups = "rowwise")
prop_catches_basins

#adding props to catches dataframe
catches_basins$prop<-prop_catches_basins$prop

#dividing Full catches (Branch 2004) by props in each basin 
final_catches<- catches_basins %>%
  left_join(TB_2004, by = "Year") %>%
  mutate(final_catch = Catches*prop)


Catches<-final_catches %>% 
  pivot_wider(names_from = Basin, values_from = final_catch)

Catches<-Catches %>%
  replace_na(list(Atlantic = 0, Pacific = 0, Indian = 0)) %>%
  group_by(Year) %>%
  summarise_at(c("Atlantic", "Indian", "Pacific"), sum)
head(Catches)
#adding Catches before 1913 in assuming all are Atlantic
Catches<-Catches %>% add_row(Year = TB_2004$Year[TB_2004$Year<1913], Atlantic = TB_2004$Catches[TB_2004$Year<1913], 
                             Indian = 0, Pacific = 0) %>% arrange(Year)
summary(Catches)

#save
write_csv(Catches, file = "FINAL Code/Data/catchcorrected.csv")


# Abundance Estimates -----------------------------------------------------
#Branch 2007 circumpolar estimates
TB_2007<-read_csv("FINAL Code/Data/TB_2007_Abundance.csv")[1:19, 1:5]

#Converting IWC estimates to estimates by basin
#function for converting CVs
sigma_tot<-function(Ns, CVs, props, Ntot){
  props_n<-props/sum(props)
  sigs<-(CVs*Ns)^2
  props2<-props_n^2
  sigtot<-sqrt(sum(sigs*props2))
  endCV<-sigtot/Ntot
  return(endCV)
}
#sigma_tot(c(25,16,219), c(0.8, 0.81, 0.61), Aprop, 91)
Abundance<-tibble("Year" = NA, "Basin" = c(rep("Atlantic", 3), rep("Indian", 3), rep("Pacific", 3)), "PopSize" = NA, "CV" = NA)
#Atlantic is all of 2, 7/60 of 1, and 20/70 of 3
#Atlantic props
Aprop<-c((7/60), 1, (2/7))

totalA1<-TB_2007 %>% filter(`IWC Area` %in% c(1,2,3) & Survey == "CPI")%>%
  mutate(FinYear = (Aprop[1]*Year[1] + Aprop[2]*Year[2] + Aprop[3]*Year[3])/(sum(Aprop))) %>% 
  mutate(FinYear = round(FinYear,0)) %>%
  summarise(Year = unique(FinYear), total = Aprop[1]*PopSize[1] + Aprop[2]*PopSize[2] + Aprop[3]*PopSize[3],
            CV = sigma_tot(PopSize, CV, Aprop, total)) %>%
  mutate(total = round(total, 0))

Abundance[1,c(1,3:4)]<-totalA1

totalA2<-TB_2007 %>% filter(`IWC Area` %in% c(1,2,3) & Survey == "CPII")%>%
  mutate(FinYear = (Aprop[1]*Year[1] + Aprop[2]*Year[2] + Aprop[3]*Year[3])/(sum(Aprop))) %>% 
  mutate(FinYear = round(FinYear,0)) %>%
  summarise(Year = unique(FinYear), total = Aprop[1]*PopSize[1] + Aprop[2]*PopSize[2] + Aprop[3]*PopSize[3],
            CV = sigma_tot(PopSize, CV, Aprop, total)) %>%
  mutate(total = round(total, 0))

Abundance[2, c(1,3:4)]<-totalA2

totalA3<-TB_2007 %>% filter(`IWC Area` %in% c(1,2,3) & Survey == "CPIII")%>%
  mutate(FinYear = (Aprop[1]*Year[1] + Aprop[2]*Year[2] + Aprop[3]*Year[3])/(sum(Aprop))) %>% 
  mutate(FinYear = round(FinYear,0)) %>%
  summarise(Year = unique(FinYear), total = Aprop[1]*PopSize[1] + Aprop[2]*PopSize[2] + Aprop[3]*PopSize[3],
            CV = sigma_tot(PopSize, CV, Aprop, total)) %>%
  mutate(total = round(total, 0))


Abundance[3, c(1,3:4)]<-totalA3

#Indian is 5/7 of Area 3, all of 4, and 17/60 of 5
Iprop<-c((5/7), 1, (17/60))

totalI1<-TB_2007 %>% filter(`IWC Area` %in% c(3,4,5) & Survey == "CPI")%>%
  mutate(FinYear = (Iprop[1]*Year[1] + Iprop[2]*Year[2] + Iprop[3]*Year[3])/(sum(Iprop))) %>% 
  mutate(FinYear = round(FinYear,0)) %>%
  summarise(Year = unique(FinYear), total = Iprop[1]*PopSize[1] + Iprop[2]*PopSize[2] + Iprop[3]*PopSize[3],
            CV = sigma_tot(PopSize, CV, Iprop, total)) %>%
  mutate(total = round(total, 0))

Abundance[4, c(1,3:4)]<-totalI1

totalI2<-TB_2007 %>% filter(`IWC Area` %in% c(3,4,5) & Survey == "CPII")%>%
  mutate(FinYear = (Iprop[1]*Year[1] + Iprop[2]*Year[2] + Iprop[3]*Year[3])/(sum(Iprop))) %>% 
  mutate(FinYear = round(FinYear,0)) %>%
  summarise(Year = unique(FinYear), total = Iprop[1]*PopSize[1] + Iprop[2]*PopSize[2] + Iprop[3]*PopSize[3],
            CV = sigma_tot(PopSize, CV, Iprop, total)) %>%
  mutate(total = round(total, 0))

Abundance[5, c(1,3:4)]<-totalI2

totalI3<-TB_2007 %>% filter(`IWC Area` %in% c(3,4,5) & Survey == "CPIII")%>%
  mutate(FinYear = (Iprop[1]*Year[1] + Iprop[2]*Year[2] + Iprop[3]*Year[3])/(sum(Iprop))) %>% 
  mutate(FinYear = round(FinYear,0)) %>%
  summarise(Year = unique(FinYear), total = Iprop[1]*PopSize[1] + Iprop[2]*PopSize[2] + Iprop[3]*PopSize[3],
            CV = sigma_tot(PopSize, CV, Iprop, total)) %>%
  mutate(total = round(total, 0))

Abundance[6, c(1,3:4)]<-totalI3

#Pacific is 53/60 of Area 1, 43/60 of Area 5 and all of Area 6
Pprop<-c((53/60), (43/60), 1)

totalP1<-TB_2007 %>% filter(`IWC Area` %in% c(1,5,6) & Survey == "CPI")%>%
  mutate(FinYear = (Pprop[1]*Year[1] + Pprop[2]*Year[2] + Pprop[3]*Year[3])/(sum(Pprop))) %>% 
  mutate(FinYear = round(FinYear,0)) %>%
  summarise(Year = unique(FinYear), total = Pprop[1]*PopSize[1] + Pprop[2]*PopSize[2] + Pprop[3]*PopSize[3],
            CV = sigma_tot(PopSize, CV, Pprop, total)) %>%
  mutate(total = round(total, 0))

Abundance[7, c(1,3:4)]<-totalP1

totalP2<-TB_2007 %>% filter(`IWC Area` %in% c(1,5,6) & Survey == "CPII")%>%
  mutate(FinYear = (Pprop[1]*Year[1] + Pprop[2]*Year[2] + Pprop[3]*Year[3])/(sum(Pprop))) %>% 
  mutate(FinYear = round(FinYear,0)) %>%
  summarise(Year = unique(FinYear), total = Pprop[1]*PopSize[1] + Pprop[2]*PopSize[2] + Pprop[3]*PopSize[3],
            CV = sigma_tot(PopSize, CV, Pprop, total)) %>%
  mutate(total = round(total, 0))

Abundance[8, c(1,3:4)]<-totalP2

totalP3<-TB_2007 %>% filter(`IWC Area` %in% c(1,5,6) & Survey == "CPIII")%>%
  mutate(FinYear = (Pprop[1]*Year[1] + Pprop[2]*Year[2] + Pprop[3]*Year[3])/(sum(Pprop))) %>% 
  mutate(FinYear = round(FinYear,0)) %>%
  summarise(Year = unique(FinYear), total = Pprop[1]*PopSize[1] + Pprop[2]*PopSize[2] + Pprop[3]*PopSize[3],
            CV = sigma_tot(PopSize, CV, Pprop, total)) %>%
  mutate(total = round(total, 0))

Abundance[9, c(1,3:4)]<-totalP3

#Matsuoka and Hakamada 2014
#IIIE = 35E-70E
#VIW = 170W-145W

Mat_Abund_2014<-read_csv("FINAL Code/Data/Matsouka_Abundance.csv")

#excluding years where only areas IV and V were surveyed
Mat_Abund_filtered<-Mat_Abund_2014 %>% filter(Season >= 1995)

Ind_Mat_Prop<-c(1,1,(17/60)) #all of IIIE, IV, and 17/60 of V, using odd years since IIIE and IV are odd years 
Pac_Mat_Prop<-c((43/60), 1) #43/60 of V, and all of VIW, using even years since V and VIW are even years

totalI<-Mat_Abund_filtered %>% filter(Area %in% c("IIIE", "IV", "V"))%>%
  mutate(Year = ifelse(Season %% 2 != 1, Season -1, Season)) %>% #using odd years
  group_by(Year) %>%
  summarise(Year = unique(Year), total = Ind_Mat_Prop[1]*N[1] + Ind_Mat_Prop[2]*N[2] + Ind_Mat_Prop[3]*N[3], 
            CV = sigma_tot(N, CV, Ind_Mat_Prop, total)) %>%
  mutate(total = round(total, 0)) %>% 
  rename(PopSize = total) %>%
  add_column(Basin = rep("Indian", nrow(.))) %>%
  add_column(Source = rep("Matsouka", nrow(.)))

totalP<-Mat_Abund_filtered %>% filter(Area %in% c("V", "VIW"))%>%
  mutate(Year = ifelse(Season %% 2 != 1, Season, Season +1)) %>% #using even years
  group_by(Year) %>%
  summarise(Year = unique(Year), total = Pac_Mat_Prop[1]*N[1] + Pac_Mat_Prop[2]*N[2], 
            CV = sigma_tot(N, CV, Pac_Mat_Prop, total)) %>%
  mutate(total = round(total, 0)) %>% 
  rename(PopSize = total) %>%
  add_column(Basin = rep("Pacific", nrow(.))) %>%
  add_column(Source = rep("Matsouka", nrow(.)))

Abundance<-Abundance %>% add_column(Source = rep("Branch", nrow(.)))
AbundanceNew_Mat<-bind_rows(Abundance, totalI, totalP)

#save
write_csv(AbundanceNew_Mat, file = "FINAL Code/Data/basin_abundances_Mat.csv")


#Hamabe 2023 Abundance estimates
#Just covers areas IV and V

HambeAbund<-read_csv("FINAL Code/Data/Hambe2023Abundance.csv")

Ind_H_Prop<-c(1,(17/60)) #all IV, and 17/60 of V, using odd years IV are odd years 
Pac_H_Prop<-c((43/60)) #43/60 of V even years since V  are even years

totalI<-HambeAbund %>% filter(Area %in% c("IV", "V"))%>%
  mutate(Year = ifelse(Season %% 2 != 1, Season -1, Season)) %>% #using odd years
  group_by(Year) %>%
  summarise(Year = unique(Year), total = Ind_H_Prop[1]*N[1] + Ind_H_Prop[2]*N[2], 
            CV = sigma_tot(N, CV, Ind_H_Prop, total)) %>%
  mutate(total = round(total, 0)) %>% 
  rename(PopSize = total) %>%
  add_column(Basin = rep("Indian", nrow(.))) %>%
  add_column(Source = rep("Hamabe", nrow(.)))

totalP<-HambeAbund %>% filter(Area == "V")%>%
  mutate(Year = ifelse(Season %% 2 != 1, Season, Season +1)) %>% #using even years
  group_by(Year) %>%
  summarise(Year = unique(Year), total = Pac_H_Prop[1]*N[1], 
            CV = sigma_tot(N, CV, Pac_H_Prop, total)) %>%
  mutate(total = round(total, 0)) %>% 
  rename(PopSize = total) %>%
  add_column(Basin = rep("Pacific", nrow(.))) %>%
  add_column(Source = rep("Hamabe", nrow(.)))


AbundanceNew_Hamabe<-bind_rows(Abundance, totalI, totalP)

write_csv(AbundanceNew_Hamabe, file = "FINAL Code/Data/basin_abundances_withHamabe.csv")


