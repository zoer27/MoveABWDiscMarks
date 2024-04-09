#Data Summary
#Zoe Rand
#Last Updated: 10/17/23
#sets up mark recovery and catch data for use in models as well as looks as summary stats about them
#note that raw mark recovery data is not available in this repository
#catches and Discovery mark data are available statistics@iwc.int

library(tidyverse)


# Functions ---------------------------------------------------------------
source("FINAL Code/Functions/Data_manip_functions.R")

# Mark Recovery Different Seasons -----------------------------------------------------------
ABW_releases<-read_csv("FINAL Code/Data/mark_releases.csv") #table of releases
ABW_recoveries<-read_csv("FINAL Code/Data/mark_recoveries_3groups.csv") #table of recoveries
ABW_mark_recaps<-read_csv("FINAL Code/Data/ABW_marks_recaps_diffseasons.csv")

#releases
colSums(ABW_releases[,2:4]) #number released in each basin
sum(colSums(ABW_releases[,2:4])) #total released
colSums(ABW_releases[,2:4])/sum(colSums(ABW_releases[,2:4])) #p(basin)

#recoveries
colSums(ABW_recoveries[,2:10]) #total in each group
sum(colSums(ABW_recoveries[,2:10])) #total recoveries
A<-sum(colSums(ABW_recoveries[,2:4])) #total marked in Atlantic
B<-sum(colSums(ABW_recoveries[,5:7])) #total marked in Indian
C<-sum(colSums(ABW_recoveries[,8:10])) #total marked in Pacific

colSums(ABW_recoveries[,2:4])/A #proportions recovered in each basin
sum(colSums(ABW_recoveries[,2:4]))/colSums(ABW_releases[,2:4])[1]
colSums(ABW_recoveries[,5:7])/B
sum(colSums(ABW_recoveries[,5:7]))/colSums(ABW_releases[,2:4])[2]
colSums(ABW_recoveries[,8:10])/C
sum(colSums(ABW_recoveries[,8:10]))/colSums(ABW_releases[,2:4])[3]

#summary information
summary(ABW_mark_recaps$`Duration(days)_rec`) #range of mark duration
#doesn't make much snese in days so converting to years
datesMark<-ymd(paste0(1900 + ABW_mark_recaps$Yr_mark, "-", ABW_mark_recaps$Mon_mark, "-", ABW_mark_recaps$Day_mark))
datesRec<- ymd(paste0(1900 + ABW_mark_recaps$Yr_rec, "-", ABW_mark_recaps$Mon_rec, "-", ABW_mark_recaps$Day_rec))

yrs<-interval(datesMark, datesRec)/years(1)
summary(yrs)

#range of distance: same basin

ABW_mark_recaps %>% filter(MARK_basin == recMARK_basin) %>% 
  summary(`Dist(km)_rec`)

#range of distance: different basin

ABW_mark_recaps %>% filter(MARK_basin != recMARK_basin) %>% 
  summary(`Dist(km)_rec`)

#more info on longest: 
ABW_mark_recaps %>% filter(`Dist(km)_rec` == 6250)

# Mark Recovery Same Seasons ----------------------------------------------
same_seasons<-read_csv("FINAL Code/Data/ABW_marks_recaps_sameseasons.csv")
#table of recoveries
#organizing same season data like above
marking_years<-seq(1926, 1972, by = 1)
ss_releases<-tibble("Year"= marking_years)


#number of releases in each year in each basin
temp<-same_seasons %>%
  group_by(Seasn_mark, MARK_basin) %>%
  summarise(count = n()) %>%
  pivot_wider(names_from = MARK_basin, values_from = count) %>%
  replace_na(list(A = 0, B = 0, C = 0)) %>%
  filter(Seasn_mark <= 1972) %>%
  rename(Year = Seasn_mark)

summary(temp)

temp1<-same_seasons %>%
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
ss_releases<-left_join(ss_releases,temp_all)
ss_releases<-replace_na(ss_releases, list(A = 0, B = 0, C = 0))
head(ss_releases)
summary(ss_releases)
#sum(ss_releases[,-1])

#addingn years with 0 releases
all_years<-seq(1926, 1972, by = 1)
ss_no_releases<-all_years[which(all_years %!in% ss_releases$Year)]
ss_releases<-add_row(ss_releases, Year = ss_no_releases, A = 0,B = 0, C = 0)

#recoveries, grouped by where they were released
ss_recoveries<-tibble(Yr_rec = unique(ss_releases$Year))
ss_recoveries<-filter(ss_recoveries, ss_recoveries$Yr_rec != 1926) #recoveries start in 1927 because first marks were 1926 and can't be recovered in the same season

temp2<-same_seasons %>%
  group_by(Seasn_rec, MARK_basin, recMARK_basin) %>%
  summarise(count = n()) %>%
  pivot_wider(names_from = c(MARK_basin,recMARK_basin), values_from = count, names_sep = "") %>%
  add_column(AB = 0, CA = 0, BC = 0) %>%
  replace_na(list(AA = 0, AB = 0, AC = 0,BB = 0,BA = 0,BC = 0, CC = 0, CA = 0, CB = 0))%>%
  dplyr::select(sort(colnames(.)))%>%
  relocate(Seasn_rec, .before = AA) %>%
  rename(Yr_rec = Seasn_rec)

temp2$Yr_rec<-temp2$Yr_rec + 1900

summary(temp2)

ss_recoveries<-left_join(ss_recoveries, temp2)
ss_recoveries<-replace_na(ss_recoveries, list(AA = 0, AB = 0, AC = 0,BB = 0,BA = 0,BC = 0, CC = 0, CA = 0, CB = 0))
head(ss_recoveries)
summary(ss_recoveries)



#adding in years where there were no recoveries
all_years<-seq(1927, 1972, by = 1)
ss_no_recoveries<-all_years[which(all_years %!in% ss_recoveries$Yr_rec)]
ss_recoveries<-add_row(ss_recoveries, Yr_rec = ss_no_recoveries, AA = 0, AB = 0, AC = 0,BB = 0,BA = 0,BC = 0, CC = 0, CA = 0, CB = 0)
tail(ss_recoveries)
head(ss_recoveries)
#sum(ss_recoveries[,-1])

#summary stats
#releases
colSums(ss_releases[,2:4])
sum(colSums(ss_releases[,2:4]))
colSums(ss_releases[,2:4])/sum(colSums(ss_releases[,2:4]))

#recoveries

colSums(ss_recoveries[,2:10])
sum(colSums(ss_recoveries[,2:10]))
A<-sum(colSums(ss_recoveries[,2:4]))
B<-sum(colSums(ss_recoveries[,5:7]))
C<-sum(colSums(ss_recoveries[,8:10]))
print(c(A, B, C))
colSums(ss_recoveries[,2:4])/A

colSums(ss_recoveries[,5:7])/B

colSums(ss_recoveries[,8:10])/C

#summary information about the marks
summary(same_seasons)

hist(same_seasons$`Duration(days)_rec`)
summary(same_seasons$`Duration(days)_rec`) #range of mark duration
same_seasons %>% filter(`Duration(days)_rec` == 114) #more info about the longest mark

#range of distance
summary(same_seasons$`Dist(km)_rec`)


# Comparing catches with and without location data ------------------------
TB_2004<-read_csv("FINAL Code/Data/TB2004_catches.csv")
catch_basins<-read_csv("FINAL Code/Data/catch_basins_not_corrected.csv")

catch_location_year<-catch_basins %>% group_by(Year) %>% summarise(Tot = sum(n))
head(catch_location_year)

TB_2004$Year<-TB_2004$Year-1 #Trevor uses end of whaling season as his years but I use start

catches_together<-TB_2004 %>% left_join(catch_location_year) %>% 
  mutate(Prop_location = ifelse(Tot/Catches > 1, 1, Tot/Catches))

ggplot(catches_together) + geom_bar(aes(x = Year, y = Prop_location), stat = "identity") + 
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) + 
  scale_x_continuous(limits = c(1912, 1973), breaks = seq(1912, 1973, by = 5)) + 
  labs(x = "Season", y = "Proportion of catches with locations") + theme_classic()

#proportion of catches without location data
1-sum(catches_together$Tot, na.rm = TRUE)/sum(catches_together$Catches, na.rm = TRUE)


