#Functions for data summary and cleaning

releases_tab_newgrps<-function(groups){
  ##groups are indicated as numbers and are IWC areas
  releases<-tibble("Year"= marking_years)
  
  #number of releases in each year in each basin
  temp<-ABW_marks_norecaps %>%
    mutate(Area = recode(IWCarea_mark, !!!groups)) %>%
    group_by(Seasn_mark, Area) %>%
    summarise(count = n()) %>%
    pivot_wider(names_from = Area, values_from = count) %>%
    replace_na(list(A = 0, B = 0, C = 0)) %>%
    filter(Seasn_mark <= 1972) %>%
    rename(Year = Seasn_mark)
  
  temp1<-ABW_marks_recaps %>%
    mutate(Area = recode(IWCarea_mark, !!!groups)) %>%
    group_by(Seasn_mark, Area) %>%
    summarise(count = n()) %>%
    pivot_wider(names_from = Area, values_from = count) %>%
    replace_na(list(A = 0, B = 0, C = 0)) %>%
    filter(Seasn_mark <= 1972) %>%
    rename(Year = Seasn_mark) %>%
    mutate(Year = Year + 1900)
  
  temp_all<-bind_rows(temp, temp1) %>%
    group_by(Year) %>%
    summarise_all(sum, na.rm = TRUE)
  
  #summary(temp_all)
  
  #adding to releases tibble
  releases<-left_join(releases,temp_all)
  releases<-replace_na(releases, list(A = 0, B = 0, C = 0))
  
  #adding years with 0 releases
  all_years<-seq(1926, 1972, by = 1)
  no_releases<-all_years[which(all_years %!in% releases$Year)]
  releases<-add_row(releases, Year = no_releases, A = 0,B = 0, C = 0)
  
  return(releases)
}

recoveries_tab_newgrps<-function(groups){
  cols<-c("AA", "AB", "AC", "BA", "BB", "BC", "CA", "CB", "CC")
  recoveries<-tibble(Yr_rec = all_years)
  recoveries<-filter(recoveries, recoveries$Yr_rec != 1926) #recoveries start in 1927 because first marks were 1926 and can't be recovered in the same season
  
  temp2<-ABW_marks_recaps %>%
    mutate(Area = recode(IWCarea_mark, !!!groups), AreaRec = recode(IWCarea_rec, !!!groups)) %>%
    group_by(Seasn_rec, Area, AreaRec) %>%
    summarise(count = n()) %>%
    pivot_wider(names_from = c(Area,AreaRec), values_from = count, names_sep = "")
    
  idx<-which(!(cols %in% colnames(temp2)[-1]))
  
  if(length(idx) > 0){
    newcols<-cols[idx]
    #print(newcols)
    for(i in 1:length(newcols)){
      #print(newcols[i])
      temp2[newcols[i]]<-0
    }
  }
  #print(temp2)
  temp2<-temp2 %>% 
    replace_na(list(AA = 0, AB = 0, AC = 0,BB = 0,BA = 0,BC = 0, CC = 0, CA = 0, CB = 0))%>%
    dplyr::select(sort(colnames(.)))%>%
    relocate(Seasn_rec, .before = AA) %>%
    rename(Yr_rec = Seasn_rec)
  
  temp2$Yr_rec<-temp2$Yr_rec + 1900
  
  #summary(temp2)
  
  recoveries<-left_join(recoveries, temp2)
  recoveries<-replace_na(recoveries, list(AA = 0, AB = 0, AC = 0,BB = 0,BA = 0,BC = 0, CC = 0, CA = 0, CB = 0))
  return(recoveries)
}

cor_catches_by_group<-function(brks, opt){
  #number of catches per year in each ocean group
  #originally:
  #A = Atlantic, B = Indian, C = Pacific 
  #Atlantic is -67.26 to 20 Longitudes
  #Indian is 20 to 146.9167 Longitude
  #Pacific is everything else (146.9 to -67) 
  #now need to input specific breaks
  
  #cutpoints<-c(-67.267, 20, 146.9167)
  if(opt == 1){ #can keep same order
    cutpoints<-brks
    catches_basins<-ABW_Catches %>%
      rename(OldYear = Year, Year = Seasn) %>%
      group_by(Year) %>%
      count(cut(Longitude, breaks = cutpoints, labels = c('A', 'B')))%>%
      rename(Basin = 'cut(Longitude, breaks = cutpoints, labels = c("A", "B"))') %>%
      mutate(Basin = fct_explicit_na(Basin, na_level = "C"))
    
  } else if(opt == 2){ #need to reorder boundaries so they work
    cutpoints<-c(brks[3], brks[1:2])
    catches_basins<-ABW_Catches %>%
      rename(OldYear = Year, Year = Seasn) %>%
      group_by(Year) %>%
      count(cut(Longitude, breaks = cutpoints, labels = c('C', 'A')))%>%
      rename(Basin = 'cut(Longitude, breaks = cutpoints, labels = c("C", "A"))') %>%
      mutate(Basin = fct_explicit_na(Basin, na_level = "B"))
    
  }
  
  total_year_catch<-catches_basins %>%
    group_by(Year) %>%
    summarise(total = sum(n))
  
  prop_catches_basins<- catches_basins %>% 
    left_join(total_year_catch, by = "Year") %>%
    summarise(prop = n/total, .groups = "rowwise")
  
  
  catches_basins$prop<-prop_catches_basins$prop
  
  final_catches<- catches_basins %>%
    left_join(TB_2004, by = "Year") %>%
    mutate(final_catch = Catches*prop)
  
  
  Catches<-final_catches %>% 
    pivot_wider(names_from = Basin, values_from = final_catch)
  
  Catches<-Catches %>%
    replace_na(list(A = 0, B = 0, C = 0)) %>%
    group_by(Year) %>%
    summarise_at(c("A", "B", "C"), sum)
  
  #adding in catches before 1913, assuming all are in the "A" group for both (Atlantic is mostly area 2)
  Catches<-Catches %>% add_row(Year = TB_2004$Year[TB_2004$Year<1913], A = TB_2004$Catches[TB_2004$Year<1913], 
                               B = 0, C = 0) %>% arrange(Year)
  return(Catches)
}

sigma_IWC<-function(Ns, CVs, Ntot){
  #same as for basins but props are both 1 because whole area is included
  sigs<-(CVs*Ns)^2 #converting CVs into sigmas
  sigtot<-sqrt(sum(sigs))
  endCV<-sigtot/Ntot
  return(endCV)
}

abundance_dat_by_IWC<-function(groups, Mat = T, opt = 1){
 Pop_TB<-TB_2007 %>% mutate(Area = recode(`IWC Area`, !!!groups)) %>%
    filter(Survey != "CPIII_star") %>%
    group_by(Area, Survey) %>%
    mutate(Year = floor(mean(Year))) %>% #using first year 
    ungroup() %>%
    group_by(Area, Year) %>%
    summarise(Total = sum(PopSize), CV = sigma_IWC(PopSize, CV, Total)) %>%
    add_column(Source = rep("Branch", nrow(.)))
 #JARPA abundance
    #In option 1: Matsuoka has B (3,4) and both were in odd years, and C(5,6) and both were in even years
    #In option 2: Matsuoka has A (3) in odd years, B(4,5) in both even and odd and C(6) in even years, assuming even years for B
    if(Mat){ #True if using Matsuoka
      if(opt == 1){ #option 1
        Pop_Mat<-Mat_Abund_filtered %>% mutate(Area = recode(AreaNum, !!!groups)) %>%
          #group_by(Area) %>%
          #mutate(Year = ifelse(Season %% 2 != 1, Season, Season +1)) %>% #using even years
          #ungroup() %>%
          mutate(Year = Season) %>%
          group_by(Area, Year) %>%
          summarise(Total = sum(N), CV = sigma_IWC(N, CV, Total)) %>%
          add_column(Source = rep("Matsouka", nrow(.)))
      } else{ #option 2
        Pop_Mat<-Mat_Abund_filtered %>% mutate(Area = recode(AreaNum, !!!opt2)) %>%
          group_by(Area) %>%
          mutate(Year = ifelse(Area == "B", ifelse(Season %% 2 != 1, Season, Season +1), Season)) %>% #using even years for group B and actual years everywhere else
          ungroup() %>%
          group_by(Area, Year) %>%
          summarise(Total = sum(N), CV = sigma_IWC(N, CV, Total)) %>%
          add_column(Source = rep("Matsouka", nrow(.)))
      }
     PopTotal<-bind_rows(Pop_TB, Pop_Mat)
    } else{#using Hamabe
      if(opt == 1){ #B is odd and C is even
      Pop_Ham<-HambeAbund %>% mutate(Area = recode(AreaNum, !!!groups)) %>%
        mutate(Year = Season) %>%
        group_by(Area, Year) %>%
        summarise(Total = sum(N), CV = sigma_IWC(N, CV, Total)) %>%
        add_column(Source = rep("Hamabe", nrow(.)))
      } else{ #option 2, B is only one represented and it's even and odd so using even
        Pop_Ham<-HambeAbund %>% mutate(Area = recode(AreaNum, !!!opt2)) %>%
          group_by(Area) %>%
          mutate(Year = ifelse(Area == "B", ifelse(Season %% 2 != 1, Season, Season +1), Season)) %>% #using even years for group B and actual years everywhere else
          ungroup() %>%
          group_by(Area, Year) %>%
          summarise(Total = sum(N), CV = sigma_IWC(N, CV, Total)) %>%
          add_column(Source = rep("Hamabe", nrow(.)))
    }
   PopTotal<-bind_rows(Pop_TB, Pop_Ham) 
  }
  return(PopTotal)
}



