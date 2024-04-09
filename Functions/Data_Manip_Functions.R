### Data Manipulation Functions for Custom Mark Recapture Models
#Zoe Rand

#%!in% finds values that are not equal to any in a vector
'%!in%' <- function(x,y)!('%in%'(x,y))

#get_direction() takes the mark data that is missing the E/W, N/S columns,
#and the MarkMaster data that has it and adds in the columns to the mark data
get_direction<-function(missing_df, main_df){
  for (i in 1:nrow(missing_df)){
    t<-which(missing_df$MarkID[i] == main_df$"No.")
    missing_df$"N/S"[i]<-main_df$N_S[t]
    missing_df$"E/W"[i]<-main_df$E_W[t]
  }
  return(missing_df)
}



#isEmpty
#tests whether a vector is empty (for use in other functions)
isEmpty<-function(x){
  return(identical(x, numeric(0)))
}

#stay fills in the probability of staying a location as 1-probability of moving
stay<-function(move_mat){
  #fills in probability of staying with 1-prob of moving
  totals<-colSums(move_mat, na.rm = TRUE)
  for (i in 1:nrow(move_mat)){
    for (j in 1:ncol(move_mat)){
      if (is.na(move_mat[i,j])){
        move_mat[i,j]<-(1-totals[i])
      }
    }
  }
  return(move_mat)
}



# Adapted Trevor's useful functions ------------------------------------
#Get subspecies for mark data
#Based on Function by Trevor Branch in Useful functions
#changed to use actual longitude and latitude columns 
#Removed the pygmy blue whale designation south of 60 and the ones from the catch database that are > 80 ft
get.subspecies.ZR <- function(Lat, Lon) { #give it actual latitude and longitude columns, assumes decimal lat long already
  latlon <-data.frame(Lat = Lat, Lon = Lon)
  
  subspecies <- vector(length=NROW(latlon))
  subspecies[] <- "Antarctic" #start by assuming A
  
  #Catches around Durban, which are unknown
  subspecies[latlon$Lon>=20 & latlon$Lon<39 & 
               latlon$Lat> -36 & latlon$Lat<= 17] <- "Unknown Durban"
  
  #assign catches to pygmy blue whales
  subspecies[latlon$Lon>=20 & latlon$Lon<30 & 
               latlon$Lat>=-46 & latlon$Lat<= -36] <- "pygmy"
  subspecies[latlon$Lon>=30 & latlon$Lon<39 & 
               latlon$Lat>=-52 & latlon$Lat<= -36] <- "pygmy"
  subspecies[latlon$Lon>=39 & latlon$Lon<70 & 
               latlon$Lat>=-52 & latlon$Lat<= 30] <- "pygmy"
  subspecies[latlon$Lon>=70 & latlon$Lon<80 & 
               latlon$Lat>=-53 & latlon$Lat<= 30] <- "pygmy"
  subspecies[latlon$Lon>=80 & latlon$Lon<100 & 
               latlon$Lat>=-52 & latlon$Lat<= 30] <- "pygmy"
  subspecies[latlon$Lon>=100 & latlon$Lon<=180 & 
               latlon$Lat>=-52 & latlon$Lat<= 0] <- "pygmy"
  
  #61 catches that were coded as pygmy blue whales between 53-57S and 44-87E
  #were coded as pygmy in the IWC catch database and are assumed here to be pygmy
  #but five coded as pygmy south of 60S were assumed to be Antarctic
  #subspecies[latlon$Lon>=40 & latlon$Lon<=90 & 
               #latlon$Lat > -58 & latlon$Lat <= -54 & 
               #bluedata$Sp==15] <- "pygmy"
  
  #SE PACIFIC  
  subspecies[latlon$Lon>=-120 & latlon$Lon<=-69 & 
               latlon$Lat>=-50 & latlon$Lat<= 2] <- "SE Pacific"
  
  subspecies[latlon$Lon>=-180 & latlon$Lon<=-80 & 
               latlon$Lat>2 & latlon$Lat<= 89] <- "North Pacific"
  subspecies[latlon$Lon>=100 & latlon$Lon<=180 & 
               latlon$Lat>2 & latlon$Lat<= 89] <- "North Pacific"
  
  subspecies[latlon$Lon>=-80 & latlon$Lon<=30 & 
               latlon$Lat>=0 & latlon$Lat<= 89] <- "North Atlantic"
  subspecies[latlon$Lon>=30 & latlon$Lon<=100 & 
               latlon$Lat>=35 & latlon$Lat<= 89] <- "North Atlantic"
  
  subspecies[latlon$Lon==0 | latlon$Lat==0] <- "Unknown"
  subspecies[is.na(latlon$Lon) | is.na(latlon$Lat)] <- "Unknown"
  
  #coded as pygmy blue whales but longer than 80.5 ft
  #are re-coded as Antarctic
  #DecFeet <- convert.to.decfeet(lenvec=bluedata$Len, codevec=bluedata$L.u)
  #subspecies[subspecies=="pygmy" & DecFeet >= 80.01] <- "Antarctic"
  
  invisible(subspecies)   
}


#Catch data using tidyverse
is.blue.whale.z <- function(sppvec) {
  return(sppvec %in% c(4, 15))
}
#spp <- is.blue.whale(sppvec=XX$Sp)
#table(spp)
#table(XX$Sp)

#=========get.blue.data========================
get.blue.data.z <- function(files) {
  XX <- read_csv(file=files[1])[,-40] #last column was causing issues and is just notes so took it out
  if (length(files) > 1) {
    for (i in 2:length(files)) { 
      XX <- rbind(XX,read_csv(file=files[i])[,-40])
    }
  }
  #extract just blue whales
  bluedata <- XX[is.blue.whale.z(sppvec=XX$Sp),]
  
  return(bluedata)
}
