### Data Manipulation Functions for Custom Mark Recapture Models
#Zoe Rand

#%!in% finds values that are not equal to any in a vector
'%!in%' <- function(x,y)!('%in%'(x,y))

#get_direction() takes the mark data that is missing the E/W, N/S columns,
#and the MarkMaster data that has it and adds in the columns to the mark data
get_direction<-function(missing_df, main_df){
  for (i in 1:nrow(missing_df)){
    t<-which(missing_df$MarkID[i] == main_df$"No.")
    missing_df$"N/S"[i]<-main_df$X8[t]
    missing_df$"E/W"[i]<-main_df$X10[t]
  }
  return(missing_df)
}

#adapted from Trevor's Useful Functions
get.subspecies.withlatlon <- function(Lon, Lat) {
  latlon <- data.frame(Lon = Lon, Lat = Lat)
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
  
  invisible(subspecies)   
}


#isEmpty
#tests whether a vector is empty (for use in other functions)
isEmpty<-function(x){
  return(identical(x, numeric(0)))
}

#stay fills in the probability of staying a location as 1-probability of moving
stay<-function(move_mat){
  #fills in probability of staying with 1-prob of moving
  totals<-rowSums(move_mat, na.rm = TRUE)
  for (i in 1:nrow(move_mat)){
    for (j in 1:ncol(move_mat)){
      if (is.na(move_mat[i,j])){
        move_mat[i,j]<-(1-totals[i])
      }
    }
  }
  return(move_mat)
}
