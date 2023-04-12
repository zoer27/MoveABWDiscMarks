////Integrated Discovery Mark Model--SIMULATION VERSION
//Negative binomial likelihood for marks, fixed survival, estimated r
//explicit incorporation of survival in the popualtion model
//includes  abundance estimates and q parameters for Indian and Pacific estimated with MLE methods
//K in each basin determined by the stationary distribution of m
//movement between basins is constrained
//MSY paramterization of theta-logistic
//No partial surveys for simulations, so no need for q esimates and all abundance estimates
//are included in the variables called "Br" (no separate JARPA info)
//Zoe Rand 
//12/9/22 
data {
  int<lower=0> ABr; //number data points for abundance, Branch
  int<lower=0> Nbasin; //number of basins
  int<lower = 0> Ibasin[Nbasin]; //index for basins
  int<lower=0> NyearAbund; //number of years for abundance model(total years of trajectory)
  int<lower=0> Nyearcatch; //number of years with catch data
  int<lower=0> Nyeartag; //number of years with tagging data (recoveries)
  vector[ABr] AEstBr; //Abundance Estimates, Branch
  vector[ABr] ACVBr; //Abundance CV, Branch
  int<lower=0> NAYearBrA; //Number of years with abundance estimates in each Atlantic 
  int<lower = 0> NAYearBrI;//Number of years with abundance estimates in each Indian
  int<lower = 0> NAYearBrP; //Number of years with abundance estimates in each Pacific
  int<lower = 0> AYear[NAYearBrA]; //index indicating the years of the trajectory that have Abudance estiamtes in Atlantic
  int<lower = 0> IYearBr[NAYearBrI]; //index indicating the years of the trajectory that have Abudance estiamtes in Indian, Branch
  int<lower = 0> PYearBr[NAYearBrP]; //index indicating the years of the trajectory that have Abudance estiamtes in Pacific, Branch
  matrix[Nyearcatch, Nbasin] Catch;//Matrix of catches in each basin and year
  int<lower = 0> Rec[Nyeartag-1, Nbasin*Nbasin]; //recoveries, needs to be an array for negative binomial likelihood
  matrix[Nyeartag, Nbasin] Rel; //releases
  real ub; //upper bound for movement probability
  real s;
}

transformed data{
  int<lower = 0> NYNocatch = NyearAbund - Nyearcatch;
  matrix[NyearAbund, Nbasin] CatchAll;
  matrix[NYNocatch, Nbasin] NoCatch;
  
  
  NoCatch = rep_matrix(0,NYNocatch, Nbasin);
  CatchAll = append_row(Catch, NoCatch); //adding rows of 0 catch
  

}
parameters {
  real<lower = 0> r; //r for abundance model, fixing this for now
  real<lower = 0> lnMSY; // MSY for Antarctic
  real<lower = 0, upper = ub> MAB; //movement from A to B
  real<lower = 0, upper = ub> MAC; //movement from A to C
  real <lower = 0, upper = ub> MBA; //movement from B to A
  real<lower = 0, upper = ub> MBC; //movement from B to C
  real<lower = 0, upper = ub> MCA; //movement from C to A
  real<lower = 0, upper = ub> MCB; //movement from C to B
  real<lower = 0, upper = 1> tl; //tag loss parameter
  real<lower = 0.001> inv_theta; //overdispersion--but actually 1/sqrt(theta), the bound makes it halfnormal
}

transformed parameters{
  real lnK; //K in log space for prior
  real <lower = 0> K[Nbasin];//carrying capacity for each basin antarctic
  real<lower = 0> MSY[Nbasin]; //MSY for each basin
  real<lower = 0> MAA; //stay in A
  real<lower = 0> MBB; //stay in B
  real<lower = 0> MCC; //stay in C
  matrix<lower = 0>[Nbasin,Nbasin] m; //converting parameter to matrix form for use in model
  matrix[Nbasin,Nbasin] stat_m; //equilibrium matrix
  row_vector[Nbasin] dist; //equilibirium proportions 
  real<lower = 0> theta; //transformation of overdispersion
  real r_star; //transformed r for explicitly incorporating survival
  vector[ABr] sig; //vector for sigma values for likelihood
  
  
  //transformation into r_star
  r_star = r/(1-s);
  
   //by subtraction
  MAA = 1 - MAB - MAC;
  MBB = 1 - MBA - MBC;
  MCC = 1 - MCA - MCB;
  
 //matrix
 m[1,1] = MAA;
 m[1,2] = MAB;
 m[1,3] = MAC;
 m[2,1] = MBA;
 m[2,2] = MBB;
 m[2,3] = MBC;
 m[3,1] = MCA;
 m[3,2] = MCB;
 m[3,3] = MCC;
 
 
  stat_m = matrix_power(m,100);
  dist = stat_m[1,]; //one row of stationary dist

  

  lnK = log(((3.39)^((1/2.39) + 1)) *exp(lnMSY)/((1-s)*r_star*2.39));

  for(n in 1:Nbasin){
     K[n] = exp(lnK)*dist[n]; //carrying capacity from equilibrium movement structure
     MSY[n] = ((1-s)*r_star*K[n]*2.39)/((3.39)^((1/2.39) + 1));
  }
  //converting to sigma for lognormal from cvs 
    for(i in 1:(ABr)){
      sig[i] = sqrt(log((ACVBr[i])^2 + 1)); 
    }
    
      //conversion of theta for negative binomial
  theta = 1/(inv_theta^2);
  
}

model {
  matrix[NyearAbund, Nbasin] Npop;//population 
  vector[ABr] PredAbundBr; //vector for subsetting population predictions (Branch)
  vector[NAYearBrA] PredAbundAt; // population predictions for Atlantic
  vector[NAYearBrI] PredAbundInd; //population predictions for Indian
  vector[NAYearBrP] PredAbundPac; //population predictions for Pacific
  matrix[Nyeartag, Nbasin*Nbasin] NT; //predicted number of tags
  row_vector[Nbasin] new_tags; //vector of new tags for each year
  matrix[Nyeartag, Nbasin] h; //harvest rates for remaining in population
  matrix[(Nyeartag-1), Nbasin] PredH; //harvest rates for recovery
  matrix[(Nyeartag-1), Nbasin*Nbasin] NTR; //predicted number of tags recovered
  matrix[(Nyeartag-1), Nbasin*Nbasin] PredRec;//predicted recoveries that match up with  Data

  
  //priors
  r ~ uniform(0, 0.114); //prior for r, based on Monahan paper
  lnK ~ uniform(9, 13);
  tl ~ beta(1,1);
  MAB ~ uniform(0,ub);
  MAC ~ uniform(0, ub);
  MBA ~ uniform(0, ub);
  MBC ~ uniform(0, ub);
  MCA ~ uniform(0, ub);
  MCB ~ uniform(0,ub);
  
  inv_theta ~ normal(0,1); //prior for overdisperion (from prior choice recommendations on stan-dev wiki)

  //model predictions
    //pop model
   for(I in Ibasin){
     Npop[1,I] = K[I]; //first year is carrying capacity
     //print(Npop[1,I]);
    } //pop model explicitly includes survival
    for (i in 2:NyearAbund){
      for(I in Ibasin){
    Npop[i,I] = s*Npop[i-1,I] + (1-s)*Npop[i-1,I]*(1+r_star*(1-(((1-s)*r_star*2.39*Npop[i-1,I])/(MSY[I]*(3.39^(1/2.39 +1))))^2.39));
    if(CatchAll[i-1,I] <= (Npop[i,I])){ //limiting harvest rate so can't be greater than max set
      Npop[i,I] = Npop[i,I] - CatchAll[i-1,I];
     
    } else{
      Npop[i,I] = 1; //if catch is greater than abundance turns pop into 1
    }
    Npop[i,I] = fmax(Npop[i,I], 1); //if population goes below 1 this becomes 1
      }
    Npop[i,] = Npop[i,] * m;
    }
    
    //harvest rates
    for (i in 1:(Nyeartag)){ //harvest rates start in 1926 and ends in 1972
      for (I in Ibasin){
        h[i,I] = CatchAll[i+13,I]/Npop[i+13,I];
        h[i,I] = fmin(h[i,I], 1); //harvest rate can't be more than 1 (if pop is 1 then all whales are harvested and harvest rate is 1)
        }
      }
      
    PredH = h[2:Nyeartag,]; //no harvest rates before 1927

    //mark recovery model
    NT[1,1] = Rel[1,1]; //first year releases
    NT[1, 2:4] = rep_row_vector(0,3);
    NT[1,5] = Rel[1,2];
    NT[1,6:8] = rep_row_vector(0,3);
    NT[1,9] = Rel[1,3];
    
    for(i in 2:Nyeartag){ //recoveries for all years
      //print(i);
        new_tags = Rel[i,];
        NT[i,1] = (NT[i-1, 1]*(1-h[i-1,1])*s*m[1,1]) + (NT[i-1,2]*(1-h[i-1,2])*s*m[2,1]) + (NT[i-1,3]*(1-h[i-1,3])*s*m[3,1]);
        NT[i,2] = (NT[i-1, 1]*(1-h[i-1,1])*s*m[1,2]) + (NT[i-1,2]*(1-h[i-1,2])*s*m[2,2]) + (NT[i-1,3]*(1-h[i-1,3])*s*m[3,2]);
        NT[i,3] = (NT[i-1, 1]*(1-h[i-1,1])*s*m[1,3]) + (NT[i-1,2]*(1-h[i-1,2])*s*m[2,3]) + (NT[i-1,3]*(1-h[i-1,3])*s*m[3,3]);
        NT[i,4] = (NT[i-1, 4]*(1-h[i-1,1])*s*m[1,1]) + (NT[i-1,5]*(1-h[i-1,2])*s*m[2,1]) + (NT[i-1,6]*(1-h[i-1,3])*s*m[3,1]);
        NT[i,5] = (NT[i-1, 4]*(1-h[i-1,1])*s*m[1,2]) + (NT[i-1,5]*(1-h[i-1,2])*s*m[2,2]) + (NT[i-1,6]*(1-h[i-1,3])*s*m[3,2]);
        NT[i,6] = (NT[i-1, 4]*(1-h[i-1,1])*s*m[1,3]) + (NT[i-1,5]*(1-h[i-1,2])*s*m[2,3]) + (NT[i-1,6]*(1-h[i-1,3])*s*m[3,3]);
        NT[i,7] = (NT[i-1, 7]*(1-h[i-1,1])*s*m[1,1]) + (NT[i-1,8]*(1-h[i-1,2])*s*m[2,1]) + (NT[i-1,9]*(1-h[i-1,3])*s*m[3,1]);
        NT[i,8] = (NT[i-1, 7]*(1-h[i-1,1])*s*m[1,2]) + (NT[i-1,8]*(1-h[i-1,2])*s*m[2,2]) + (NT[i-1,9]*(1-h[i-1,3])*s*m[3,2]);
        NT[i,9] = (NT[i-1, 7]*(1-h[i-1,1])*s*m[1,3]) + (NT[i-1,8]*(1-h[i-1,2])*s*m[2,3]) + (NT[i-1,9]*(1-h[i-1,3])*s*m[3,3]);
      
        for(I in Ibasin){
          NTR[i-1,I] = NT[i,I] * PredH[i-1, I];
          NTR[i-1, I+3] = NT[i,I + 3] * PredH[i-1, I];
          NTR[i-1, I + 6] = NT[i,I + 3] * PredH[i-1, I];
        }
       NT[i,1] = NT[i, 1] + new_tags[1];
       NT[i,5] = NT[i, 5] + new_tags[2];
       NT[i, 9] = NT[i,9] + new_tags[3];
      }
    
    NTR = NTR*(1-tl); //recovered tags are not lost
    
    for(i in 1:(Nyeartag-1)){
      for(k in 1:9){
        NTR[i,k] = fmax(NTR[i,k], 0.001);//if 0 is in predictions replace with 0.001 so likelihood works
        }
    }

    
     //getting abundance predictions for same year and basins as data--all
     
    for(i in 1:NAYearBrA){
      PredAbundAt[i] = Npop[AYear[i], Ibasin[1]];
    }
    for(k in 1:NAYearBrI){
      PredAbundInd[k] = Npop[IYearBr[k], Ibasin[2]];
    }
    for(j in 1:NAYearBrP){
      PredAbundPac[j] = Npop[PYearBr[j], Ibasin[3]];
    }
    
    PredAbundBr = append_row(append_row(PredAbundAt, PredAbundInd), PredAbundPac);
    
  
    //likelihood for abundance 
    
    //lognormal likelihood -- Branch 
    AEstBr ~ lognormal(log(PredAbundBr), sig[1:ABr]);
    
  
   //likelihood for tags
   for (i in 1:(Nyeartag-1)){
     Rec[i,] ~ neg_binomial_2(NTR[i,], theta);
   }
}

generated quantities{
  matrix[NyearAbund, Nbasin] Npop;//population
  matrix[(Nyeartag-1), Nbasin*Nbasin] NTR;//predicted tags recovered
  matrix[Nyeartag, Nbasin] h; //harvest rates

   {
  
    matrix[Nyeartag, Nbasin*Nbasin] NT; //predicted number of tags
    row_vector[Nbasin] new_tags; //vector of new tags for each year
    matrix[(Nyeartag-1), Nbasin] PredH;
    matrix[(Nyeartag-1), Nbasin*Nbasin] PredRec; //predicted recoveries that match up with data
     //pop model
   for(I in Ibasin){
     Npop[1,I] = K[I]; //first year is carrying capacity
     //print(Npop[1,I]);
    }
    for (i in 2:NyearAbund){
      for(I in Ibasin){ //pop model explicitly incorporates survival
     Npop[i,I] = s*Npop[i-1,I] + (1-s)*Npop[i-1,I]*(1+r_star*(1-(((1-s)*r_star*2.39*Npop[i-1,I])/(MSY[I]*(3.39^(1/2.39 +1))))^2.39));
    if(CatchAll[i-1,I] <= (Npop[i,I])){ //limiting harvest rate so can't be greater than max set
      Npop[i,I] = Npop[i,I] - CatchAll[i-1,I];
     
    } else{
      Npop[i,I] = 1; //if catch is greater than abundance turns pop into 1
    }
    Npop[i,I] = fmax(Npop[i,I], 1); //if population goes below 1 this becomes 1
      }
    Npop[i,] = Npop[i,] * m;
    }
  
  //tag recovery model
  
  
    //harvest rates
    for (i in 1:(Nyeartag)){ //harvest rates start in 1926 and ends in 1972
      for (I in Ibasin){
        h[i,I] = CatchAll[i+13,I]/Npop[i+13,I];
        h[i,I] = fmin(h[i,I], 1); //harvest rate can't be more than 1 (if pop is 1 then all whales are harvested and harvest rate is 1)
        }
    }
      
    PredH = h[2:Nyeartag,]; //no harvest rates before 1927
   //mark recovery model
    NT[1,1] = Rel[1,1]; //first year releases
    NT[1, 2:4] = rep_row_vector(0,3);
    NT[1,5] = Rel[1,2];
    NT[1,6:8] = rep_row_vector(0,3);
    NT[1,9] = Rel[1,3];
    
   for(i in 2:Nyeartag){ //recoveries for all years
      //print(i);
        new_tags = Rel[i,];
        NT[i,1] = (NT[i-1, 1]*(1-h[i-1,1])*s*m[1,1]) + (NT[i-1,2]*(1-h[i-1,2])*s*m[2,1]) + (NT[i-1,3]*(1-h[i-1,3])*s*m[3,1]);
        NT[i,2] = (NT[i-1, 1]*(1-h[i-1,1])*s*m[1,2]) + (NT[i-1,2]*(1-h[i-1,2])*s*m[2,2]) + (NT[i-1,3]*(1-h[i-1,3])*s*m[3,2]);
        NT[i,3] = (NT[i-1, 1]*(1-h[i-1,1])*s*m[1,3]) + (NT[i-1,2]*(1-h[i-1,2])*s*m[2,3]) + (NT[i-1,3]*(1-h[i-1,3])*s*m[3,3]);
        NT[i,4] = (NT[i-1, 4]*(1-h[i-1,1])*s*m[1,1]) + (NT[i-1,5]*(1-h[i-1,2])*s*m[2,1]) + (NT[i-1,6]*(1-h[i-1,3])*s*m[3,1]);
        NT[i,5] = (NT[i-1, 4]*(1-h[i-1,1])*s*m[1,2]) + (NT[i-1,5]*(1-h[i-1,2])*s*m[2,2]) + (NT[i-1,6]*(1-h[i-1,3])*s*m[3,2]);
        NT[i,6] = (NT[i-1, 4]*(1-h[i-1,1])*s*m[1,3]) + (NT[i-1,5]*(1-h[i-1,2])*s*m[2,3]) + (NT[i-1,6]*(1-h[i-1,3])*s*m[3,3]);
        NT[i,7] = (NT[i-1, 7]*(1-h[i-1,1])*s*m[1,1]) + (NT[i-1,8]*(1-h[i-1,2])*s*m[2,1]) + (NT[i-1,9]*(1-h[i-1,3])*s*m[3,1]);
        NT[i,8] = (NT[i-1, 7]*(1-h[i-1,1])*s*m[1,2]) + (NT[i-1,8]*(1-h[i-1,2])*s*m[2,2]) + (NT[i-1,9]*(1-h[i-1,3])*s*m[3,2]);
        NT[i,9] = (NT[i-1, 7]*(1-h[i-1,1])*s*m[1,3]) + (NT[i-1,8]*(1-h[i-1,2])*s*m[2,3]) + (NT[i-1,9]*(1-h[i-1,3])*s*m[3,3]);
        for(I in Ibasin){
          NTR[i-1,I] = NT[i,I] * PredH[i-1, I];
          NTR[i-1, I+3] = NT[i,I + 3] * PredH[i-1, I];
          NTR[i-1, I + 6] = NT[i,I + 3] * PredH[i-1, I];
        }
       NT[i,1] = NT[i, 1] + new_tags[1];
       NT[i,5] = NT[i, 5] + new_tags[2];
       NT[i, 9] = NT[i,9] + new_tags[3];
      }
    
    NTR = NTR*(1-tl); //recovered tags are not lost
    
    for(i in 1:(Nyeartag-1)){
      for(k in 1:9){
        NTR[i,k] = fmax(NTR[i,k], 0.001);//if 0 is in predictions replace with 0.001 so likelihood works
        }
    }
    
   }
   
  }

