////Updated Model fitting to actual Abundance Data--New groups Option 1
//Poisson likelihood for tags, fixed survival, estimated r
//explicit incorporation of survival in the popualtion model
//includes Matsouka abundance estimates and q parameters for Indian and Pacific estimated with MLE methods
//K in each basin determined by the stationary distribution of m
//movement between basins is constrained
//MSY paramterization of theta-logistic
//Updated by Zoe Rand 5/11/22
//Updated 10/12/22 to include posterior predictive sampling and linear algebra
//Updated by Zoe Rand 10/20/22 -- fixed MSY mistake, and changed location of new tag releases in mark model
//Updated ZR 11/7/22 to use prior on K instead of MSY. looks like should be a linear transformation so do not need jacobian
//updated 11/8/22 to use informative prior for S instead of fixing it
//updated 11/21/22 back to fixed s, and includes log likelihood generated quantities for model selection
//11/29/22 added negative binomial back in for mark recovery data
//8/21/23 adding option so I can re-run this with new groups
//9/20/23 splitting new groups into two different stan files
data {
  int<lower=0> ABr; //number data points for abundance, Branch
  int<lower = 0> AMat; //number of data points for Indian abdundance, Matsouka
  int<lower = 0> AMat2; //number of data points in Pacific basin Matsouka
  int<lower = 0> AMat3; //1 if option 1, 0 if otherwise, to account for the 0 that isn't in the Indian
  int<lower=0> Nbasin; //number of basins
  int<lower = 0> Ibasin[Nbasin]; //index for basins
  int<lower=0> NyearAbund; //number of years for abundance model(total years of trajectory)
  int<lower=0> Nyearcatch; //number of years with catch data
  int<lower=0> Nyeartag; //number of years with tagging data (recoveries)
  vector[ABr] AEstBr; //Abundance Estimates, Branch
  vector[ABr] ACVBr; //Abundance CV, Branch
  int<lower = 0> NCV; //Total number of CVs in Matusoka/Hamabe
  vector[AMat] AEstMatInd; //Abundance Estimates Matsouka--Indian
  vector[AMat2] AEstMatPac; //Abundance Estimates Matsouka -- Pacific
  vector[NCV] ACVMat; //Abundance CV Matsouka
  array[AMat] int IndCVIdx; //Index for CVs from indian--remove CV from 0
  array[AMat2] int PacCVIdx; //Index for CVs from Pacific
  array[AMat3] int ZeroCVIdx; //Index for CVs from Zero estimate
  int<lower=0> NAYearBr; //Number of years with abundance estimates in each basin (from branch)
  int<lower=0> NAYearMatI; //Number of years with abundance estimates in each basin (from Matsouka)
  int<lower=0> NAYearMatP; //Number of years with abundance estimates in each basin (from Matsouka)
  int<lower = 0> AYear[NAYearBr]; //index indicating the years of the trajectory that have Abudance estiamtes in Atlantic
  int<lower = 0> IYearBr[NAYearBr]; //index indicating the years of the trajectory that have Abudance estiamtes in Indian, Branch
  int<lower = 0> PYearBr[NAYearBr]; //index indicating the years of the trajectory that have Abudance estiamtes in Pacific, Branch
  int<lower = 0>IYearMat[NAYearMatI];//index indicating the years of the trajectory that have Abudance estiamtes in Indian, Matsouka
  int<lower=0>PYearMat[NAYearMatP]; //index indicating the years of the trajectory that have Abudance estiamtes in Pacific, Matsouka
  int<lower=0>NormYear[AMat3]; //index for year that should use a normal distribution in option 1
  matrix[Nyearcatch, Nbasin] Catch;//Matrix of catches in each basin and year
  int<lower = 0> Rec[Nyeartag-1, Nbasin*Nbasin]; //recoveries, needs to be an array for negative binomial likelihood
  matrix[Nyeartag, Nbasin] Rel; //releases
  real s; //fixed s
  real ub; //upper bound for movement probability
}

transformed data{
  int<lower = 0> NYNocatch = NyearAbund - Nyearcatch;
  matrix[NyearAbund, Nbasin] CatchAll;
  matrix[NYNocatch, Nbasin] NoCatch;
  vector[ABr + NCV] ACV; //all the abundance estimate CVs
  
  ACV = append_row(ACVBr, ACVMat); 
  //real r_star; //transformed r for explicitly incorporating survival
  
  NoCatch = rep_matrix(0,NYNocatch, Nbasin);
  CatchAll = append_row(Catch, NoCatch); //adding rows of 0 catch
  

}
parameters {
  real<lower = 0> r; //r for abundance model, fixing this for now
  real lnMSY; // MSY for Antarctic
  //real logAddCV; //additional CV for abundance (in log space), not including right now
  real<lower = 0, upper = ub> MAB; //movement from A to B
  real<lower = 0, upper = ub> MAC; //movement from A to C
  real <lower = 0, upper = ub> MBA; //movement from B to A
  real<lower = 0, upper = ub> MBC; //movement from B to C
  real<lower = 0, upper = ub> MCA; //movement from C to A
  real<lower = 0, upper = ub> MCB; //movement from C to B
  real<lower = 0, upper = 1> tl; //tag loss parameter
  real<lower = 0.001> inv_theta; //overdispersion--but actually 1/sqrt(theta), the bound makes it halfnormal
  //real<lower = 0.93, upper = 0.99> s; //natural survival truncated normal prior
  //real<lower = 0, upper = 1> q[2]; //adding in q parameter as parameter instead of calculation
}

transformed parameters{
  //real<lower = 0> AddCV; //additional CV for abundance
  //real <lower=0> MSY_SO; //MSY for Southern Ocean
  //real<lower = 0> K_SO; //K for soutehrn ocean
  real lnK; //K in log space for prior
  real <lower = 0> K[Nbasin];//carrying capacity for each basin antarctic
  real<lower = 0> MSY[Nbasin]; //MSY for each basin
  real<lower = 0> MAA; //stay in A
  real<lower = 0> MBB; //stay in B
  real<lower = 0> MCC; //stay in C
  matrix<lower = 0>[Nbasin,Nbasin] m; //converting parameter to matrix form for use in model
  matrix[Nbasin,Nbasin] stat_m; //equilibrium matrix
  row_vector[Nbasin] dist; //equilibirium proportions 
  //complex_matrix[Nbasin, Nbasin] stat_m; //equilibrium distribution
  //vector[Nbasin] dist; //equilibirium proportions 
  real<lower = 0> theta; //transformation of overdispersion
  real r_star; //transformed r for explicitly incorporating survival
  vector[ABr+NCV] sig; //vector for sigma values for likelihood
  
  
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
 
 //stationary distribution using generalized inverse function
 //{complex_matrix[Nbasin, Nbasin] left_eig;
  
  
    //left_eig = eigenvectors(m'); //stan gives right (column) eigen vectors so need left eigen vectors by using transpose of the matrix
    //print(left_eig);
    //stat_m = left_eig./sum(left_eig); //normalize so sum to 1
 //} 
 //stat_m = matrix_power(m, 100); //old version
 //dist = get_real(stat_m[,1]); 
  stat_m = matrix_power(m,100);
  dist = stat_m[1,]; //one row of stationary dist

  //MSY in each basin
  //MSY_SO = exp(lnMSY);

  lnK = log(((3.39)^((1/2.39) + 1)) *exp(lnMSY)/((1-s)*r_star*2.39));
  //lnK = log(K_SO);
  for(n in 1:Nbasin){
     K[n] = exp(lnK)*dist[n]; //carrying capacity from equilibrium movement structure
     MSY[n] = ((1-s)*r_star*K[n]*2.39)/((3.39)^((1/2.39) + 1));
  }
  //converting to sigma for lognormal from cvs (for both Mat and Branch)
    for(i in 1:(ABr+NCV)){
      sig[i] = sqrt(log((ACV[i])^2 + 1)); 
    }
    //conversion of theta for negative binomial
  theta = 1/(inv_theta^2);
}

model {
  matrix[NyearAbund, Nbasin] Npop;//population 
  vector[ABr] PredAbundBr; //vector for subsetting population predictions (Branch)
  vector[AMat] PredAbundMatInd; //vector for subsetting population predictions, Indian (Matsouka)
  //vector[AMat] PredAbundMatInd2; //vector including the 0 (group B)
  vector[AMat2] PredAbundMatPac; //vector for population predictions, Pacific, (Matsouka)
  vector[1]PredAbundNorm; //for one year that has a normal distribution in option 1
  vector[2] q; //q for Matsouka abundance, one for each basin (Indian and Pacific)
  matrix[Nyeartag, Nbasin*Nbasin] NT; //predicted number of tags
  row_vector[Nbasin] new_tags; //vector of new tags for each year
  matrix[Nyeartag, Nbasin] h; //harvest rates for remaining in population
  matrix[(Nyeartag-1), Nbasin] PredH; //harvest rates for recovery
  matrix[(Nyeartag-1), Nbasin*Nbasin] NTR; //predicted number of tags recovered
  matrix[(Nyeartag-1), Nbasin*Nbasin] PredRec;//predicted recoveries that match up with  Data

  
  //priors
  r ~ uniform(0, 0.114); //prior for r, based on Monahan paper
  lnK ~ uniform(9, 13);
  //lnMSY~ uniform(0,10); //based on r and previous prior for K
   // lower bound is minimum number needed to not go extinct based on catch series
  //AddCV ~ uniform(0,1); 
  tl ~ beta(1,1);
  MAB ~ uniform(0,ub);
  MAC ~ uniform(0, ub);
  MBA ~ uniform(0, ub);
  MBC ~ uniform(0, ub);
  MCA ~ uniform(0, ub);
  MCB ~ uniform(0,ub);
  //for(i in 1:2){
    //q[i] ~ beta(1,1);
  //}
  
  //s ~ normal(0.96, 0.02);
  
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
        h[i,I] = CatchAll[i+22,I]/Npop[i+22,I];
        h[i,I] = fmin(h[i,I], 0.999); //harvest rate can't be more than 1 (if pop is 1 then all whales are harvested and harvest rate is 1)
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
        //NT[i,1] = NT[i, 1] + new_tags[1];
        NT[i,2] = (NT[i-1, 1]*(1-h[i-1,1])*s*m[1,2]) + (NT[i-1,2]*(1-h[i-1,2])*s*m[2,2]) + (NT[i-1,3]*(1-h[i-1,3])*s*m[3,2]);
        NT[i,3] = (NT[i-1, 1]*(1-h[i-1,1])*s*m[1,3]) + (NT[i-1,2]*(1-h[i-1,2])*s*m[2,3]) + (NT[i-1,3]*(1-h[i-1,3])*s*m[3,3]);
        NT[i,4] = (NT[i-1, 4]*(1-h[i-1,1])*s*m[1,1]) + (NT[i-1,5]*(1-h[i-1,2])*s*m[2,1]) + (NT[i-1,6]*(1-h[i-1,3])*s*m[3,1]);
        NT[i,5] = (NT[i-1, 4]*(1-h[i-1,1])*s*m[1,2]) + (NT[i-1,5]*(1-h[i-1,2])*s*m[2,2]) + (NT[i-1,6]*(1-h[i-1,3])*s*m[3,2]);
        //NT[i,5] = NT[i, 5] + new_tags[2];
        NT[i,6] = (NT[i-1, 4]*(1-h[i-1,1])*s*m[1,3]) + (NT[i-1,5]*(1-h[i-1,2])*s*m[2,3]) + (NT[i-1,6]*(1-h[i-1,3])*s*m[3,3]);
        NT[i,7] = (NT[i-1, 7]*(1-h[i-1,1])*s*m[1,1]) + (NT[i-1,8]*(1-h[i-1,2])*s*m[2,1]) + (NT[i-1,9]*(1-h[i-1,3])*s*m[3,1]);
        NT[i,8] = (NT[i-1, 7]*(1-h[i-1,1])*s*m[1,2]) + (NT[i-1,8]*(1-h[i-1,2])*s*m[2,2]) + (NT[i-1,9]*(1-h[i-1,3])*s*m[3,2]);
        NT[i,9] = (NT[i-1, 7]*(1-h[i-1,1])*s*m[1,3]) + (NT[i-1,8]*(1-h[i-1,2])*s*m[2,3]) + (NT[i-1,9]*(1-h[i-1,3])*s*m[3,3]);
        //NT[i, 9] = NT[i,9] + new_tags[3];
        for(I in Ibasin){
          NTR[i-1,I] = NT[i,I] * PredH[i-1, I];
          NTR[i-1, I+3] = NT[i,I + 3] * PredH[i-1, I];
          NTR[i-1, I + 6] = NT[i,I + 3] * PredH[i-1, I];
        }
       NT[i,1] = NT[i, 1] + new_tags[1];
       NT[i,5] = NT[i, 5] + new_tags[2];
       NT[i, 9] = NT[i,9] + new_tags[3];
      }
    //PredRec = NT[2:Nyeartag,]; //no first year recoveries
    //PredH = h[2:Nyeartag,]; //no harvest rates before 1927
    //for (I in Ibasin){
     //NTR[,I] = PredRec[,I] .* PredH[,I]; //recovered tags are harvested
     //NTR[,I+3] = PredRec[,I+3] .* PredH[,I];
     //NTR[,I+6] = PredRec[,I+6] .* PredH[,I];
   //}
    
    NTR = NTR*(1-tl); //recovered tags are not lost
    
    for(i in 1:(Nyeartag-1)){
      for(k in 1:9){
        NTR[i,k] = fmax(NTR[i,k], 0.001);//if 0 is in predictions replace with 0.001 so likelihood works
        }
    }

    
     //getting abundance predictions for same year and basins as data--Branch
    
    for(i in 1:NAYearBr){
      PredAbundBr[i] = Npop[AYear[i], Ibasin[1]];
    }
    for(k in (NAYearBr + 1):(NAYearBr + NAYearBr)){
      PredAbundBr[k] = Npop[IYearBr[k-3], Ibasin[2]];
    }
    for(j in (NAYearBr + NAYearBr + 1): (NAYearBr + NAYearBr + NAYearBr)){
      PredAbundBr[j] = Npop[PYearBr[j-6], Ibasin[3]];
    }
    
    //getting abundance predictions Matsouka
    for(i in 1:(NAYearMatI)){
      PredAbundMatInd[i] = Npop[IYearMat[i], Ibasin[2]];
    }
      PredAbundNorm = Npop[NormYear, Ibasin[2]]; //get the prediction for the 0
      //PredAbundMatInd2 = append_row(PredAbundMatInd, PredAbundNorm); //adding 0 to end for q
      
    for(i in 1:NAYearMatP){
      PredAbundMatPac[i] = Npop[PYearMat[i], Ibasin[3]];
      }
    
   
    
    
    //likelihood for abundance 
    
    //lognormal likelihood -- Branch 
    AEstBr ~ lognormal(log(PredAbundBr), sig[1:ABr]);
    //print(PredAbundMatInd);
    //qs for Mat abundance estimates --from MLE
      q[1] = exp(mean(log(AEstMatInd./PredAbundMatInd))); //Indian
      q[2] = exp(mean(log(AEstMatPac./PredAbundMatPac))); //Pacific
      
      //q[1] = fmax(q[1], 0.0001);
      //print(q);
      q[1] = fmin(q[1], 1);//if q > 1, brings it back down to 1 -should make this not fit the data well and move away from this space
      q[2] = fmin(q[2], 1);
      //print(q);
    //lognormal likelihood--Matsouka
    AEstMatInd ~ lognormal(log(q[1] * PredAbundMatInd), sig[IndCVIdx]); //Indian not including last one
    0 ~ normal(q[1]*PredAbundNorm, sig[ZeroCVIdx]); //for 0 in 1995
    AEstMatPac ~ lognormal(log(q[2]*PredAbundMatPac), sig[PacCVIdx]);//Pacific

   //likelihood for tags
   for (i in 1:(Nyeartag-1)){
     //Rec[i,] ~ poisson(NTR[i,]);
     Rec[i,] ~ neg_binomial_2(NTR[i,], theta);
   }
}

generated quantities{
  matrix[NyearAbund, Nbasin] Npop;//population
  matrix[(Nyeartag-1), Nbasin*Nbasin] NTR;//predicted tags recovered
  matrix[Nyeartag, Nbasin] h; //harvest rates
  vector[2] q; //q parameter for Matsouka abundance estimates
  array[ABr] real PostPredBr; //Posterior predictive for Abundance Estimates, Branch
  array[AMat] real PostPredMatInd; //Posterior predictive Abundance Estimates Matsouka--Indian
  array[AMat3] real PostPredZero; //Posterior predictive for the zero abundance estimate
  array[AMat2] real PostPredMatPac; //Posterior predictive Abundance Estimates Matsouka -- Pacific
  array[Nyeartag-1, Nbasin*Nbasin] int<lower =0> PostPredRec; //Posterior Predictive mark recoveries
  vector[ABr+AMat + AMat2 +((Nyeartag -1)*Nbasin*Nbasin)] log_like; //vector of loglikelihood values for model selection, Rec values are created in a matrix and then returned in column-major order
   {
    vector[ABr] PredAbundBr; //vector for subsetting population predictions (Branch)
    vector[AMat] PredAbundMatInd; //vector for subsetting population predictions, Indian (Matsouka)
    vector[AMat2] PredAbundMatPac; //vector for population predictions, Pacific, (Matsouka)
    vector[AMat3] PredAbundNorm;
    vector[AMat] PredAbundMatInd2; //vector including the 0 (group B)
    matrix[Nyeartag, Nbasin*Nbasin] NT; //predicted number of tags
    row_vector[Nbasin] new_tags; //vector of new tags for each year
    matrix[(Nyeartag-1), Nbasin] PredH;
    matrix[(Nyeartag-1), Nbasin*Nbasin] PredRec; //predicted recoveries that match up with data
    vector[ABr] loglikeAbundBr; //vectors for storing individual log likelihoods
    vector[AMat2] loglikeAbundMatInd;
    //vector[AMat2] loglikeAbundMatPac;
    matrix[Nyeartag-1, Nbasin*Nbasin]loglikeRec;
    
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

     //getting abundance predictions for same year and basins as data--Branch
    
    for(i in 1:NAYearBr){
      PredAbundBr[i] = Npop[AYear[i], Ibasin[1]];
    }
    for(k in (NAYearBr + 1):(NAYearBr + NAYearBr)){
      PredAbundBr[k] = Npop[IYearBr[k-3], Ibasin[2]];
    }
    for(j in (NAYearBr + NAYearBr + 1): (NAYearBr + NAYearBr + NAYearBr)){
      PredAbundBr[j] = Npop[PYearBr[j-6], Ibasin[3]];
    }
    
    //getting abundance predictions Matsouka
    for(i in 1:(NAYearMatI)){
      PredAbundMatInd[i] = Npop[IYearMat[i], Ibasin[2]];
    }
      PredAbundNorm = Npop[NormYear, Ibasin[2]]; //get the prediction for the 0
      //PredAbundMatInd2 = append_row(PredAbundMatInd, PredAbundNorm); //adding 0 to end for q
      
    for(i in 1:NAYearMatP){
      PredAbundMatPac[i] = Npop[PYearMat[i], Ibasin[3]];
      }
   
    
    //calculating q
    q[1] = exp(mean(log(AEstMatInd./PredAbundMatInd)));
    q[2] = exp(mean(log(AEstMatPac./PredAbundMatPac))); //Pacific
    //q[1] = fmax(q[1], 0.000001);
    //print(q);
    q[1] = fmin(q[1], 1);//if q > 1, brings it back down to 1 -should make this not fit the data well and move away from this space
    q[2] = fmin(q[2], 1);
  //tag recovery model
  
  
    //harvest rates
    for (i in 1:(Nyeartag)){ //harvest rates start in 1926 and ends in 1972
      for (I in Ibasin){
        h[i,I] = CatchAll[i+22,I]/Npop[i+22,I];
        h[i,I] = fmin(h[i,I], 0.999); //harvest rate can't be more than 1 (if pop is 1 then all whales are harvested and harvest rate is 1)
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
        //NT[i,1] = NT[i, 1] + new_tags[1];
        NT[i,2] = (NT[i-1, 1]*(1-h[i-1,1])*s*m[1,2]) + (NT[i-1,2]*(1-h[i-1,2])*s*m[2,2]) + (NT[i-1,3]*(1-h[i-1,3])*s*m[3,2]);
        NT[i,3] = (NT[i-1, 1]*(1-h[i-1,1])*s*m[1,3]) + (NT[i-1,2]*(1-h[i-1,2])*s*m[2,3]) + (NT[i-1,3]*(1-h[i-1,3])*s*m[3,3]);
        NT[i,4] = (NT[i-1, 4]*(1-h[i-1,1])*s*m[1,1]) + (NT[i-1,5]*(1-h[i-1,2])*s*m[2,1]) + (NT[i-1,6]*(1-h[i-1,3])*s*m[3,1]);
        NT[i,5] = (NT[i-1, 4]*(1-h[i-1,1])*s*m[1,2]) + (NT[i-1,5]*(1-h[i-1,2])*s*m[2,2]) + (NT[i-1,6]*(1-h[i-1,3])*s*m[3,2]);
        //NT[i,5] = NT[i, 5] + new_tags[2];
        NT[i,6] = (NT[i-1, 4]*(1-h[i-1,1])*s*m[1,3]) + (NT[i-1,5]*(1-h[i-1,2])*s*m[2,3]) + (NT[i-1,6]*(1-h[i-1,3])*s*m[3,3]);
        NT[i,7] = (NT[i-1, 7]*(1-h[i-1,1])*s*m[1,1]) + (NT[i-1,8]*(1-h[i-1,2])*s*m[2,1]) + (NT[i-1,9]*(1-h[i-1,3])*s*m[3,1]);
        NT[i,8] = (NT[i-1, 7]*(1-h[i-1,1])*s*m[1,2]) + (NT[i-1,8]*(1-h[i-1,2])*s*m[2,2]) + (NT[i-1,9]*(1-h[i-1,3])*s*m[3,2]);
        NT[i,9] = (NT[i-1, 7]*(1-h[i-1,1])*s*m[1,3]) + (NT[i-1,8]*(1-h[i-1,2])*s*m[2,3]) + (NT[i-1,9]*(1-h[i-1,3])*s*m[3,3]);
        //NT[i, 9] = NT[i,9] + new_tags[3];
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
   
    
    
  //posterior predictive-- Abundance
    PostPredBr = lognormal_rng(log(PredAbundBr), sig[1:ABr]);
    PostPredMatInd = lognormal_rng(log(q[1] * PredAbundMatInd), sig[IndCVIdx]); //Indian 
    PostPredZero = normal_rng(q[1]*PredAbundNorm, sig[ZeroCVIdx]);
    PostPredMatPac = lognormal_rng(log(q[2]*PredAbundMatPac), sig[PacCVIdx]);//Pacific
    
    
    //posterior predictive -- Mark recoveries
    for (i in 1:(Nyeartag-1)){
      //PostPredRec[i,] = poisson_rng(NTR[i,]);
      PostPredRec[i,] = neg_binomial_2_rng(NTR[i,], theta);
   }
   //loglikelihoods for model selection--NOTE THIS DOESN'T WORK FOR THE NEW GROUPS
   //for(j in 1:ABr){
     //loglikeAbundBr[j] = lognormal_lpdf(AEstBr[j]|log(PredAbundBr[j]), sig[j]);
   //}
   //for(j in 1:AMat2){
     //loglikeAbundMatInd[j] = lognormal_lpdf(AEstMatInd[j]|log(q[1]*PredAbundMatInd[j]), sig[j + ABr]);
     //loglikeAbundMatPac[j] = lognormal_lpdf(AEstMatPac[j]| log(q[2]*PredAbundMatPac[j]), sig[j + AMat2 + ABr]);
   //}
   //for(t in 1:(Nyeartag-1)){
     //for(k in 1:(Nbasin*Nbasin)){
       //loglikeRec[t,k] = poisson_lpmf(Rec[t,k]|NTR[t,k]);
       //loglikeRec[t,k] = neg_binomial_2_lpmf(Rec[t,k]|NTR[t,k], theta);
     //}
     
   }
   
   //log_like[1:ABr] = loglikeAbundBr;
   //log_like[(ABr+1):(AMat2+ABr)] = loglikeAbundMatInd;
   //log_like[(ABr+AMat2+1):(ABr+AMat2+AMat2)] = loglikeAbundMatPac;
   //log_like[(ABr+AMat2+AMat2+1):(ABr+AMat2+AMat2 + ((Nyeartag-1)*Nbasin*Nbasin))] = to_vector(loglikeRec);
   
   //}
   
}







