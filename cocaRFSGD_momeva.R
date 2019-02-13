##################################################
## Publication:    "Quantifying regional fresh Submarine Groundwater Discharge with the lumped modeling approach CoCa-RFSGD" 
## Script purpose: Calculate mean Fresh Submarine Groundwater Discharge (FSGD) with global catchment model parameter mean values (momeva)
## Date:           2019/01/03
## License:        GNU GPLv3
## Author:         Mithra-Christin Hajati
## Co-Authors:     Edwin Sutanudjaja and Nils Moosdorf
##################################################

setwd("/Users/<username>/Desktop/cocaRFSGD")

#----> LOAD Input Data <----#
Co_ca_LRM   <- read.csv(paste0(getwd(),"/model_input/Coca.csv"))                                         # Model Parameters lumped for each coastal catchment
ET          <- read.csv(paste0(getwd(),"/model_input/ET.csv"));             ET$Date <- as.Date(ET$Date)  # Evapotranspiration lumped for each coastal catchment
P           <- read.csv(paste0(getwd(),"/model_input/P.csv"));              P$Date  <- as.Date(P$Date)   # Precipitation lumped for each coastal catchment
In.cond     <- read.csv(paste0(getwd(),"/model_input/In.cond_momeva.csv")); In.cond$cat <- Co_ca_LRM$cat # Model Initial Conditions for top soil, sub soil and aquifer

#----> Model Settings <----#
date_start  <- as.Date("2001/01/01",format='%Y/%m/%d')                      # Choose Start date of the model [days]
date_end    <- as.Date("2014/12/30",format='%Y/%m/%d')                      # Choose End date of the model [days]
timeStep    <- 1                                                            # Time step [days]
t_0         <- 1                                                            # Start date [days]
days_LRM    <- seq.Date(date_start,date_end,timeStep)                       # Time period [days]
t_max       <- length(days_LRM)                                             # Amount of modeled days [days]
ET_LRM      <- ET[ET$Date>=days_LRM[1] & ET$Date<=days_LRM[t_max],]         # Precipitation of the model period
P_LRM       <- P[P$Date>=days_LRM[1] & P$Date<=days_LRM[t_max],]            # Evapotranspiration of the model period
cocas       <- Co_ca_LRM$cat                                                # Set which coastal catchments to be modeled

#----> Create Model Result Files <----#
S3     <- seq(0,0,length.out=t_max)                                          # Aquifer water volume
S2     <- seq(0,0,length.out=t_max)                                          # Sub soil water volume
S1     <- seq(0,0,length.out=t_max)                                          # Top soil water volume
SGD    <- P_LRM                                                              # FSGD
OF     <- data.frame(Date=days_LRM); OF[,2:134]      <- (ET_LRM[,2:134] * 0) # Overflow 
ET_act <- data.frame(Date=days_LRM); ET_act[,2:134]  <- (ET_LRM[,2:134] * 0) # actual Evaporation
R1     <- data.frame(Date=days_LRM); R1[,2:134]      <- (ET_LRM[,2:134] * 0) # Recharge top to sub soil
cr1    <- data.frame(Date=days_LRM); cr1[,2:134]     <- (ET_LRM[,2:134] * 0) # Capillary rise
R2     <- data.frame(Date=days_LRM); R2[,2:134]      <- (ET_LRM[,2:134] * 0) # Recharge sub soil to aquifer
momeva <- data.frame(Co_Ca=Co_ca_LRM$cat, Lith=Co_ca_LRM$litho_median,       # Summary of all results
                     SGDsd=NA, SGDmean=NA, SGDtot=NA, OFtot=NA, CRtot=NA, 
                     R1tot=NA, R2tot=NA, k1=NA, k2=NA, k3=NA, K_r3=NA, 
                     S1max=NA, S2max=NA, D3=NA,  P=NA, ET_act=NA, ET_pot=NA)

k<-0
for(j in Co_ca_LRM$cat){ 
  k<-k+1
  
  # Recieving Soil and Aquifer Parameters 
  wc1_sat   <- Co_ca_LRM$S1_WCsat_average[Co_ca_LRM$cat==j]/10000                                   # TOP SOIL saturated water content [-]
  wc1_res   <- Co_ca_LRM$S1_WCres_average[Co_ca_LRM$cat==j]/10000                                   # TOP SOIL residual water content [-]
  S1_max    <- Co_ca_LRM$D1_average[Co_ca_LRM$cat==j]                                               # TOP SOIL max water volume [mm]
  V1tot     <- S1_max/wc1_sat                                                                       # TOP SOIL total Volume [mm]
  S1_min    <- wc1_res*V1tot                                                                        # TOP SOIL min water volume [mm]
  n1        <- Co_ca_LRM$S1_n_average[Co_ca_LRM$cat==j]/10000                                       # TOP SOIL van Genuchten parameter n [-]
  m1        <- 1 - 1/n1                                                                             # TOP SOIL van Genuchten parameter m [-]
  alpha1    <- 1/10 * (Co_ca_LRM$S1_alpha_average[Co_ca_LRM$cat==j]/10000)                          # TOP SOIL van Genuchten parameter aplha [1/mm]
  k1        <- 10*Co_ca_LRM$S1_k_average[Co_ca_LRM$cat==j]/10000                                    # TOP SOIL saturated hydraulic conductivity [mm/day]
  wc2_sat   <- Co_ca_LRM$S2_WCsat_average[Co_ca_LRM$cat==j]/10000                                   # SUB SOIL saturated water content [-]
  wc2_res   <- Co_ca_LRM$S2_WCres_average[Co_ca_LRM$cat==j]/10000                                   # SUB SOIL residual water content [-]
  V2tot     <- 700                                                                                  # SUB SOIL total Volume [mm]
  S2_max    <- wc2_sat*V2tot                                                                        # SUB SOIL max water volume [mm]
  S2_min    <- wc2_res*V2tot                                                                        # SUB SOIL min water volume [mm]
  k2        <- 10*Co_ca_LRM$S2_k_average[Co_ca_LRM$cat==j]/10000                                    # SUB SOIL saturated hydraulic conductivity [mm/day]
  n2        <- Co_ca_LRM$S2_n_average[Co_ca_LRM$cat==j]/10000                                       # SUB SOIL van Genuchten parameter n [-]
  m2        <- 1 - 1/n2                                                                             # SUB SOIL van Genuchten parameter m [-]
  alpha2    <- 1/10 * (Co_ca_LRM$S2_alpha_average[Co_ca_LRM$cat==j]/10000)                          # SUB SOIL van Genuchten parameter alpha [1/mm]
  D3        <- max(1000,(1000*Co_ca_LRM$Dtot_average[Co_ca_LRM$cat==j])-(V1tot + V2tot))            # AQUIFER Depth [mm]
  epsilon3  <- Co_ca_LRM$S3_porosity_average[Co_ca_LRM$cat==j]                                      # AQUIFER Porosity [-]
  k3        <- (10^(Co_ca_LRM$S3_log_perm_avg[Co_ca_LRM$cat==j])*1000*9.81/0.001)*(1000*(60*60*24)) # AQUIFER saturated hydraulic conductivity [mm/day]
  L         <- 1000 * Co_ca_LRM$L[Co_ca_LRM$cat==j]/2                                               # AQUIFER Average flow width [mm]
  K_r3      <- (pi^2 * k3 * D3) / (epsilon3 * L^2)                                                  # AQUIFER Recession Parameter [1/d] (Kraaijenhoff van de Leur, 1958)
  
  # INCREASE spinning up parameters
  x3    <- 30
  diff2 <- seq(0,0,length.out=t_max)
  i    <- 3 
  
  ### --- START WATER BALANCE MODEL --- ###
  repeat{
    i<-i+1 # Increases spinning up process
    # ASSIGN Initial Conditions
    S1[1]  <- In.cond$S1_0[grep(which(In.cond$cat==j) , rownames(In.cond))][1]                   # Set top soil water volume of day 1 same as last day
    S2[1]  <- In.cond$S2_0[grep(which(In.cond$cat==j) , rownames(In.cond))][1]                   # Set sub soil water volume of day 1 same as last day
    S3[1]  <- In.cond$S3_0[grep(which(In.cond$cat==j) , rownames(In.cond))][1] + x3 * diff2[i-1] # Set aquifer water volume of day 1 same as last day
    
    for(l in seq(t_0, t_max, by=timeStep)){
      # Top Soil Input P - ET
      S1[l+1]                     <- S1[l] + P_LRM[l,paste0("CoCa.",j)]                                     
      ET_act[l,paste0("CoCa.",j)] <- min(ET_LRM[l,paste0("CoCa.",j)],S1[l+1]-S1_min)
      S1[l+1]                     <- S1[l+1] - ET_act[l,paste0("CoCa.",j)]                                     
      # Overflow?! If yes, subtract from top soil water
      OF[l,paste0("CoCa.",j)]     <- max(0,S1[l+1] - S1_max)                                                   
      S1[l+1]                     <- S1[l+1] - OF[l,paste0("CoCa.",j)]                                         
      # Calculate Sub Soil Recharge R1 and add to Water Balance
      if((S1[l+1]/V1tot)<=wc1_res){
        K_r1 <- 0
      }else if((S1[l+1]/V1tot)>=wc1_sat){
        h1=0
        K_r1                      <- k1 * ((1 - ((alpha1*h1)^(n1-1)) * (1+(((alpha1*h1))^n1))^(-m1))^2)/((1+(alpha1*h1)^n1)^(m1/2))
      }else{
        h1                      <- (((((S1[l+1]/V1tot)-wc1_res)/(wc1_sat-wc1_res))^(-1/m1) -1)/alpha1^n1)^(1/n1)
        K_r1                    <- k1 * ((1 - ((alpha1*h1)^(n1-1)) * (1+((alpha1*h1)^n1))^(-m1))^2)/((1+(alpha1*h1)^n1)^(m1/2))
      }
      R1[l,paste0("CoCa.",j)]     <- K_r1                                                            
      R1[l,paste0("CoCa.",j)]     <- min(R1[l,paste0("CoCa.",j)],S1[l+1]-S1_min)
      S1[l+1]                     <- S1[l+1] - R1[l,paste0("CoCa.",j)]                                       
      S2[l+1]                     <- S2[l] + R1[l,paste0("CoCa.",j)]                                          
      # Overflow?! If yes, subtract from top soil water
      OF2_new                     <- max(0,S2[l+1] - S2_max)
      S2[l+1]                     <- S2[l+1] - OF2_new
      S1[l+1]                     <- S1[l+1] + OF2_new
      OF1_new                     <- max(0,S1[l+1] - S1_max)
      S1[l+1]                     <- S1[l+1] - OF1_new
      OF[l,paste0("CoCa.",j)]     <- OF[l,paste0("CoCa.",j)]  + OF1_new                                      
      # Calculate Aquifer Recharge R2 and add to Water Balance
      if((S2[l+1]/V2tot)>wc2_res & (S2[l+1]/V2tot)<wc2_sat){
        h2                        <- (((((S2[l+1]/V2tot)-wc2_res)/(wc2_sat-wc2_res))^(-1/m2) -1)/alpha2^n2)^(1/n2)
      }else{h2=0}
      K_r2                        <- k2 * ((1 - ((alpha2*h2)^(n2-1)) * (1+((alpha2*h2)^n2))^(-m2))^2)/((1+(alpha2*h2)^n2)^(m2/2))
      R2[l,paste0("CoCa.",j)]     <- K_r2                                                                     
      R2[l,paste0("CoCa.",j)]     <- min(R2[l,paste0("CoCa.",j)],S2[l+1]-S2_min)
      S2[l+1]                     <- S2[l+1] - R2[l,paste0("CoCa.",j)]                                        
      S3[l+1]                     <- S3[l] + R2[l,paste0("CoCa.",j)]                                         
      # Calculate Capillary Rice CR and add to Water Balance
      if((S2[l+1]/V2tot)>wc2_res & (S2[l+1]/V2tot)<wc2_sat){
        h2                        <- (((((S2[l+1]/V2tot)-wc2_res)/(wc2_sat-wc2_res))^(-1/m2) -1)/alpha2^n2)^(1/n2)
      }else{h2=0}
      K_r2                        <- k2 * ((1 - ((alpha2*h2)^(n2-1)) * (1+((alpha2*h2)^n2))^(-m2))^2)/((1+(alpha2*h2)^n2)^(m2/2))
      cr1[l,paste0("CoCa.",j)]    <- K_r2 * (1-(S1[l+1]/S1_max))                                               
      cr1[l,paste0("CoCa.",j)]    <- min(cr1[l,paste0("CoCa.",j)],S2[l+1]-S2_min)
      S1[l+1]                     <- S1[l+1] + cr1[l,paste0("CoCa.",j)]                                     
      S2[l+1]                     <- S2[l+1] - cr1[l,paste0("CoCa.",j)]                                      
      # Overflow?! If yes, subtract from top soil water 
      OF1_new                     <- max(0,S1[l+1] - S1_max)
      OF[l,paste0("CoCa.",j)]     <- OF[l,paste0("CoCa.",j)] + OF1_new                                        
      S1[l+1]                     <- S1[l+1] - OF1_new                                                      
      # Calculate FSGD and add to Water Balance
      SGD[l,paste0("CoCa.",j)]    <- K_r3*timeStep*S3[l+1]/2 # = SGD
      SGD[l,paste0("CoCa.",j)]    <- min(S3[l+1],SGD[l,paste0("CoCa.",j)])                                    
      S3[l+1]                     <- S3[l+1] - SGD[l,paste0("CoCa.",j)]                                     
    }
    
    # Temporary plots to check the Model results
    par(mfrow=c(3,1))
    plot(seq(t_0, t_max, by=timeStep),S1[-1], typ='l', ann=FALSE)
    plot(seq(t_0, t_max, by=timeStep),S2[-1], typ='l', ann=FALSE)
    plot(seq(t_0, t_max, by=timeStep),S3[-1], typ='l', ann=FALSE)
    
    # SET New Initial Conditions
    In.cond$S1_0[grep(which(In.cond$cat==j) , rownames(In.cond))][1] <- S1[length(S1)]
    In.cond$S2_0[grep(which(In.cond$cat==j) , rownames(In.cond))][1] <- S2[length(S2)]
    In.cond$S3_0[grep(which(In.cond$cat==j) , rownames(In.cond))][1] <- S3[length(S3)]
    
    # INCREASE spinning up by choosing parameters (1.5,4.5)
    diff2[i] <- S3[length(S3)] - S3[1]
    if(abs(diff2[i])-abs(diff2[i-2])>0){ 
      x3<-max(0,x3/1.5) 
    }else if(abs(diff2[i])<abs(diff2[i-2]*0.35)){ 
      x3<-x3 
    }else{ 
      x3<-x3*4
    }
    
    # Initiate END of LOOP if roughly the difference of SGD@starttime = SGD@endtime = Aquifer reached equilibrium
    if(abs(diff2[i]) < 1 | abs(diff2[i]) == abs(diff2[i-2])){
    
    # SAVE results in momeva  
    momeva[k,3:19] <- c(sd(unlist(SGD[paste0("CoCa.",j)])),
                        mean(unlist(SGD[paste0("CoCa.",j)])),
                        sum(SGD[paste0("CoCa.",j)])/14,
                        sum(OF[paste0("CoCa.",j)])/14,
                        sum(cr1[paste0("CoCa.",j)])/14,
                        sum(R1[paste0("CoCa.",j)])/14,
                        sum(R2[paste0("CoCa.",j)])/14,
                        k1,k2,k3,K_r3,S1_max,S2_max,D3,sum(P_LRM[,paste0("CoCa.",j)])/14,sum(ET_act[,paste0("CoCa.",j)])/14,sum(ET_LRM[,paste0("CoCa.",j)])/14)

    # Set Up Water Balance 
    S_end   <- S1[length(S3)] + S2[length(S3)] + S3[length(S3)]
    S_beg   <- S1[1] + S2[1] + S3[1]
    input   <- sum(P_LRM[,paste0("CoCa.",j)])
    output  <- sum(ET_act[,paste0("CoCa.",j)]) + sum(SGD[,paste0("CoCa.",j)]) + sum(OF[,paste0("CoCa.",j)])
    deltaS  <- S_end - S_beg
    iszero  <- deltaS - input + output
    
    # DO NOT SAVE RESULTS IF WATER BALANCE IS SMALLER 1
    if(-1 > iszero | iszero > 1 | abs(diff2[i]) == abs(diff2[i-2])){
        print(paste0("ERROR: Water Balance in ", colnames(SGD[paste0("CoCa.",j)]), " is ", format(round(iszero, 15), nsmall = 2)))
        momeva[k,3:21] <- c(NA,NA,NA,NA,NA,NA,NA,k1,k2,k3,K_r3,S1_max,D3,sum(P_LRM[,paste0("CoCa.",j)]),sum(ET_act[,paste0("CoCa.",j)]),sum(ET_LRM[,paste0("CoCa.",j)]))
    }
    
    # PRINT Coastal Catchment SGD mean and WB    
    cat(paste0(colnames(SGD[paste0("CoCa.",j)]),": ",round(diff2[i]), "mm difference in aquifer storage (spinupfactor = ", round(x3,digits=2), "), WB = ", iszero, ", SGD_mean = ", round(mean(unlist(SGD[paste0("CoCa.",j)])),digits=2),"\n"))
    
    # stop loop
    break
    }
  }
} 

# Save results
write.csv(In.cond , file=paste0(getwd(),"/model_input/In.cond_momeva.csv"), row.names=FALSE)
write.csv(momeva , file=paste0(getwd(),"/model_output/",format(Sys.Date(), '%Y%m%d'),"_momeva.csv"), row.names=FALSE)
save.image(paste0(getwd(),"/model_output/",format(Sys.Date(), '%Y%m%d'),"_momeva.RData"))

