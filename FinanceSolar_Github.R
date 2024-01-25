# Default Probability Estimation via Closed-form formula and Monte Carlo Simulation
# Code by: Hossein Jadidi (Email: hossein.jadidi@gmail.com)

# This code corresponds to the manuscript submitted to the Energy Economics Journal titled
# "Risk Mitigation in Project Finance for Utility-Scale Solar PV Projects."

# The code aims to estimate the default probability of project finance loans for utility-scale solar PV projects.
# Detailed information about the methodology is available in the manuscript.

# The output of the code (df) represents the cumulative default probability for each year of the loan life
# at the specified leverage ratio (L).

rm(list=ls())
library(stats4)
library(miscFuncs)
library(somebm)
library(ggplot2)
library(lognorm)
library(magrittr) 
library(dplyr)    
library(FinCal)
library(LSMRealOptions)
it <- 10000                 #iteration for Monte Carlo simulation
cap <- 10                   # Capacity of Solar PV plant in MW
LoanLife <- 20              #Loan life including grace period
ProjectLife <- 30           #Project life 
EnergyYield <- 15822        #Annual Energy yield of the plant in MW (The output of a simulation software for example PVsyst)
UC <- 9 * EnergyYield/100;  #Uncertainty of Energy yield
NL <- 9                     #number of Leverage ratios for analysis
WorkingCapital <- 0                  
Capex <- 0
Derv_cost <- 0
Degradation <- 0.5/100      #solar panel degradation in each year
Cost <- 600000              #Cost of Capital for Each MW
TotalCost <- cap*Cost       #Total cost of capital
Opex <- 0.5/100*TotalCost   #Yearly Operating Expenditure
rd <- 0.04                  #Cost of Debt
rf <- 0.03                  #risk free rate
re <- 0.1                   #Cost of equity
GracePeriod <- 1            #in Years
ConstPeriod <- 1            #in Years
Price <- 45                 #Market Spot Price EUR per MWh
CorTax <- 0.1               #Tax rate
#Computation
Production <- rnorm(it, EnergyYield,UC)       #Energy yield as a random variable
Opex <- rlnorm(it, log(Opex),0.03*log(Opex))  #Opex as a random variable
#Computation
  Fpp <- matrix(0, nrow=LoanLife+1 , ncol=NL )
  Leland <- matrix(0, nrow=LoanLife+1 , ncol=NL )
  for (ll in 1:NL){                     #In this loop I change the leverage ratio in each itereation
      Leverage <-0.7                    #Debt rate
      EquityRate <- 1-Leverage          #Equity rate
      LoanValue <- TotalCost*Leverage   #Total Debt value 
      wacc <- EquityRate * re + Leverage * rd * (1-CorTax)

#Debt Service computation
     
       TD <- LoanValue*(1+rd)^(ConstPeriod+GracePeriod)                 #Compound debt 
      DS <- (TD * rd)/(1-(1+rd)^-(LoanLife-(ConstPeriod+GracePeriod))) #Annuity formula
      OT <- matrix(0,nrow = LoanLife , ncol= 1)
      IT <- matrix(0,nrow = LoanLife , ncol= 1)
      PT <- matrix(0,nrow = LoanLife , ncol= 1)
      IT[2] <- LoanValue*(1+rd)^(ConstPeriod+GracePeriod)*rd
      PT[2] <- DS - IT[2] 
      OT[2]=LoanValue *(1+rd)^(ConstPeriod+GracePeriod)-PT[2]
      for (t in (2:LoanLife)) {
        IT[t] <- OT[t-1]*rd
        PT[t] <- DS-IT[t]
        OT[t]=OT[t-1]-PT[t]
      }
      
#Monte Carlo simulation to estimate the distribution parameters of DSCR for each year of the loan life 
      
      EY <- matrix(0,nrow = LoanLife , ncol= it)
      EBITDA <- matrix(0,nrow = LoanLife , ncol= it)
      Tax <- matrix(0,nrow = LoanLife , ncol= it)
      UFCF <- matrix(NA,nrow = LoanLife , ncol= it)
      DiffCash <- matrix(0,nrow = LoanLife , ncol= it)
      ExcessCash <- matrix(0,nrow = LoanLife , ncol= it)
      DSCR <- matrix(0,nrow = LoanLife , ncol= it)
      DSCR_dist_param <- matrix(0,nrow = LoanLife , ncol= 2)
      for (y in 2:LoanLife){
        EY[y,] <- Production*(1-Degradation)^y
        EBITDA[y,] <- EY[y,]*Price-Opex
        Tax[y,] <- CorTax*EBITDA[y,] 
        UFCF[y,] <- EBITDA[y,]-Tax[y,]-WorkingCapital-Capex+IT[y]*CorTax-Derv_cost
        DiffCash[y,] <- UFCF[y,]-DS 
        ExcessCash[y,] <- (DiffCash[y,]>0)*(DiffCash[y,]>0)*DiffCash[y,]
        UFCF[y,] <- UFCF[y,]+ExcessCash[y-1,]*(1+rf)  
      }
      DSCR <- UFCF/DS
      for (y in 2:LoanLife){
        estimation<-estimateParmsLognormFromSample(DSCR[y,],na.rm = FALSE)
        DSCR_dist_param[y,1] <- estimation[1]             #First parameter of the Lognormal distribution
        DSCR_dist_param[y,2] <- estimation[2]  }          #Second parameter of the Lognormal distribution
        output <- matrix(0,nrow = LoanLife , ncol= 7)
      for (i in 2:LoanLife){
        output[i,1] <- i
        output[i,2] <- DSCR_dist_param[i,1]
        output[i,3] <- DSCR_dist_param[i,2]
        output[i,4] <- log(abs(DSCR_dist_param[i,1]/DSCR_dist_param[2,1]))/i
        output[i,5] <- sqrt((log(1+(DSCR_dist_param[i,2]/(DSCR_dist_param[2,1]^2*exp(2*output[i,4]*i))))/i))
        output[i,7] <- DSCR_dist_param[2,1]^2*exp(2*output[i,4]*i)*(exp(output[i,5]^2*i)-1) #for checking
        output[i,6] <- DSCR_dist_param[2,1]*exp(i*output[i,4]) #for checking
      }
      mu <- mean(output[2:LoanLife,4])              #Estimated parameter (Drift) for Equation 6 of the paper 
      sigma <- mean(output[2:LoanLife,5])           #Estimated parameter (Noise) for Equation 6 of the paper

        
#Second part: Calculation of Cumulative Default Probability from the proposed closed form formula (Equation 6) 
  tt <- seq(0,LoanLife,by=1)
  D <- 1                          #Barrier level (I assume 1 as a material default in Project Finance)
  k0 <- log(exp(DSCR_dist_param[2,1])/D)
  t2 <- seq(2,LoanLife,by=1)
  kt <- log(exp(DSCR_dist_param[t2,1])/D)
  kt <- c(0,kt)
  d1 <- (-k0-(mu-0.5*sigma^2)*(tt))/(sigma*sqrt(tt))
  d1 <- d1[2:21]
  d2 <- (-k0+(mu-0.5*sigma^2)*(tt))/(sigma*sqrt(tt))
  d2 <- d2[2:21]
  d3 <- (-2*kt*(mu-0.5*sigma^2))/sigma^2
  Default_Leland <- pnorm(d1)+exp(d3)*pnorm(d2)
  Default_Leland <- c(0,Default_Leland)
  
#Check the Cumulative default probability by MCS method
  
  FirstYear <- matrix(NA, nrow=it , ncol=LoanLife+1 )
  Initial <- exp(DSCR_dist_param[2,1])
  countyear <- matrix(NA, nrow=it , ncol=LoanLife+1 )
  fpp <- matrix(0, nrow=LoanLife+1 , ncol=1 )
  obj <- GBM_simulate(it,LoanLife,mu,sigma,Initial,1)
  obj <- t(obj)
  count_default <- obj[] < D
  
  for (kk in 1:it) {
    aa <- cumsum(count_default[kk,])
    FirstYear[kk,] <- (count_default[kk,]==1)*(count_default[kk,]==1) #remove more than one default
    while (sum(FirstYear[kk,]) > 1){                            #The default can be happened just 
      bb <- FirstYear[kk,] == 1                                 #one time (First passage probability)
      countyear[kk,] <- cumsum(bb)
      FirstYear[kk,] <- (countyear[kk,]==1)*(countyear[kk,]==1) #hold just the first default
    }
  }
  
  for (n in 1:LoanLife+1) {
    fpp[n] <-sum(FirstYear[,n])/it                              #first passage probability
  }
  Fpp[,ll]<-t(cumsum(fpp))                      #Output of MCS 
  Leland[,ll]<-t(Default_Leland)                #Output of Closed-form
}

df_Fpp <- data.frame(tt,Fpp*100,Leland*100)

colnames(df_Fpp) <- c("tt",sprintf("Monte Carlo_L[%d]",seq(50,90,by=5)),
                      sprintf("Closed-form_L[%d]",seq(50,90,by=5))
                      )
LvGraph <- "L[70]"       #change[Leverage number] to show the final graph based on the input number
df <- df_Fpp %>%         # For instance change 70 to 85 to see the results for 85% leverage ratio
  select(tt,ends_with(LvGraph)) %>%
  tidyr::gather(key = "Method", value = "value", -tt)

a <- ggplot(df, aes(x=tt, y=value ,col=Method))   
a+
  geom_line(aes(linetype=Method, color=Method),size=1.3)+
  geom_point(aes(color=Method,shape=Method),size=4)+
  labs(x = "Loan year",y = "Cumulative default probability (%)")+
  theme_bw()+
  theme(legend.position = c(0.8, 0.2))+
  theme(axis.text = element_text(size = 16))+
  theme(text=element_text(family="serif", size=16))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(text=element_text(family="serif", size=16)) 


