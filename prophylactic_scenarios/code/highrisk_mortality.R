##########################################
## This script generates the data for 
## the bottom panel of figure 5 which 
## considers heterogeneous risk structure 
## among populations with equal size 
## allowing individuals to be at higher risk
## of death
##########################################

###################################
## Import libraries
####################################
library(deSolve)
library(grid)
library(gridExtra)
library(viridis)
library(tidyverse)
options(scipen = 999)
####################################
## Define Parameters
####################################
N1 <- 1e6 ## total population size 
N2 <- 1e6 ## set to be same 
phi1 <- 0 ## proportion immune in pop 1 
phi2 <- 0 ## proportion immune in pop 2
ve <- 0.95 ## vaccine efficacy
tau <- 1-ve ## vaccine failure rate 
f1 <- 1 ## prop high-risk we want to vaccinate in pop 1 (and we give the rest to low risk)
f2 <- f1 ## prop high-risk we want to vaccinate in pop 2 (and we give the rest to low risk)
years <- 3 
dt <- seq(0,365*years, 1) ## daily for 3 years

SEIR_AN_HIGHLOW_DEATHS <- function(t, x, params){
  
  with(as.list(c(params,x)),{
    
    N1H = S1H + E1H + I1H + R1H + V1H
    N2H = S2H + E2H + I2H + R2H + V2H
    
    N1L = S1L + E1L + I1L + R1L + V1L
    N2L = S2L + E2L + I2L + R2L + V2L
    
    
    dS1H = -(beta1HH * S1H * I1H) - 
      (beta1HL * S1H * I1L)
    
    dV1H = 0
    
    dE1H = (beta1HH * S1H * I1H) + 
      (beta1HL * S1H * I1L) - 
      (sigma * E1H)
    
    dI1H = (sigma * E1H) - 
      (nu * I1H) 
    
    dR1H = sh *(nu * I1H)
    
    dC1H = (sigma *E1H) 
    
    dD1H = (1-sh) *(nu * I1H) 
    
    dS1L = -(beta1LL * S1L * I1L) - 
      (beta1LH * S1L * I1H)
    
    dV1L = 0
    
    dE1L = (beta1LL * S1L * I1L) + 
      (beta1LH * S1L * I1H) - 
      (sigma * E1L)
    
    dI1L = (sigma *E1L) - 
      (nu * I1L)
    
    dR1L = sl * (nu * I1L)
    
    dC1L = (sigma *E1L) 
    
    dD1L = (1-sl) * (nu * I1L)
    
    dS2H = -(beta2HH * S2H * I2H) - 
      (beta2HL * S2H * I2L)
    
    dV2H = 0
    
    dE2H = (beta2HH * S2H * I2H) + 
      (beta2HL * S2H * I2L) - 
      (sigma * E2H)
    
    dI2H = (sigma *E2H) - 
      (nu * I2H)
    
    dR2H = sh * (nu * I2H)
    
    dC2H = (sigma *E2H)  
    
    dD2H = (1-sh) * (nu * I2H)
    
    dS2L =  -(beta2LL * S2L * I2L) - 
      (beta2LH * S2L * I2H)
    
    dV2L = 0
    
    dE2L = (beta2LL * S2L * I2L) + 
      (beta2LH * S2L * I2H) - 
      (sigma * E2L)
    
    dI2L = (sigma *E2L) - 
      (nu * I2L)
    
    dR2L = sl * (nu * I2L)
    
    dC2L = (sigma *E2L) ## tracks cumulative infections
    
    dD2L = (1-sl) * (nu * I2L)
    
    states = c(dS1H, dV1H, dE1H, dI1H, dR1H, dC1H, dD1H,
               dS1L, dV1L, dE1L, dI1L, dR1L, dC1L, dD1L,
               dS2H, dV2H, dE2H, dI2H, dR2H, dC2H, dD2H,
               dS2L, dV2L, dE2L, dI2L, dR2L, dC2L, dD2L)
    
    list(states)
  })
}

###########################
## Make Run
###########################

highrisk_mortality <- data.frame()

## proportion vaccinated pv_list
## vaccine seq(0, 1000, by = 100)
## proportion high v low risk seq(0,1,0.1) ph_list
## underlying immunity phi_list
v_list <- seq(0, 1e6, 1e5)
pv_list <- seq(0,1, 0.1)
ph_list <- c(0.25, 0.5)
r_list <- c(2,4,8,16)
pct <- N1 * (1/1000)

for(v in v_list){
  for(pv1 in pv_list){
    for(ph1 in ph_list){
      for(RHH1 in r_list){
        
        pct_hr <- (N1 *ph1) * (1/1000)
        pct_lr <- (N1 *(1-ph1)) * (1/1000)
        
        ph2 = ph1
        
        ## This will vary for population 1 
        
        RHL1 = RHH1
        RLH1 = RHH1
        RLL1 = RHH1 
        
        ## Fixed at baseline for population 2 
        RHH2 = RHH1  
        RHL2 = RHH1
        RLH2 = RHH1
        RLL2 = RHH1
        
        ## use params from kissler et al. 
        params <- c(beta1HH = RHH1/(5*N1*ph1), ## per contact probability of transmission
                    beta1HL = RHL1/(5*N1*ph1),
                    beta1LH = RLH1/(5*N1*(1-ph1)),
                    beta1LL = RLL1/(5*N1*(1-ph1)),
                    beta2HH = RHH2/(5*N2 *ph2),
                    beta2HL = RHL2/(5*N2*ph2),
                    beta2LH = RLH2/(5*N2*(1-ph2)),
                    beta2LL = RLL2/(5*N2*(1-ph2)),
                    sigma = 1/3, ## 1/latent period
                    nu = 1/5, ## 1/infectious period)
                    mu = 1/(80*365),
                    sh = 0.995,
                    sl = 0.999) 
        
        ## want to continuously give the high risk 
        ## the number of vaccines available 
        
        ## if there are more high risk than the 
        ## number of vaccines then all vaccines 
        ## go to the high risk 
        
        ## if are are less high risk people 
        ## than vaccines, then give as many
        ## vaccines as you can to high risk
        ## then left overs go to low risk 
        
        if((N1*ph1*f1) > (v*pv1)){
          V1H = v*pv1*(1-tau)
          V1L = 0
        } else{
          V1H = N1*f1*ph1* (1-tau) 
          V1L = v*pv1*(1-tau) - N1*f1*ph1*(1-tau)
        }
        # 
        # ## Initital conditions for population 1 
        if(ph1 == 0){
          S1H = 0
          I1H = 0
          C1H = 0
        } else {
          S1H = (N1*ph1) - V1H - pct_hr
          I1H = pct_hr
          C1H = pct_hr
        }
        if(ph1 == 1){
          S1L = 0
          I1L = 0
          C1L = 0
        } else {
          S1L = (N1* (1-ph1)) - V1L - pct_lr 
          I1L = pct_lr
          C1L = pct_lr
          
        }
        
        if((N2*ph2*f2) > (v*(1-pv1))){
          V2H = v*(1-pv1)*(1-tau)
          V2L = 0
        } else{
          V2H = N2*f2*ph2*(1-tau) 
          V2L = v*(1-pv1)*(1-tau) - N2*f2*ph2*(1-tau)
        }
        
        if(ph2 == 0){
          S2H = 0
          I2H = 0
          C2H = 0
          
        } else {
          S2H = (N2 * ph2) - V2H - pct_hr
          I2H = pct_hr
          C2H = pct_hr
        }
        if(ph2 == 1){
          S2L = 0
          I2L = 0
          C2L = 0
        } else {
          S2L = (N2 * (1-ph2)) - V2L - pct_lr
          I2L = pct_lr
          C2L = pct_lr
          
        }
        
        inits<-c(
          S1H = S1H,
          V1H = V1H, 
          E1H = 0, 
          I1H = I1H,
          R1H = 0,
          C1H = C1H, 
          D1H = 0,
          
          S1L = S1L, 
          V1L = V1L,
          E1L = 0,
          I1L = I1L, 
          R1L = 0,
          C1L = C1L, 
          D1L = 0,
          
          S2H = S2H,
          V2H = V2H, 
          E2H = 0, 
          I2H = I2H,
          R2H = 0,
          C2H = C2H, 
          D2H = 0,
          
          S2L = S2L,
          V2L = V2L,
          E2L = 0, 
          I2L = I2L,
          R2L = 0,
          C2L = C2L,
          D2L = 0
        )
        
        output2 <- lsoda(y = inits,
                         times = dt,
                         func = SEIR_AN_HIGHLOW_DEATHS,
                         parms = params)
        
        out2.df <- data.frame(output2)[length(dt),]
        out2.df$v = v
        out2.df$pv = pv1
        out2.df$ph1 = ph1
        out2.df$ph2 = ph2
        out2.df$phi1 = phi1
        out2.df$RHH1 = RHH1
        highrisk_mortality = rbind(highrisk_mortality,out2.df)
        
        if(nrow(highrisk_mortality) %%100 == 0){
          print(nrow(highrisk_mortality))}
        
      } 
    }
  }
}

### check what is going on with the Ns 
highrisk_mortality <- highrisk_mortality %>% 
  rowwise() %>% 
  mutate(N1 = sum(c(S1H, V1H, E1H, I1H, R1H, D1H,
                    S1L, V1L, E1L, I1L, R1L, D1L)),
         N1H = sum(c(S1H, V1H, E1H, I1H, R1H, D1H)),
         N1L = sum(c(S1L, V1L, E1L, I1L, R1L, D1L)),
         N2 = sum(c(S2H, V2H, E2H, I2H, R2H,D2H, 
                    S2L, V2L, E2L, I2L, R2L, D2L)),
         N2H = sum(c(S2H, V2H, E2H, I2H, R2H, D2H)),
         N2L = sum(c(S2L, V2L, E2L, I2L, R2L, D2L))) %>%
  mutate(vtot = V1H + V1L + V2H + V2L,
         vtau = v * (1-tau),
         vtotH = V1H + V2H,
         vtotL = V1L + V2L)

save(highrisk_mortality, file =
        "prophylactic_scenarios/data/highrisk_mortality.Rdata")
