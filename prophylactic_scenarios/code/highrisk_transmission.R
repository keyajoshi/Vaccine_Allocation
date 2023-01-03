##########################################
## This script generates the data for 
## the bottom panel of figure 5 which 
## considers heterogeneous risk structure 
## among populations with equal size 
## allowing individuals to be at higher risk
## of transmission
##########################################

###################################
## Import libraries
####################################
library(deSolve)
library(grid)
library(gridExtra)
library(viridis)
library(tidyverse)
options(scipen =999)

####################################
## Define Parameters
####################################
N1 <- 1e6 ## total population size 
N2 <- N1 ## set to be same 
phi1 <- 0 ## proportion immune in pop 1 
phi2 <- 0 ## proportion immune in pop 2
ve <- 0.95
tau <- 1-ve ## vaccine failure rate 
f1 <- 1 ## prop high-risk we want to vaccinate in pop 1 (and we give the rest to low risk)
f2 <- f1 ## prop high-risk we want to vaccinate in pop 2 (and we give the rest to low risk)
years <- 3

dt <- seq(0,365*years, 1) ## daily for 3 years


###################
## function
###################
SEIR_AN_HIGHLOW <- function(t, x, params){
  
  with(as.list(c(params,x)),{

    dS1H = -(beta1HH * S1H * I1H) - 
      (beta1HL * S1H * I1L)
    
    dV1H = 0
    
    dE1H = (beta1HH * S1H * I1H) + 
      (beta1HL * S1H * I1L) - 
      (sigma * E1H)
    
    dI1H = (sigma * E1H) - (nu * I1H) 
    
    dR1H = (nu * I1H) 
    
    dC1H = (sigma *E1H) 
    
    dS1L = -(beta1LL * S1L * I1L) - 
      (beta1LH * S1L * I1H)
    
    dV1L = 0
    
    dE1L = (beta1LL * S1L * I1L) + 
      (beta1LH * S1L * I1H) - 
      (sigma * E1L)
    
    dI1L = (sigma *E1L) - (nu * I1L)
    
    dR1L = (nu * I1L)
    
    dC1L = (sigma *E1L)
    
    dS2H = -(beta2HH * S2H * I2H) - 
      (beta2HL * S2H * I2L)
    
    dV2H = 0
    
    dE2H = (beta2HH * S2H * I2H) + 
      (beta1HL * S2H * I2L) - 
      (sigma * E2H)
    
    dI2H = (sigma *E2H) - (nu * I2H)
    
    dR2H = (nu * I2H)
    
    dC2H = (sigma *E2H) ## tracks cumulative infections 
    
    dS2L = -(beta2LL * S2L * I2L) - 
      (beta2LH * S2L * I2H)
    
    dV2L = 0
    
    dE2L = (beta2LL * S2L * I2L) + 
      (beta2LH * S2L * I2H) - (sigma * E2L)
    
    dI2L = (sigma *E2L) - (nu * I2L)
    
    dR2L = (nu * I2L)
    
    dC2L = (sigma *E2L) ## tracks cumulative infections 
    
    states = c(dS1H, dV1H, dE1H, dI1H, dR1H, dC1H, 
               dS1L, dV1L, dE1L, dI1L, dR1L, dC1L,
               dS2H, dV2H, dE2H, dI2H, dR2H, dC2H, 
               dS2L, dV2L, dE2L, dI2L, dR2L, dC2L)
    
    list(states)
  })
}

###########################
## Make Run
###########################

highrisk_transmission <- data.frame()

v_list <- seq(0, 1e6, 1e5)
pv_list <- seq(0,1, 0.1)
ph_list <- c(0.25, 0.5)
r_list <- c(2,4,8,16)
pct <- N1 * (1/1000)


for(v in v_list){
  for(pv1 in pv_list){
    for(ph1 in ph_list){
      for(r.val in r_list){
        ## increasing proportion of high risk 
        pct_hr <- (N1 *ph1) * (1/1000)
        pct_lr <- (N1 *(1-ph1)) * (1/1000)
        ph2 = ph1
        
        RHH = 8 *  ph1
        RHL = 4 * ph1
        RLH = 4 * (1-ph1)
        RLL = 2 * (1-ph1) 
  
        #### re-scale betas to make global R0 = 2
        
        scale.f <- (RHH + RLL + sqrt((RLL-RHH)^2 + 4 *(RHL*RLH)))/2
        ## to check 
        test <- (RHH + RLL + sqrt((-RHH-RLL)^2 - 4 *(RLL*RHH - RHL*RLH)))/2
        ## check again
        r0 = matrix(data = c(RHH, RLH, RHL, RLL), nrow = 2, ncol = 2, byrow = FALSE)
        eigen(r0)
        
        r.val/scale.f
        
        RHH = RHH * (r.val/scale.f)
        RHL = RHL * (r.val/scale.f)
        RLH = RLH * (r.val/scale.f)
        RLL = RLL * (r.val/scale.f)
        
        params <- c(beta1HH = RHH/(5*N1*ph1), ## per contact probability of transmission
                    beta1HL = RHL/(5*N1*ph1),
                    beta1LH = RLH/(5*N1*(1-ph1)),
                    beta1LL = RLL/(5*N1*(1-ph1)),
                    beta2HH = RHH/(5*N2 *ph2),
                    beta2HL = RHL/(5*N2*ph2),
                    beta2LH = RLH/(5*N2*(1-ph2)),
                    beta2LL = RLL/(5*N2*(1-ph2)),
                    sigma = 1/3, ## 1/latent period
                    nu = 1/5)
        
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
          S1L = (N1* (1-ph1)) - V1L -pct_lr
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
          S2H = (N2 * ph2) - V2H -pct_hr
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
          
          S1L = S1L, 
          V1L = V1L,
          E1L = 0,
          I1L = I1L, 
          R1L = 0,
          C1L = C1L, 
          
          S2H = S2H,
          V2H = V2H, 
          E2H = 0, 
          I2H = I2H,
          R2H = 0,
          C2H = C2H, 
          
          S2L = S2L,
          V2L = V2L,
          E2L = 0, 
          I2L = I2L,
          R2L = 0,
          C2L = C2L
        )
        
        output2 <- lsoda(y = inits,
                         times = dt,
                         func = SEIR_AN_HIGHLOW,
                         parms = params)
        
        out2.df <- data.frame(output2)[length(dt),]
        out2.df$v = v
        out2.df$pv = pv1
        out2.df$ph1 = ph1
        out2.df$ph2 = ph2
        out2.df$R = r.val
        highrisk_transmission = rbind(highrisk_transmission,out2.df)
        
        if(nrow(highrisk_transmission) %%100 == 0){
          print(nrow(highrisk_transmission))}
      }
    }
  }
}

## to check Ns 
highrisk_transmission <- highrisk_transmission%>% 
  rowwise() %>% 
  mutate(N1 = sum(c(S1H, V1H, E1H, I1H, R1H,
                    S1L, V1L, E1L, I1L, R1L)),
         N1H = sum(c(S1H, V1H, E1H, I1H, R1H)),
         N1L = sum(c(S1L, V1L, E1L, I1L, R1L)),
         N2 = sum(c(S2H, V2H, E2H, I2H, R2H,
                    S2L, V2L, E2L, I2L, R2L)),
         N2H = sum(c(S2H, V2H, E2H, I2H, R2H)),
         N2L = sum(c(S2L, V2L, E2L, I2L, R2L))) %>%
  mutate(vtot = V1H + V1L + V2H + V2L,
         vtau = v * (1-tau),
         vtotH = V1H + V2H,
         vtotL = V1L + V2L)

save(highrisk_transmission, file = "prophylactic_scenarios/data/highrisk_transmission.RData")
