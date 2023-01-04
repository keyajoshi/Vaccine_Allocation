##########################################
## This script generates figures S14
## This scenario considers rollout of vaccines
## in a population with heterogeneous risk structure
## where we have individuals at higher risk of death.
## We allow timing of rollout to vary between the two 
## populations 
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
N1 <- 1e6 ## population 1 size 
N2 <- 1e6 ## population 2 size 
ve <- 0.95 ## vaccine efficacy
tau <- 1-ve ## vaccine failure rate 
phi1 <- 0 ## proportion immune in pop 1 
phi2 <- 0 ## proportion immune in pop 2
f1 <- 1 ## prop high-risk we want to vaccinate in pop 1 (and we give the rest to low risk)
f2 <- f1 ## prop high-risk we want to vaccinate in pop 2 (and we give the rest to low risk)
years <- 3 ## duration of simulation 

dt <- seq(0,365*years, 1) ## daily for 3 years

################################
## Function for high risk mortality rollout vary time
################################

SEIR_AN_HIGHLOW_DEATHS_ro_vt <- function(t, x, params){
  
  with(as.list(c(params,x)),{
    
    N1H = S1H + E1H + I1H + R1H + V1H + D1H 
    N2H = S2H + E2H + I2H + R2H + V2H + D2H
    
    N1L = S1L + E1L + I1L + R1L + V1L + D1L
    N2L = S2L + E2L + I2L + R2L + V2L + D2L
    
    dS1H = -(beta1HH * S1H * I1H) - 
      (beta1HL * S1H * I1L)
    
    dE1H = (beta1HH * S1H * I1H) + 
      (beta1HL * S1H * I1L) - 
      (sigma * E1H)
    
    dI1H = (sigma * E1H) - 
      (nu * I1H)
    
    dR1H = sh* (nu * I1H)
    
    dC1H = (sigma *E1H) 
    
    dV1H = 0
    
    dVC1H = 0  
    
    dD1H = (1-sh) * (nu * I1H)
    
    
    dS1L =  -(beta1LL * S1L * I1L) - 
      (beta1LH * S1L * I1H)
    
    dE1L = (beta1LL * S1L * I1L) + 
      (beta1LH * S1L * I1H) - 
      (sigma * E1L)
    
    dI1L = (sigma *E1L) - 
      (nu * I1L) 
    
    dR1L = sl * (nu * I1L) 
    
    dC1L = (sigma *E1L) 
    
    dV1L = 0
    
    dVC1L = 0
    
    dD1L = (1-sl) * (nu*I1L)
    
    
    dS2H = -(beta2HH * S2H * I2H) - 
      (beta2HL * S2H * I2L)
    
    dE2H = (beta2HH * S2H * I2H) + 
      (beta2HL * S2H * I2L) - 
      (sigma * E2H)
    
    dI2H = (sigma *E2H) - 
      (nu * I2H)
    
    dR2H = sh * (nu * I2H)
    
    dC2H = (sigma *E2H) 
    
    dV2H = 0
    
    dVC2H = 0 
    
    dD2H = (1-sh) * (nu * I2H)
    
    
    dS2L = -(beta2LL * S2L * I2L) - 
      (beta2LH * S2L * I2H)
    
    dE2L = (beta2LL * S2L * I2L) + 
      (beta2LH * S2L * I2H) - 
      (sigma * E2L)
    
    dI2L = (sigma *E2L) - 
      (nu * I2L)
    
    dR2L = sl * (nu * I2L)
    
    dC2L = (sigma *E2L)
    
    dV2L = 0
    
    dVC2L = 0
    
    dD2L = (1-sl) * (nu * I2L)
    
    
    
    if(t > time.rollout) {
      
      ## now we are in roll-out 
      ## Option 1: we have more vaccine than number at high risk 
      ## Loop 1: Do we have less high risk individuals then vaccine? 
      ## Loop 2: Once all high risk get vaccine, give the left over vaccines
      ## to the low risk 
      ## Loop 3: stop vaccinating low risk when we run out of (effective) vaccines 
      
      if(N1*ph1*f1 <= v*pv1){
        
        dS1H  = -(beta1HH * S1H * I1H) -
          (beta1HL * S1H * I1L) - 
          (f * S1H)
        
        dR1H = sh * (nu * I1H) - 
          (f* R1H) 
        
        dVC1H = (f * S1H) + (f* R1H)
        
        if(VC1H > N1*ph1*f1*(1-tau)){
          
          ## stop vaccinating people in the high risk category 
          dS1H  = -(beta1HH * S1H * I1H) -
            (beta1HL * S1H * I1L)
          
          dR1H = sh * (nu * I1H)
          
          dVC1H = 0
          
          ## vaccinate people in low risk category 
          dS1L = -(beta1LL * S1L * I1L) - 
            (beta1LH * S1L * I1H) - 
            (f*S1L)
          
          dR1L = sl * (nu * I1L) - 
            (f*R1L)
          
          dVC1L = (f*S1L) + (f*R1L)
          
          if((VC1H + VC1L) > (v * pv1)*(1-tau)){
            
            ## stop vaccinating the low risk 
            
            dS1L = -(beta1LL * S1L * I1L) -
              (beta1LH * S1L * I1H)
            
            dR1L = sl * (nu * I1L)

            dVC1L = 0
          }
        }
      }
      
      ## Option 2: we have less vaccine than number at high risk 
      ## Loop 1: check we have more high risk than vaccine 
      ## Loop 2: did we use up all of our vaccines 
      
      if(v*pv1 < (N1* ph1 * f1)){
        
        dS1H  = -(beta1HH * S1H * I1H) -
          (beta1HL * S1H * I1L) - 
          (f * S1H)
        
        dR1H = sh * (nu * I1H) - 
          (f * R1H)

        dVC1H = (f * S1H) + (f * R1H)
        
        if(VC1H > v*pv1*(1-tau)){
          
          dS1H  = -(beta1HH * S1H * I1H) -
            (beta1HL * S1H * I1L)
          
          dR1H = sh * (nu * I1H)
          
          dVC1H = 0
          
        }
      }
      
      if(t > (time.rollout *1.5)) {
        
        ## now we are in rollout 
        ## Option 1: we have more vaccine than number at high risk 
        ## Loop 1: Do we have less high risk individuals then vaccine? 
        ## Loop 2: Once all high risk get vaccine, give the left over vaccines
        ## to the low risk 
        ## Loop 3: stop vaccinating low risk when we run out of (effective) vaccines 
        
        if(N2*ph2*f2 <= v*(1-pv1)){
          
          dS2H  = -(beta2HH * S2H * I2H) -
            (beta2HL * S2H * I2L) - 
            (f * S2H)
          
          dR2H = sh * (nu * I2H) - 
            (f * R2H)
          
          dVC2H = (f * S2H) + (f * R2H)
          
          if(VC2H > N2*ph2*f2*(1-tau)){
            
            dS2H  = -(beta2HH * S2H * I2H) -
              (beta2HL * S2H * I2L)
            
            dR2H = sh * (nu * I2H)
            
            dVC2H = 0
            
            dS2L = -(beta2LL * S2L * I2L) - 
              (beta2LH * S2L * I2H) - 
              (f*S2L) 
            
            dR2L =  sl* (nu * I2L) - 
              (f*R2L)
            
            dVC2L = (f*S2L) + (f*R2L) 
            
            if((VC2H + VC2L) > (v * (1-pv1)*(1-tau))){
              
              dS2L = -(beta2LL * S2L * I2L) -
                (beta2LH * S2L * I2H)
              
              dR2L = sl * (nu * I2L)
              
              dVC2L = 0
            }
          }
        }
        
        ## Option 2: we have less vaccine than number at high risk 
        ## Loop 1: check we have more high risk than vaccine 
        ## Loop 2: did we use up all of our vaccines 
        
        if(v*(1-pv1) < (N2* ph2 * f2)){
          
          dS2H  = -(beta2HH * S2H * I2H) -
            (beta2HL * S2H * I2L) - 
            (f * S2H)

          dR2H = sh * (nu * I2H) - 
            (f * R2H)
          
          dVC2H = (f * S2H) + (f * R2H)
          
          if(VC2H > v*(1-pv1)*(1-tau)){
            
            dS2H  = -(beta2HH * S2H * I2H) -
              (beta2HL * S2H * I2L)
            
            dR2H = sh * (nu * I2H)
            
            dVC2H = 0
            
          }
        }
      } 
      
    }
    
    states = c(dS1H, dE1H, dI1H, dR1H, dC1H, dV1H, dVC1H, dD1H,
               dS1L, dE1L, dI1L, dR1L, dC1L, dV1L, dVC1L, dD1L, 
               dS2H, dE2H, dI2H, dR2H, dC2H, dV2H, dVC2H, dD2H,
               dS2L, dE2L, dI2L, dR2L, dC2L, dV2L, dVC2L, dD2L)
    
    list(states)
  })
}

###########################
## Make Run
###########################

highrisk_mortality_rollout_varytime <- data.frame()

## proportion vaccinated pv_list
## vaccine seq(0, 1000, by = 100)
## proportion high v low risk seq(0,1,0.1) ph_list
## underlying immunity phi_list
v_list <- c(2e5, 4e5, 6e5, 8e5,10e5)
pv_list <- seq(0,1, 0.1)
rollout_list <- c(1, 10,30,50,100)
f_list <- c(0.01, 0.02, 0.03)
r_list <- c(2,4,8,16)
ph_list <- c(0.25, 0.5)
pct <- N1 * (1/1000)


for(v in v_list){
  for(pv1 in pv_list){
    for(time.rollout in rollout_list){
      for(f in f_list){
        for(ph1 in ph_list){
          for(r.val in r_list){
            
            pct_hr <- (N1 *ph1) * (1/1000)
            pct_lr <- (N1 *(1-ph1)) * (1/1000)
            
            ph2 = ph1
            
            RHH = r.val 
            RHL = r.val
            RLH = r.val
            RLL = r.val
            
            ## use params from kissler et al. 
            params <- c(beta1HH = RHH/(5*N1*ph1), 
                        beta1HL = RHL/(5*N1*ph1),
                        beta1LH = RLH/(5*N1*(1-ph1)),
                        beta1LL = RLL/(5*N1*(1-ph1)),
                        beta2HH = RHH/(5*N2 *ph2),
                        beta2HL = RHL/(5*N2*ph2),
                        beta2LH = RLH/(5*N2*(1-ph2)),
                        beta2LL = RLL/(5*N2*(1-ph2)),
                        sigma = 1/3, ## 1/latent period
                        nu = 1/5, ## 1/infectious period)
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
          
            inits<-c(
              S1H = (N1 * ph1) -pct_hr,
              E1H = 0, 
              I1H = pct_hr,
              R1H = 0,
              C1H = pct_hr, 
              V1H = 0,
              VC1H = 0,
              D1H = 0,
              
              S1L = (N1 * (1-ph1)) - pct_lr, 
              E1L = 0,
              I1L = pct_lr, 
              R1L = 0,
              C1L = pct_lr, 
              V1L = 0,
              VC1L = 0,
              D1L = 0,
              
              S2H = (N2* ph2) - pct_hr,
              E2H = 0, 
              I2H = pct_hr,
              R2H = 0, 
              C2H = pct_hr, 
              V2H = 0,
              VC2H = 0,
              D2H = 0,
              
              S2L = (N2 *(1-ph2)) - pct_lr,
              E2L = 0, 
              I2L = pct_lr,
              R2L = 0, 
              C2L = pct_lr,
              V2L = 0,
              VC2L = 0,
              D2L = 0
            )
            
            output2 <- ode(y = inits,
                           times = dt,
                           func = SEIR_AN_HIGHLOW_DEATHS_ro_vt,
                           method = "adams",
                           parms = params)
            
            out2.df <- data.frame(output2)[length(dt),]
            out2.df$v = v
            out2.df$pv = pv1
            out2.df$ph1 = ph1
            out2.df$ph2 = ph2
            out2.df$ro = time.rollout
            out2.df$f = f 
            out2.df$R = r.val
            highrisk_mortality_rollout_varytime = rbind(highrisk_mortality_rollout_varytime,out2.df)
            
            if(nrow(highrisk_mortality_rollout_varytime) %%100 == 0){
              print(nrow(highrisk_mortality_rollout_varytime))
            }
          }
        }
      }
    }
  }
}

### check what is going on with the Ns 
highrisk_mortality_rollout_varytime <- highrisk_mortality_rollout_varytime%>% 
  rowwise() %>% 
  mutate(N1 = sum(c(S1H, VC1H, E1H, I1H, R1H, D1H, 
                    S1L, VC1L, E1L, I1L, R1L, D1L)),
         N1H = sum(c(S1H, VC1H, E1H, I1H, R1H, D1H)),
         N1L = sum(c(S1L, VC1L, E1L, I1L, R1L, D1L)),
         N2 = sum(c(S2H, VC2H, E2H, I2H, R2H, D2H, 
                    S2L, VC2L, E2L, I2L, R2L, D2L)),
         N2H = sum(c(S2H, VC2H, E2H, I2H, R2H, D2H)),
         N2L = sum(c(S2L, VC2L, E2L, I2L, R2L, D2L))) %>%
  mutate(vtot = VC1H + VC1L + VC2H + VC2L,
         vtau = v * (1-tau),
         vtotH = VC1H + VC2H,
         vtotL = VC1L + VC2L)

save(highrisk_mortality_rollout_varytime , 
     file= "rollout_scenarios/data/highrisk_mortality_rollout_varytime.RData")

############################
## load & subset the data 
############################
load("rollout_scenarios/data/highrisk_mortality_rollout_varytime.RData")

highrisk_mortality_rollout_varytime$f.2 <- factor(highrisk_mortality_rollout_varytime$f, 
                                                     levels = c('0.03', '0.02', '0.01'))

highrisk_mortality_rollout_varytime %>% 
  filter(time == max(time)) %>% 
  mutate(D1D2 = (D1H + D2H + D1L + D2L)) %>% 
  mutate(D1D2H = (D1H + D2H), 
          D1D2L = (D1L + D2L))  -> dat.highrisk.m.rollout.vt

dat.highrisk.m.rollout.vt %>%
  count(R)

max_list_ro_hd.vt <- length(v_list) *
  length(pv_list) *
  length(rollout_list) * 
  length(f_list) *
  length(ph_list)

dat.highrisk.m.rollout.vt$xaxis <- rep(1:max_list_ro_hd.vt, 
                     each = length(r_list))

#############################
## Figure S14
#############################
#my.col <- c( "#6BAED6", "#4292C6", "#2171B5",  "#08519C",  "#08306B", "#969696", "#737373", "#525252", "#252525",  "#000000")
my.col.update <- c("#FF410DFF", "#6EE2FFFF", "#5B1A18", "#748AA6FF", "#F79D1EFF", "#D0DFE6FF",  "#02401B", "#FD6467", "#C27D38", "#C6CDF7" )

speed <- as_labeller(c(`0.03` = "Populations 1 & 2: 3% vaccinated/day", 
                       `0.02` = "Populations 1 & 2: 2% vaccinated/day", 
                       `0.01` = "Populations 1 & 2: 1% vaccinated/day"))


rollout <- as_labeller(c(`1` = "Population 1: 1 days \n Population 2: 1.5 days", 
                         `10` = "Population 1: 10 days \n Population 2: 15 days", 
                         `30` = "Population 1: 30 days \n Population 2: 45 days", 
                         `50` = "Population 1: 50 days \n Population 2: 75 days", 
                         `100` = "Population 1: 100 days \n Population 2: 150 days"))


rollout_hr_m.vt <- dat.highrisk.m.rollout.vt %>% 
  group_by(v, pv, D1D2, ro, f.2, R, xaxis) %>% 
  filter(ph1 == 0.5) %>% 
  filter(ro != 100) %>%
  summarize() %>%
  ggplot(aes(x = xaxis, 
             y = D1D2, 
             col = as.factor(v),  
             linetype = as.factor(R))) +
  facet_wrap(ro ~ f.2, 
             labeller = labeller(f.2 = as_labeller(speed), 
                                 ro = as_labeller(rollout)),
             ncol = 3) + 
  geom_line() + 
  theme_light() +
  labs(x = "Combination of total number of vaccine doses & proportion given to population 1", 
       y = "Cumulative number of deaths \n in population 1 & 2", 
       color = "Vaccine Doses",
       linetype = "R0 Value") + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.box = "vertical") + 
  ggtitle("50% High Mortality in both Populations") + 
  coord_cartesian(ylim = c(0, 6000), clip = "off") + 
  scale_color_manual(values = my.col.update) + 
  scale_y_continuous(breaks= seq(0,6000,2000), labels = seq(0,6000,2000)) 

rollout_hr_m.vt

dev.copy2pdf(file = "rollout_scenarios/figures/FigureS14.pdf", useDingbats=FALSE, 
             width = 8, height = 11)