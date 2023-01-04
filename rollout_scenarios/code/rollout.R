##########################################
## This script generates figure 4 in the 
## main text. In this scenario we consider
## timing and speed of vaccine rollout. 
## In this script both the timing and 
## speed are symmetric across populations. 
##########################################

#########################################
## legend function
#########################################
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

##########################################
## Load libraries 
##########################################
library(deSolve)
library(tidyverse)
library(grid)
library(gridExtra)
library(viridis)
library(ggpubr)
options(scipen = 999)

##########################################
## Parameters 
##########################################
N1 <- 1e6 ## population 1 size 
N2 <- N1 ## population 2 size
ve <- 0.95 ## vaccine efficacy
tau <- 1-ve ## failure 
phi1 <- 0 ## underlying immunity population 1
phi2 <- 0 ## underlying immunity population 2 

years <- 3 ## duration of simulation 
dt <- seq(0, 365*years, 1)

##########################################
## Function for two populations 
##########################################
SEIR.twopop.aon.ro <- function(t, x, params){
  
  with(as.list(c(params,x)),{
    
    N1 = S1 + VC1 +  E1 + I1 + R1

    dS1  = -(beta1* S1* I1)
    dV1 = 0 
    dE1 = (beta1*S1*I1) - (sigma*E1)
    dI1 = (sigma*E1) - (nu*I1)
    dR1 = (nu*I1)
    dC1 = (sigma*E1)
    dVC1 = 0  
    
    if(t > time.rollout){
      dS1  = -(beta1* S1* I1) - (f * S1)
      
      dR1 = (nu*I1) - (f * R1)
      
      dVC1 = (f * S1) + (f * R1)
      
      if(VC1 > v * pv1* (1-tau)){
        
        dS1  = -(beta1* S1* I1)
        
        dR1 = (nu*I1)
        
        dVC1 = 0
      }
    }
    
    N2 = S2 + VC2 + E2 + I2 + R2 
    
    dS2  = -(beta2* S2* I2) 
    dV2 = 0 
    dE2 = (beta2*S2*I2) - (sigma*E2)
    dI2 = (sigma*E2) - (nu*I2)
    dR2 = (nu*I2)
    dC2 = (sigma*E2)
    dVC2 = 0 
    
    if(t > time.rollout){
      
      dS2  = -(beta2* S2* I2) - (f * S2)
      
      dR2 = (nu*I2) - (f * R2)
      
      dVC2 = (f * S2) + (f *R2)
      
      if(VC2 > v *(1-pv1) * (1-tau)){
        
        dS2  = -(beta2* S2* I2)
        
        dR2 = (nu*I2)
        
        dVC2 = 0 
      }
    }
    
    
    states = c(dS1, dV1, dE1, dI1, dR1, dC1, dVC1, 
               dS2, dV2, dE2, dI2, dR2, dC2, dVC2)
    
    return(list(states))
    
  })
}

##########################################
## Make Run
##########################################

rollout <- data.frame()

v_list <- seq(0, 1e6, 1e5)
pv_list <- seq(0, 1, 0.1)
rollout_list <- c(1,10,30,50, 100)
f_list <- c(0.01, 0.02, 0.03)
r_list <- c(2,4,8,16)
pct <- N1 * (1/1000)

for(v in v_list){
  for(pv1 in pv_list){
    for(time.rollout in rollout_list){
      for(f in f_list){
        for(R01 in r_list){
          
          R02 = R01 ## R0 is beta*N/gamma
          
          ## use params from kissler et al. 
          params <- c(sigma = 1/3, ## 1/latent period 
                      nu = 1/5, ## 1/infectious period
                      beta1 = R01/(N1*5), ## per contact probability of transmission
                      beta2 = R01/(N2*5)) 
          inits <-c(
            S1 = N1 - pct, 
            V1 = 0,
            E1 = 0, 
            I1 = pct,
            R1 = 0,  
            C1 = pct,
            VC1 = 0, 
            
            S2 = N2 - pct,
            V2 = 0,
            E2 = 0,
            I2 = pct, 
            R2 = 0,
            C2 = pct,
            VC2 =0)
          
          output2 <- ode(y = inits,
                         times = dt,
                         func = SEIR.twopop.aon.ro,
                         parms = params,
                         method = "adams")
          
          out2.df <- data.frame(output2)[length(dt),]
          out2.df$pv = pv1
          out2.df$v = v
          out2.df$phi1 = phi1
          out2.df$phi2 = phi2
          out2.df$ro = time.rollout
          out2.df$f = f 
          out2.df$R01 = R01
          rollout = rbind(rollout,out2.df)
          
          if(nrow(rollout) %%100 == 0){
            print(nrow(rollout))
          }
        }
      }
    }
    
  }
}

## check Ns
rollout <- rollout %>% 
  rowwise() %>% 
  mutate(N1 = sum(c(S1, VC1, E1, I1, R1)),
         N2 = sum(c(S2, VC2, E2, I2, R2))) %>%
  mutate(vtot = VC1 + VC2,
         vtau = v * (1-tau)) %>%
  mutate(C1C2 = C1 + C2)

save(rollout, file = "rollout_scenarios/data/rollout.RData")

######################################
## load & subset data 
######################################
load("rollout_code/data/varyingR0/rollout.RData")

rollout$f.2 <- factor(rollout$f, 
                      levels = c('0.03', '0.02', '0.01'))

dat.rollout <- rollout %>% 
  data.frame() %>%
  filter(time == max(time)) %>% 
  mutate(C1C2 = C1 + C2)

dat.rollout %>%
  count(R01)

max_list_ro <- length(v_list)*
  length(pv_list)*
  length(rollout_list)*
  length(f_list)

dat.rollout$xaxis <- rep(1:max_list_ro, each = length(r_list))

dat.rollout %>% 
filter(v == 2e5| 
         v == 4e5 | 
         v == 6e5 | 
         v == 8e5 | 
         v == 10e5) -> dat.rollout.sub

#######################################
## Figure 4
#######################################
## setting colors 
#my.col <- c( "#6BAED6", "#4292C6", "#2171B5",  "#08519C",  "#08306B", "#969696", "#737373", "#525252", "#252525",  "#000000")
my.col.update <- c("#FF410DFF", "#6EE2FFFF", "#5B1A18", "#748AA6FF", "#F79D1EFF", "#D0DFE6FF",  "#02401B", "#FD6467", "#C27D38", "#C6CDF7" )

speed <- as_labeller(c(`0.01` = "1% vaccinated/day", 
                       `0.02` = "2% vaccinated/day", 
                       `0.03` = "3% vaccinated/day"))

rollout <- as_labeller(c(`10` = "10 days", 
                         `50` = "50 days", 
                         `100` = "100 days"))

fig_ro <- dat.rollout.sub %>% 
  group_by(v, pv, C1C2, ro, f.2, R01, xaxis) %>% 
  filter(ro != 1) %>%
  filter(ro != 30) %>%
  summarize() %>% 
  ggplot(aes(x = xaxis, 
             y = C1C2, 
             col = as.factor(v),
             linetype = as.factor(R01))) +
  facet_wrap(ro~f.2, labeller = labeller(f.2 = as_labeller(speed), 
                                       ro = as_labeller(rollout)), ncol = 3) +
  geom_line() + 
  theme_light() + 
  labs(x = "Combination of total number of vaccine doses &  proportion given to population 1", 
       y = "Cumulative number of cases\n in population 1 & 2", 
       color = "Vaccine Doses",
       linetype = "R0") +
  scale_color_manual(values = my.col.update) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.box="vertical") + 
  coord_cartesian(ylim = c(0,2e6), clip = "off") +
  scale_y_continuous(breaks = seq(0,2e6, 4e5), labels = seq(0,2e6, 4e5))

fig_ro

dev.copy2pdf(file = "rollout_scenarios/figures/Figure4.pdf", useDingbats=FALSE, 
             width = 8, height = 8)
