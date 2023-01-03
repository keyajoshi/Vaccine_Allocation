##########################################
## This script generates the data 
## for the right panel of Figure S1
## considering the scenario where we have
## a leaky vaccine in a homogeneous population
## with no underlying immunity 
##########################################

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
tau <- 1-ve ## vaccine failure
phi1 <- 0 ## underlying immunity population 1
phi2 <- 0 ## underlying immunity population 2
years <- 3 ## duration of simulation 

dt <- seq(0, 365 *years, 1)

##########################################
## Function for two populations 
##########################################
SEIR.twopop.l <- function(t, x, params){
  
  with(as.list(c(params,x)),{
    
    dS1 = -(beta1 * S1 * I1)
    dV1 = -(tau * beta1 * V1 * I1)
    dE1 = (beta1 * S1 * I1) - 
      (sigma * E1) + 
      (tau * beta1 * V1 * I1)
    dI1 = (sigma *E1) - (nu * I1)
    dR1 = (nu * I1) 
    dC1 = (sigma *E1) ## tracks cumulative infections
    
    dS2 = -(beta2 * S2 * I2)
    dV2 = -(tau * beta2 * V2 * I2)
    dE2 = (beta2 * S2 * I2) - 
      (sigma * E2) + 
      (tau * beta2 * V2 * I2)
    dI2 = (sigma *E2) - (nu * I2)
    dR2 = (nu * I2) 
    dC2 = (sigma *E2) ## tracks cumulative infections 
    
    states = c(dS1, dV1, dE1, dI1, dR1, dC1, 
               dS2, dV2, dE2, dI2, dR2, dC2)
    
    list(states)
  })
} 

##########################################
## Make Run
##########################################

leaky_vaccine <- data.frame()

v_list <-seq(0, 1e6, 1e5)
pv_list <- seq(0, 1, 0.1)
phi_list <- 0 #
r_list <- c(2,4,8,16)
pct <- N1 * (1/1000)

for(v in v_list){
  for(pv1 in pv_list){
    for(phi1 in phi_list){
      for(R01 in r_list){
        
        R02 = R01 ## R0 is beta*N/gamma
        
        ## use params from kissler et al. 
        params <- c(sigma = 1/3, ## 1/latent period 
                    nu = 1/5, ## 1/infectious period
                    beta1 = R01/(5*N1), ## per contact probability of transmission
                    beta2 = R02/(5*N2), 
                    mu = 1/80) 
        
        inits <-c(
          S1 = N1 * (1-phi1) - v * pv1 * (1-phi1) - pct,
          V1 = v* pv1 * (1-phi1),
          E1 = 0,
          I1 = pct,
          R1 = N1 * phi1,
          C1 = pct,
          
          S2 = N2 * (1-phi2) - v * (1-pv1) * (1-phi2) - pct,
          V2 = v *(1-pv1)* (1-phi2),
          E2 = 0,
          I2 = pct,
          R2 = N2 * phi2,
          C2 = pct)
        
        output2 <- lsoda(y = inits,
                         times = dt,
                         func = SEIR.twopop.l,
                         parms = params)
        out2.df <- data.frame(output2) [length(dt),]
        out2.df$pv = pv1
        out2.df$v = v
        out2.df$phi1 = phi1
        out2.df$phi2 = phi2
        out2.df$R01 = R01
        out2.df$R02 = R02
        leaky_vaccine = rbind(leaky_vaccine,out2.df)
        
        if(nrow(leaky_vaccine) %%100 == 0){
          print(nrow(leaky_vaccine))}
      }
    }
  }
}

save(leaky_vaccine, file = "prophylactic_scenarios/data/leaky_vaccine.RData")

######################################
## plots
######################################
## load the data: 
#load("two_pop_immunity.Rdata")

## setting colors 
my.col <- c( "#6BAED6", "#4292C6", "#2171B5",  "#08519C",  "#08306B", "#969696", "#737373", "#525252", "#252525",  "#000000")

max_time_l<- twopop.l_threeyear %>% 
  filter(time == max(time)) %>% 
  mutate(C1C2 = C1 + C2)

max_time_l %>%
  count(R01)

max_list_l <- length(pv_list) * length(v_list)

max_time_l %>%
  mutate(xaxis = rep(1:max_list_l, each = length(r_list))) -> max_time_l

max_time_l %>%
  filter(v == 2e5| 
           v == 4e5 | 
           v == 6e5 | 
           v == 8e5 | 
           v == 10e5) -> max_time_l_sub

#######################################
## Leaky Vaccine 
#######################################
fig1leaky <- max_time_l_sub %>% 
  filter(phi1 == 0) %>% 
  group_by(v, pv, C1C2, R01, xaxis) %>% 
  summarize() %>% 
  ggplot(aes(x = xaxis, 
             y = C1C2, 
             col = as.factor(v),
             linetype = as.factor(R01))) + 
  geom_line() + 
  theme_light() + 
  labs(x = "combination of total number of vaccine and \n proportion given to population 1", 
       y = "cumulative number of cases \n in populations 1 and 2", 
       title = "Leaky vaccine",
       color = "Vaccine Doses") +
  scale_color_manual(values = my.col) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.title = element_text(size = 10))  +
  scale_y_continuous(breaks = seq(0,2e6, 4e5),
                     labels = seq(0,2e6, 4e5),
                     limits = c(0, 2e6))

fig1leaky
dev.copy2pdf(file = "Leakyvaccine_20210601.pdf", useDingbats=FALSE, 
             width = 6, height = 5)


