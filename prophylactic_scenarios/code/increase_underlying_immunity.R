##########################################
## This script generates Figure S3.
## Here we increase the level of underlying
## immunity to be 50-80% in both populations
##########################################

#########################################
## useful function
#########################################
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
##########################################
## Load libraries & set options 
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
N2 <- N1 ## population 2 size (symmetric)
ve <- 0.95 ## vaccine efficacy
tau <- 1-ve ## failure
phi1 <- 0 ## underlying immunity population 1
phi2 <- 0 ## underlying immunity population 2
years <- 3 ## duration of simulation 

dt <- seq(0, 365 *years, 1)

##########################################
## Function for two populations 
##########################################
SEIR.twopop.aon <- function(t, x, params){
  
  with(as.list(c(params,x)),{
    
    
    N1 = S1 + E1 + I1 + R1 + V1
    N2 = S2 + E2 + I2 + R2 + V2

    dS1  = -(beta1* S1* I1)
    dV1 = 0 
    dE1 = (beta1*S1*I1)- (sigma*E1)
    dI1 = (sigma*E1) - (nu*I1)
    dR1 = (nu*I1)
    dC1 = (sigma*E1)
    
    dS2  = -(beta2* S2* I2) 
    dV2 = 0 
    dE2 = (beta2*S2*I2) - (sigma*E2)
    dI2 = (sigma*E2) - (nu*I2)
    dR2 = (nu*I2)
    dC2 = (sigma*E2)
    
    states = c(dS1, dV1, dE1, dI1, dR1, dC1, 
               dS2, dV2, dE2, dI2, dR2, dC2)
    
    return(list(states))
    
  })
}

##########################################
## Make Run
##########################################
increase_immunity <- data.frame()

v_list <- seq(0, 1e6, 1e5)
pv_list <- seq(0, 1, 0.1)
phi_list <- seq(0.5,0.9, 0.1)
r_list <- c(2,4,8,16)
pct <- N1 * (1/1000)

for(v in v_list){
  for(pv1 in pv_list){
    for(phi1 in phi_list){
      for(R01 in r_list){
        
       R02 = R01 
        phi2 = 0 
        ## use params from kissler et al. 
        params <- c(sigma = 1/3, ## 1/latent period 
                    nu = 1/5, ## 1/infectious period
                    beta1 = R01/(5*N1), ## per contact probability of transmission
                    beta2 = R02/(5*N2), 
                    mu = 1/80) 
        
    
        inits <-c(
          S1 = N1 * (1-phi1) - v * pv1* (1-tau) * (1-phi1) - pct, 
          V1 = v* pv1 * (1-tau) * (1-phi1),
          E1 = 0, 
          I1 = pct,
          R1 = N1 * phi1,  
          C1 = pct, 
          
          S2 = N2 * (1-phi2) - v * (1-pv1) * (1-tau) * (1-phi2) - pct,
          V2 = v *(1-pv1) * (1-tau) * (1-phi2),
          E2 = 0,
          I2 = pct, 
          R2 = N2 * phi2,
          C2 = pct)

        output2 <- lsoda(y = inits,
                         times = dt,
                         func = SEIR.twopop.aon,
                         parms = params)
        out2.df <- data.frame(output2)[length(dt),]
        out2.df$pv = pv1
        out2.df$v = v
        out2.df$tau = tau
        out2.df$phi1 = phi1
        out2.df$phi2 = phi2
        out2.df$R01 = R01
        out2.df$R02 = R02
        increase_immunity = rbind(increase_immunity,out2.df)
        if(nrow(increase_immunity) %%100 == 0){
          print(nrow(increase_immunity))}
      }
    }
  }
}

increase_immunity <- increase_immunity %>% 
  rowwise() %>% 
  mutate(N1 = sum(c(S1, E1, I1, R1, V1)),
         N2 =sum(c(S2, E2, I2, R2, V2)))

save(increase_immunity, file = "prophylactic_scenarios/data/increase_underlying_immunity.RData")

######################################
## load & subset the data 
######################################
## load the data: 
load("prophylactic_scenarios/data/increase_underlying_immunity.RData")

## setting colors 
#my.col <- c( "#6BAED6", "#4292C6", "#2171B5",  "#08519C",  "#08306B", "#969696", "#737373", "#525252", "#252525",  "#000000", "#525252", "#252525")
my.col.update <- c("#FF410DFF", "#6EE2FFFF", "#5B1A18", "#748AA6FF", "#F79D1EFF", "#D0DFE6FF")

dat.immunity.R0.inc <- increase_immunity %>% 
  filter(time == max(time)) %>% ## number rows should remain same 
  mutate(C1C2 = C1 + C2)

head(dat.immunity.R0.inc)

dat.immunity.R0.inc %>% 
  count(R01)
  
max_list <- length(pv_list) * 
  length(v_list) * 
  length(phi_list)

dat.immunity.R0.inc$xaxis = rep(1:max_list, each = length(r_list))

dat.immunity.R0.inc %>%
  filter(v == 2e5| 
           v == 4e5 | 
           v == 6e5 | 
           v == 8e5 | 
           v == 10e5) -> dat.immunity.R0.inc.sub 

#######################################
## Figure S3: increase immunity
#######################################
phi50 <- dat.immunity.R0.inc.sub %>% 
  filter(phi1 == 0.5,
         phi2 == 0) %>% 
  group_by(v, pv, C1C2, R01, xaxis) %>% 
  summarize() %>% 
  ggplot(aes(x = xaxis, 
             y = C1C2, 
             col = as.factor(v),
             linetype = as.factor(R01))) + 
  geom_line() + 
  theme_light() + 
  labs(x = "Combination of total number of vaccine doses and \n proportion given to population 1", 
       y = "Cumulative number of cases \n in populations 1 and 2", 
       color = "Vaccine Doses",
       linetype = "R0 Value") +
  scale_color_manual(values = my.col.update) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.box = "vertical",
        legend.position = "bottom")  +
  ggtitle("50% immune in population 1, 0% immune population 2") + 
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold")) +
  scale_y_continuous(breaks = seq(0, 2e6, 4e5),
                     labels = seq(0,2e6, 4e5),
                     limits = c(0, 2.3e6)) + 
  coord_cartesian(ylim = c(0,2.3e6), clip = "off") 

phi50_legend <- get_legend(phi50)

phi60 <- dat.immunity.R0.inc.sub %>% 
  filter(phi1 == 0.6,
         phi2 == 0) %>% 
  group_by(v, pv, C1C2, R01, xaxis) %>% 
  summarize() %>% 
  ggplot(aes(x = xaxis, 
             y = C1C2, 
             col = as.factor(v),
             linetype = as.factor(R01))) + 
  geom_line() + 
  theme_light() + 
  labs(x = "Combination of total number of vaccine doses and \n proportion given to population 1", 
       y = "Cumulative number of cases \n in populations 1 and 2", 
       color = "Vaccine Doses",
       linetype = "R0 Value") +
  scale_color_manual(values = my.col.update) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 10),
        legend.box = "vertical")  +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold")) +
  scale_y_continuous(breaks = seq(0, 2e6, 4e5),
                     labels = seq(0,2e6, 4e5),
                     limits = c(0, 2.3e6)) + 
  ggtitle("60% immune in population 1, 0% immune population 2")
  
phi70 <- dat.immunity.R0.inc.sub %>% 
  filter(phi1 == 0.7,
         phi2 == 0) %>% 
  group_by(v, pv, C1C2, R01, xaxis) %>% 
  summarize() %>% 
  ggplot(aes(x = xaxis, 
             y = C1C2, 
             col = as.factor(v),
             linetype = as.factor(R01))) + 
  geom_line() + 
  theme_light() + 
  labs(x = "Combination of total number of vaccine doses and \n proportion given to population 1", 
       y = "Cumulative number of cases \n in populations 1 and 2", 
       color = "Vaccine Doses",
       linetype = "R0 Value") +
  scale_color_manual(values = my.col.update) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 10),
        legend.box = "vertical")  +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold")) +
  scale_y_continuous(breaks = seq(0, 2e6, 4e5),
                     labels = seq(0,2e6, 4e5),
                     limits = c(0, 2.3e6)) + 
  ggtitle("70% immune in population 1, 0% immune population 2")
  
phi80 <-dat.immunity.R0.inc.sub %>% 
  filter(phi1 == 0.8,
         phi2 == 0) %>% 
  group_by(v, pv, C1C2, R01, xaxis) %>% 
  summarize() %>% 
  ggplot(aes(x = xaxis, 
             y = C1C2, 
             col = as.factor(v),
             linetype = as.factor(R01))) + 
  geom_line() + 
  theme_light() + 
  labs(x = "Combination of total number of vaccine doses and \n proportion given to population 1", 
       y = "Cumulative number of cases \n in populations 1 and 2", 
       color = "Vaccine Doses",
       linetype = "R0 Value") +
  scale_color_manual(values = my.col.update) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 10),
        legend.box = "vertical")  +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold")) +
    ggtitle("80% immune in population 1, 0% immune population 2") + 
  scale_y_continuous(breaks = seq(0, 2e6, 4e5),
                     labels = seq(0,2e6, 4e5),
                     limits = c(0, 2.3e6))

figS3 <-ggarrange(phi50 + rremove("legend"), 
                 phi60 + rremove("legend"), 
                 phi70 + rremove("legend"), 
                 phi80 + rremove("legend"), 
                 ncol=2, nrow=2)

figS3ann <- annotate_figure(figS3, 
                bottom = text_grob("Combination of total number of vaccine doses & proportion given to population 1", 
                                   size = 12),
                left = text_grob("Cumulative number of cases in populations 1 and 2",
                                 size = 12, 
                                 rot = 90))

figS3tot <- ggarrange(figS3ann, phi50_legend$grobs[[1]],phi50_legend$grobs[[2]] ,nrow = 3, heights = c(5,0.3, 0.3))
figS3tot

dev.copy2pdf(file = "prophylactic_scenarios/figures/FigureS3.pdf", useDingbats=FALSE,
             width = 11, height = 8)
