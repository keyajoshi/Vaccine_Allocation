##########################################
## This script generates figures S17 and 
## S18. Here we allow the vaccine efficacy
## to vary fixing R0 at 2. We also allow
## for underlying immunity to be different
## between the two populations 
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
library(ggsci)
options(scipen = 999)

##########################################
## Parameters 
##########################################
N1 <- 1e6 ## population 1 size 
N2 <- N1 ## population 2 size
#ve <- 0.95
#tau <- 1-ve #0.1
#phi1 <- 0 ## 
#phi2 <- 0 
years <- 3

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

vaccine_efficacy <- data.frame()

v_list <- seq(0, 1e6, 1e5)
pv_list <- seq(0, 1, 0.01)
phi_list <- seq(0, 0.5, 0.1)
pct <- N1 * (1/1000)
ve_list <- seq(0.5, 0.9, 0.05)
R01 <- 2
R02 <- 2

for(v in v_list){
  for(pv1 in pv_list){
    for(phi1 in phi_list){
      for(ve in ve_list){
        
        tau <- 1-ve
        phi2 <- phi1
        
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
        out2.df$ve = ve
        out2.df$tau = tau
        vaccine_efficacy = rbind(vaccine_efficacy,out2.df)
        
        if(nrow(vaccine_efficacy) %%100 == 0){
          print(nrow(vaccine_efficacy))}
        
      }
    }
  }
}


save(vaccine_efficacy, file = "prophylactic_scenarios/data/vaccine_efficacy.RData")

######################################
## load & subset the data
######################################
## load the data: 
load("prophylactic_scenarios/data/vaccine_efficacy.RData")


dat.immunity.efficacy <- vaccine_efficacy %>% 
  filter(time == max(time)) %>% 
  mutate(C1C2 = C1 + C2)

head(dat.immunity.efficacy)

dat.immunity.efficacy %>% 
  count(ve)
  
max_list <- length(pv_list) * 
  length(v_list) * length(phi_list)

dat.immunity.efficacy$xaxis <- rep(1:max_list, each = length(ve_list))

## for the first figure only look at situations where there is no 
## underlying immunity for either population 

dat.immunity.efficacy %>%
  filter(v == 2e5| 
           v == 4e5 | 
           v == 6e5 | 
           v == 8e5 | 
           v == 10e5) -> dat.immunity.efficacy.sub 


#######################################
## Figure S17: no underlying immunity 
#######################################
## setting colors 
#my.col <- c( "#6BAED6", "#4292C6", "#2171B5",  "#08519C",  "#08306B", "#969696", "#737373", "#525252", "#252525",  "#000000", "#525252", "#252525")
my.col.update <- c("#FF410DFF", "#6EE2FFFF", "#5B1A18", "#748AA6FF", "#F79D1EFF", "#D0DFE6FF")

ve_phi0 <- dat.immunity.efficacy.sub %>% 
  filter(phi1 == 0) %>% 
  filter(ve %in% c(0.5,0.6, 0.7, 0.8, 0.9)) %>%
  group_by(v, pv, C1C2, ve, xaxis) %>% 
  summarize() %>% 
  ggplot(aes(x = xaxis, 
             y = C1C2, 
             col = as.factor(v),
             linetype = as.factor(ve))) + 
  geom_line() + 
  theme_light() + 
  labs(x = "Combination of total number of vaccine doses and \n proportion given to population 1", 
       y = "Cumulative number of cases \n in populations 1 and 2", 
       #title = "All-or-nothing vaccine",
       color = "Vaccine Doses",
       linetype = "Vaccine Efficacy Value") +
  scale_color_manual(values = my.col.update) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.title = element_text(size = 10),
        legend.box = "vertical")  +
  scale_y_continuous(breaks = seq(0, 2e6, 4e5),
                     labels = seq(0,2e6, 4e5),
                     limits = c(0, 2.3e6)) + 
  coord_cartesian(ylim = c(0,2.3e6), clip = "off") + 
  ## don't need to copy these over as they are specifically for the arrow 
  geom_segment(aes(x = 1200, 
                   y = 1.95e6, 
                   xend = 1800, 
                   yend = 1.95e6),
               arrow = arrow(length = unit(0.1, 
                                           "cm"), 
                             type = "closed"), 
               inherit.aes = F) + 
  annotate("text", x = c(1210,1790), y=2.1e6, label = c("0%", "100%"), size = 3) +
  annotate("text", x = 1550, y = 2.25e6, label = "population 1", size = 3)

ggarrange(ve_phi0, legend = "bottom")

dev.copy2pdf(file = "prophylactic_scenarios/figures/FigureS17.pdf",
             useDingbats=FALSE,
              width = 6, height = 5)

#######################################
## Figure S18: underlying immunity 
#######################################
ve_phi0 <- dat.immunity.efficacy.sub %>% 
  filter(phi1 == 0,
         phi2 == 0) %>% 
  filter(ve %in% seq(0.5, 0.9, 0.1)) %>%
  group_by(v, pv, C1C2, ve, xaxis) %>% 
  summarize() %>% 
  ggplot(aes(x = xaxis, 
             y = C1C2, 
             col = as.factor(v),
             linetype = as.factor(ve))) + 
  geom_line() + 
  theme_light() + 
  labs(x = "Combination of total number of vaccine doses and \n proportion given to population 1", 
       y = "Cumulative number of cases \n in populations 1 and 2", 
       color = "Vaccine Doses",
       linetype = "Vaccine Efficacy Value") +
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
  ggtitle("0% immune in populations 1 and 2") + 
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold")) +
  scale_y_continuous(breaks = seq(0, 2e6, 4e5),
                     labels = seq(0,2e6, 4e5),
                     limits = c(0, 2.3e6)) + 
  coord_cartesian(ylim = c(0,2.3e6), clip = "off") +
  geom_segment(aes(x = 1200, 
                   y = 1.95e6, 
                   xend = 1800, 
                   yend = 1.95e6),
               arrow = arrow(length = unit(0.1, 
                                           "cm"), 
                             type = "closed"), 
               inherit.aes = F) + 
  annotate("text", x = c(1210,1790), y=2.1e6, label = c("0%", "100%"), size = 3) +
  annotate("text", x = 1550, y = 2.25e6, label = "population 1", size = 3)

ve_phi0_legend <- get_legend(ve_phi0)

ve_phi1 <- dat.immunity.efficacy.sub %>% 
  filter(phi1 == 0.1,
         phi2 == 0.1) %>% 
  filter(ve %in% seq(0.5, 0.9, 0.1)) %>%
  group_by(v, pv, C1C2, ve, xaxis) %>% 
  summarize() %>% 
  ggplot(aes(x = xaxis, 
             y = C1C2, 
             col = as.factor(v),
             linetype = as.factor(ve))) + 
  geom_line() + 
  theme_light() + 
  labs(x = "Combination of total number of vaccine doses and \n proportion given to population 1", 
       y = "Cumulative number of cases \n in populations 1 and 2", 
       color = "Vaccine Doses",
       linetype = "efficacy Value") +
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
  coord_cartesian(ylim = c(0,2.3e6), clip = "off") +
  geom_segment(aes(x = 1200, 
                   y = 1.95e6, 
                   xend = 1800, 
                   yend = 1.95e6),
               arrow = arrow(length = unit(0.1, 
                                           "cm"), 
                             type = "closed"), 
               inherit.aes = F) + 
  annotate("text", x = c(1210,1790), y=2.1e6, label = c("0%", "100%"), size = 3) +
  annotate("text", x = 1550, y = 2.25e6, label = "population 1", size = 3) + 
  ggtitle("10% immune in populations 1 and 2")

ve_phi2 <- dat.immunity.efficacy.sub %>% 
  filter(phi1 == 0.2,
         phi2 == 0.2) %>% 
  filter(ve %in% seq(0.5, 0.9, 0.1)) %>%
  group_by(v, pv, C1C2, ve, xaxis) %>% 
  summarize() %>% 
  ggplot(aes(x = xaxis, 
             y = C1C2, 
             col = as.factor(v),
             linetype = as.factor(ve))) + 
  geom_line() + 
  theme_light() + 
  labs(x = "Combination of total number of vaccine doses and \n proportion given to population 1", 
       y = "Cumulative number of cases \n in populations 1 and 2", 
       color = "Vaccine Doses",
       linetype = "efficacy Value") +
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
  coord_cartesian(ylim = c(0,2.3e6), clip = "off") +
  geom_segment(aes(x = 1200, 
                   y = 1.95e6, 
                   xend = 1800, 
                   yend = 1.95e6),
               arrow = arrow(length = unit(0.1, 
                                           "cm"), 
                             type = "closed"), 
               inherit.aes = F) + 
  annotate("text", x = c(1210,1790), y=2.1e6, label = c("0%", "100%"), size = 3) +
  annotate("text", x = 1550, y = 2.25e6, label = "population 1", size = 3) + 
  ggtitle("20% immune in populations 1 and 2")

ve_phi4 <-dat.immunity.efficacy.sub %>% 
  filter(phi1 == 0.4,
         phi2 == 0.4) %>% 
  filter(ve %in% seq(0.5, 0.9, 0.1)) %>%
  group_by(v, pv, C1C2, ve, xaxis) %>% 
  summarize() %>% 
  ggplot(aes(x = xaxis, 
             y = C1C2, 
             col = as.factor(v),
             linetype = as.factor(ve))) + 
  geom_line() + 
  theme_light() + 
  labs(x = "Combination of total number of vaccine doses and \n proportion given to population 1", 
       y = "Cumulative number of cases \n in populations 1 and 2", 
       color = "Vaccine Doses",
       linetype = "efficacy Value") +
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
  coord_cartesian(ylim = c(0,2.3e6), clip = "off") +
  geom_segment(aes(x = 1200, 
                   y = 1.95e6, 
                   xend = 1800, 
                   yend = 1.95e6),
               arrow = arrow(length = unit(0.1, 
                                           "cm"), 
                             type = "closed"), 
               inherit.aes = F) + 
  annotate("text", x = c(1210,1790), y=2.1e6, label = c("0%", "100%"), size = 3) +
  annotate("text", x = 1550, y = 2.25e6, label = "population 1", size = 3) + 
  ggtitle("40% immune in populations 1 and 2")

figs17 <-ggarrange(ve_phi0 + rremove("legend"), 
                 ve_phi1 + rremove("legend"), 
                 ve_phi2 + rremove("legend"), 
                 ve_phi4 + rremove("legend"), 
                 ncol=2, nrow=2)

figs17ann <- annotate_figure(figs17, 
                bottom = text_grob("Combination of total number of vaccine doses & proportion given to population 1", 
                                   size = 12),
                left = text_grob("Cumulative number of cases in populations 1 and 2",
                                 size = 12, 
                                 rot = 90))

figs17tot <- ggarrange(figs17ann, ve_phi0_legend$grobs[[1]],ve_phi0_legend$grobs[[2]] ,nrow = 3, heights = c(5,0.3, 0.3))
figs17tot

dev.copy2pdf(file = "prophylactic_scenarios/figures/FigureS18.pdf", useDingbats=FALSE,
             width = 11, height = 8)
