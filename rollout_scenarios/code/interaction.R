##########################################
## This script generates Figure S20.
## This scenario considers interaction
## between the two population and symmetric
## vaccine roll out 
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
## Load libraries & set options 
##########################################
library(deSolve)
library(tidyverse)
library(grid)
library(gridExtra)
library(viridis)
library(ggpubr)
library(ggnewscale)
library(gtable)
options(scipen = 999)

##########################################
## Parameters 
##########################################
N1 <- 1e6 ## population 1 size
N2 <- N1 ## population 2 size
ve <- 0.95 ## vaccine efficacy
tau <- 1-ve ## vaccine failure
phi1 <- 0 ## underlying immunity population 1
phi2 <- phi1 ## underlying immunity population 2
years <- 3 ## duration of simulation 

dt <- seq(0, 365 *years, 1)

##########################################
## Function for interaction rollout
##########################################
SEIR.twopop.nb.migr.ro <- function(t, x, params){
  
  with(as.list(c(params,x)),{
    
    N1 = S1 + VC1 +  E1 + I1 + R1
    
    dS1  = -(beta1* S1* ((1-i)*I1 + i*I2))
    dV1 = 0 
    dE1 = (beta1*S1*((1-i)*I1 + i*I2))- (sigma*E1)
    dI1 = (sigma*E1) - (nu*I1)
    dR1 = (nu*I1)
    dC1 = (sigma*E1)
    dVC1 = 0  
    
    if(t > time.rollout){
      dS1  = -(beta1* S1* ((1-i)*I1 + i*I2)) - (f * S1)
      
      dR1 = (nu*I1) - (f * R1)
      
      dVC1 = (f * S1) + (f * R1)
      
      if(VC1 > v * pv1* (1-tau) * (1-phi1)){
        
        dS1  = -(beta1* S1* ((1-i)*I1 + i*I2))
        
        dR1 = (nu*I1)
        
        dVC1 = 0
      }
    }
    
    N2 = S2 + VC2 + E2 + I2 + R2 
    
    dS2  = -(beta2* S2* ((1-i)*I2 + i*I1)) 
    dV2 = 0 
    dE2 = (beta2*S2*((1-i)*I2 + i*I1)) - (sigma*E2)
    dI2 = (sigma*E2) - (nu*I2)
    dR2 = (nu*I2)
    dC2 = (sigma*E2)
    dVC2 = 0 
    
    if(t > time.rollout){
      
      dS2  = -(beta2* S2* ((1-i)*I2 + i*I1)) - (f * S2)
      
      dR2 = (nu*I2) - (f * R2)
      
      dVC2 = (f * S2) + (f *R2)
      
      if(VC2 > v *(1-pv1) * (1-tau) * (1-phi2)){
        
        dS2  = -(beta2* S2* ((1-i)*I2 + i*I1))
        
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

interaction_rollout <- data.frame()

v_list <- seq(0, 1e6, 1e5)
pv_list <- seq(0, 1, 0.1)
#phi_list <- 0 
time.rollout <- 10 
f <-  0.02 
r_list <- c(2,4,8,16)
pct <- N1 * (1/1000)
i_list <- c(0, 0.001, 0.01, 0.1, 0.3, 0.5) # i can range between 0 and 0.5


for(v in v_list){
  for(pv1 in pv_list){
    for(i in i_list){
      for(R01 in r_list){
        
        R02 = R01 
        
        ## use params from kissler et al. 
        params <- c(sigma = 1/3, ## 1/latent period 
                    nu = 1/5, ## 1/infectious period
                    beta1 = R01/(5*N1), ## per contact probability of transmission
                    beta2 = R02/(5*N2), 
                    mu = 1/80) 
        
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
          VC2 = 0)
        
        output2 <- ode(method = "adams",
                       y = inits,
                       times = dt,
                       func = SEIR.twopop.nb.migr.ro, 
                       parms = params)
        out2.df <- data.frame(output2)[length(dt),]
        out2.df$pv = pv1
        out2.df$v = v
        out2.df$tau = tau
        out2.df$phi1 = phi1
        out2.df$phi2 = phi2
        out2.df$R01 = R01
        out2.df$R02 = R02
        out2.df$i = i
        out2.df$tr = time.rollout
        interaction_rollout = rbind(interaction_rollout,out2.df)
        
        if(nrow(interaction_rollout) %%100 == 0){
          print(nrow(interaction_rollout))}
        
      }
    }
  }
}


interaction_rollout <- interaction_rollout %>% 
  rowwise() %>% 
  mutate(N1 = sum(c(S1, E1, I1, R1, V1)),
         N2 =sum(c(S2, E2, I2, R2, V2)))


 save(interaction_rollout, file = "rollout_scenarios/data/interaction_rollout.RData")


######################################
## load & subset the data 
######################################
## load the data: 
load("rollout_scenarios/data/interaction_rollout.RData")

dat.interaction.rollout <- interaction_rollout %>% 
  filter(time == max(time)) %>% 
  mutate(C1C2 = C1 + C2)

dat.interaction.rollout %>%
  count(R01)

max_list_int_ro <- length(v_list) * length(pv_list) 

dat.interaction.rollout$xaxis <- rep(1:max_list_int_ro, each = length(r_list)*length(i_list))

dat.interaction.rollout %>% 
  filter(v == 2e5| 
           v == 4e5 | 
           v == 6e5 | 
           v == 8e5 | 
           v == 10e5) %>%
  filter(phi1 == 0) %>% 
  group_by(v, pv, C1C2, R01, i, xaxis) %>% 
  summarize() -> dat.interaction.rollout.sub

dfi0 <- dat.interaction.rollout.sub %>% 
  filter(i == 0)

dfi0.001<- dat.interaction.rollout.sub %>% 
  filter(i == 0.001)

dfi0.01 <- dat.interaction.rollout.sub %>% 
  filter(i == 0.01)

dfi0.1  <- dat.interaction.rollout.sub %>% 
  filter(i == 0.1)

dfi0.3  <- dat.interaction.rollout.sub %>% 
  filter(i == 0.3)

dfi0.5  <- dat.interaction.rollout.sub %>% 
  filter(i == 0.5)

###############
## Figure S20
###############

## setting colors 
#my.col <- c( "#6BAED6", "#4292C6", "#2171B5",  "#08519C",  "#08306B", "#969696", "#737373", "#525252", "#252525",  "#000000")
my.col.update <- c("#FF410DFF", "#6EE2FFFF", "#5B1A18", "#748AA6FF", "#F79D1EFF", "#D0DFE6FF",  "#02401B", "#FD6467", "#C27D38", "#C6CDF7" )
my.col.grey <- c("grey", "grey", "grey", "grey", "grey", "grey")


# Making the 'long' figure : create 5 different plots and assemble them.
plotlegend <- ggplot() + 
  geom_line(data = dfi0.001, 
            aes(x = xaxis,
                y = C1C2,
                col = as.factor(v),
                linetype = as.factor(R01))) + 
  scale_color_manual(values = my.col.update) +
  theme_light() + 
  labs(x = "Combination of total number of vaccine doses & \n proportion given to population 1", 
       # y = "Cumulative number of cases among \n in population 1 & 2",
       subtitle = "Interaction parameter = 0.001",
       color = "Vaccine doses",
       linetype = "R0 Value") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10), 
        plot.subtitle = element_text(hjust = 0.5, size = 12, face = "bold"),
        legend.title = element_text(size = 12)) + 
  coord_cartesian(ylim = c(0,2e6), clip = "off") +
  scale_y_continuous(breaks = seq(0,2e6, 4e5), labels = seq(0,2e6, 4e5)) +
  guides(linetype = guide_legend(order = 2),col = guide_legend(order = 1))

#
plot0.001 <- ggplot() + 
  geom_line(data = dfi0, 
            aes(x = xaxis, 
                y = C1C2, 
                col = as.factor(v),
                linetype = as.factor(R01))) + 
  scale_color_manual(values = my.col.grey, guide = FALSE) +
  new_scale_color() +
  geom_line(data = dfi0.001, 
            aes(x = xaxis, 
                y = C1C2, col = as.factor(v),
                linetype = as.factor(R01))) + 
  scale_color_manual(values = my.col.update, guide = FALSE) +
  theme_light() + 
  labs(x = "Combination of total number of vaccine doses & \n proportion given to population 1", 
       # y = "Cumulative number of cases among \n in population 1 & 2",
       subtitle = "Interaction parameter = 0.001",
       linetype = "R0 Value") +
  #color = "Vaccine doses") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10), 
        plot.subtitle = element_text(hjust = 0.5, size = 12, face = "bold"),
        legend.position = "none")  + 
  coord_cartesian(ylim = c(0,2e6), clip = "off") +
  scale_y_continuous(breaks = seq(0,2e6, 4e5), labels = seq(0,2e6, 4e5))

plot0.001

plot0.01 <- ggplot() + 
  geom_line(data = dfi0, 
            aes(x = xaxis, 
                y = C1C2, 
                col = as.factor(v),
                linetype = as.factor(R01))) + 
  scale_color_manual(values = my.col.grey, guide = FALSE) +
  new_scale_color() +
  geom_line(data = dfi0.01, 
            aes(x = xaxis, 
                y = C1C2, col = as.factor(v),
                linetype = as.factor(R01))) + 
  scale_color_manual(values = my.col.update, guide = FALSE) +
  theme_light() + 
  labs(x = "Combination of total number of vaccine doses & \n proportion given to population 1", 
       # y = "Cumulative number of cases among \n in population 1 & 2",
       subtitle = "Interaction parameter = 0.01",
       linetype = "R0 Value") +
  #color = "Vaccine doses") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10), 
        plot.subtitle = element_text(hjust = 0.5, size = 12, face = "bold"),
        legend.position = "none")  + 
  coord_cartesian(ylim = c(0,2e6), clip = "off") +
  scale_y_continuous(breaks = seq(0,2e6, 4e5), labels = seq(0,2e6, 4e5))

plot0.01

plot0.1 <- ggplot() + 
  geom_line(data = dfi0, 
            aes(x = xaxis, 
                y = C1C2, 
                col = as.factor(v),
                linetype = as.factor(R01))) + 
  scale_color_manual(values = my.col.grey, guide = FALSE) +
  new_scale_color() +
  geom_line(data = dfi0.1, 
            aes(x = xaxis, 
                y = C1C2, col = as.factor(v),
                linetype = as.factor(R01))) + 
  scale_color_manual(values = my.col.update, guide = FALSE) +
  theme_light() + 
  labs(x = "Combination of total number of vaccine doses & \n proportion given to population 1", 
       # y = "Cumulative number of cases among \n in population 1 & 2",
       subtitle = "Interaction parameter = 0.1",
       linetype = "R0 Value") +
  #color = "Vaccine doses") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10), 
        plot.subtitle = element_text(hjust = 0.5, size = 12, face = "bold"),
        legend.position = "none") + 
  coord_cartesian(ylim = c(0,2e6), clip = "off") +
  scale_y_continuous(breaks = seq(0,2e6, 4e5), labels = seq(0,2e6, 4e5))

plot0.1

plot0.3 <- ggplot() + 
  geom_line(data = dfi0, 
            aes(x = xaxis, 
                y = C1C2, 
                col = as.factor(v),
                linetype = as.factor(R01))) + 
  scale_color_manual(values = my.col.grey, guide = FALSE) +
  new_scale_color() +
  geom_line(data = dfi0.3, 
            aes(x = xaxis, 
                y = C1C2, col = as.factor(v),
                linetype = as.factor(R01))) + 
  scale_color_manual(values = my.col.update, guide = FALSE) +
  theme_light() + 
  labs(x = "Combination of total number of vaccine doses & \n proportion given to population 1", 
       # y = "Cumulative number of cases among \n in population 1 & 2",
       subtitle = "Interaction parameter = 0.3",
       linetype = "R0 Value") +
  #color = "Vaccine doses") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10), 
        plot.subtitle = element_text(hjust = 0.5, size = 12, face = "bold"),
        legend.position = "none") +
  coord_cartesian(ylim = c(0,2e6), clip = "off") +
  scale_y_continuous(breaks = seq(0,2e6, 4e5), labels = seq(0,2e6, 4e5))

plot0.3

plot0.5 <- ggplot() + 
  geom_line(data = dfi0, 
            aes(x = xaxis, 
                y = C1C2, 
                col = as.factor(v),
                linetype = as.factor(R01))) + 
  scale_color_manual(values = my.col.grey, guide = FALSE) +
  new_scale_color() +
  geom_line(data = dfi0.5, 
            aes(x = xaxis, 
                y = C1C2, col = as.factor(v),
                linetype = as.factor(R01))) + 
  scale_color_manual(values = my.col.update, guide = FALSE) +
  theme_light() + 
  labs(x = "Combination of total number of vaccine doses & \n proportion given to population 1", 
       # y = "Cumulative number of cases among \n in population 1 & 2",
       subtitle = "Interaction parameter = 0.5") +
  #color = "Vaccine doses") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10), 
        plot.subtitle = element_text(hjust = 0.5, size = 12, face = "bold"),
        legend.position = "none") + 
  coord_cartesian(ylim = c(0,2e6), clip = "off") +
  scale_y_continuous(breaks = seq(0,2e6, 4e5), labels = seq(0,2e6, 4e5))

plot0.5


# Join all 5 plots in one figure 

interactionlegend <- get_legend(plotlegend)
plotInteraction <- ggarrange(plot0.001, plot0.01, plot0.1, plot0.3, plot0.5, 
                             interactionlegend, nrow = 2, ncol = 3) 
#plotInteraction
annotate_figure(plotInteraction,
                left = text_grob("Cumulative number of cases among in populations 1 and 2",
                                 size = 12,  rot = 90),
                bottom = text_grob("Combination of total number of vaccine doses & proportion given to population 1", 
                                   size = 12))

dev.copy2pdf(file = "rollout_scenarios/figures/FigureS20.pdf", useDingbats=FALSE,
             width = 10, height = 7)

#










########################################
################################################################################
########################################
dat.immunity.R0 <- twopop_immunity_threeyear_R0 %>% 
  filter(time == max(time)) %>% ## number rows should remain same 
  mutate(C1C2 = C1 + C2)

head(dat.immunity.R0)

## for the first figure only look at situations where there is no 
## underlying immunity for either population 

dat.immunity.R0 %>%
  filter(v == 2e5| 
           v == 4e5 | 
           v == 6e5 | 
           v == 8e5 | 
           v == 10e5) -> dat.immunity.R0.sub 

#######################################
## Figure 1: no underlying immunity 
#######################################

tmp <- dat.immunity.R0.sub %>% 
  filter(phi1 == 0) %>% 
  group_by(v, pv, C1C2, R01) %>% 
  summarize() %>% 
  ggplot(aes(seq(0,1000, length.out = nrow(.)), 
             y = C1C2, 
             col = as.factor(v),
             linetype = as.factor(R01))) + 
  geom_line() + 
  theme_light() + 
  labs(x = "Combination of total number of vaccine doses and \n proportion given to population 1", 
       y = "Cumulative number of cases \n in populations 1 and 2", 
       #title = "All-or-nothing vaccine",
       color = "Vaccine Doses",
       linetype = "R0 Value") +
  scale_color_manual(values = my.col) +
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
  geom_segment(aes(x = 0, y = 1.95e6, xend = 180, yend = 1.95e6),
               arrow = arrow(length = unit(0.1, "cm"), type = "closed"), inherit.aes = F) + 
  annotate("text", x = c(10,170), y=2.1e6, label = c("0%", "100%"), size = 3) +
  annotate("text", x = 90, y = 2.25e6, label = "population 1", size = 3)

ggarrange(tmp, legend = "bottom")

# dev.copy2pdf(file = "Figure1_addR0.pdf", useDingbats=FALSE,
#              width = 6, height = 5)
