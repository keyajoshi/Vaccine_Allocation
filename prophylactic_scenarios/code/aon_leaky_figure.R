##########################################
## This script plots the leaky and 
## aon figure for the supplement 
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
options(scipen = 999)

#################################
## Code combining AoN and Leaky 
#################################
##load dataframes 
load("prophylactic_scenarios/data/underlying_immunity.RData")
load("prophylactic_scenarios/data/leaky_vaccine.RData")

######################################
## subset data 
######################################
r_list <- c(2,4,8,16)

#######
## aon vaccine 
#######

dat.immunity.R0 <- underlying_immunity %>% 
  filter(time == max(time)) %>% ## number rows should remain same 
  mutate(C1C2 = C1 + C2)

dat.immunity.R0 %>% 
  count(R01)

## need to manually specify because
## the length of the lists are different across 
## the two runs (pv is 0.1 v 0.01)
max_list <- length(seq(0, 1, 0.01)) * 
  length(seq(0, 1e6, 1e5)) * length(seq(0, 0.5, 0.1))

dat.immunity.R0$xaxis <- rep(1:max_list, each = length(r_list))

dat.immunity.R0 %>%
  filter(v == 2e5| 
           v == 4e5 | 
           v == 6e5 | 
           v == 8e5 | 
           v == 10e5) -> dat.immunity.R0.sub 

#######
## leaky vaccine 
#######

dat.leaky.R0 <- leaky_vaccine %>% 
  filter(time == max(time)) %>% 
  mutate(C1C2 = C1 + C2)

dat.leaky.R0 %>%
  count(R01)

max_list_l <- length(seq(0, 1, 0.1)) * 
  length(seq(0, 1e6, 1e5))

dat.leaky.R0$xaxis <- rep(1:max_list_l, each = length(r_list))

dat.leaky.R0 %>%
  filter(v == 2e5| 
           v == 4e5 | 
           v == 6e5 | 
           v == 8e5 | 
           v == 10e5) -> dat.leaky.R0.sub

##############################
## plotting
##############################
## set colors 
#my.col <- c( "#6BAED6", "#4292C6", "#2171B5",  "#08519C",  "#08306B", 
#"#969696", "#737373", "#525252", "#252525",  "#000000")

my.col.update <- c("#FF410DFF", "#6EE2FFFF", "#5B1A18", "#748AA6FF", "#F79D1EFF", "#D0DFE6FF",  "#02401B", "#FD6467", "#C27D38", "#C6CDF7" )

## generate leaky and AoN plots
fig1aon <- dat.immunity.R0.sub %>% 
  filter(phi1 == 0) %>%
  group_by(v, pv, C1C2, R01, xaxis) %>% 
  summarize() %>% 
  ggplot(aes(x = xaxis, 
             y = C1C2, 
             col = as.factor(v),
             linetype = as.factor(R01))) + 
  geom_line() + 
  theme_light() + 
  labs(title = "All-or-nothing vaccine",
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
        legend.position = "bottom",
        legend.box = "vertical") + 
  scale_y_continuous(breaks = seq(0,2e6, 4e5),
                     labels = seq(0,2e6, 4e5),
                     limits = c(0, 2e6)) + 
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold")) 


fig1leaky <-  dat.leaky.R0.sub %>% 
  filter(phi1 == 0) %>% 
  group_by(v, pv, C1C2, R01, xaxis) %>% 
  summarize() %>% 
  ggplot(aes(x = xaxis, 
             y = C1C2, 
             col = as.factor(v),
             linetype = as.factor(R01))) + 
  geom_line() + 
  theme_light() + 
  labs(title = "Leaky vaccine",
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
        legend.position = "bottom",
        legend.box = "vertical") + 
  scale_y_continuous(breaks = seq(0,2e6, 4e5),
                     labels = seq(0,2e6, 4e5),
                     limits = c(0, 2e6)) + 
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold")) 


aon_legend <- get_legend(fig1aon)

figS1 <-ggarrange(fig1aon + rremove("legend"), 
                 fig1leaky + rremove("legend"),
                 ncol=2)

figS1ann <- annotate_figure(figS1, 
                           bottom = text_grob("Combination of total number of vaccine doses & proportion given to population 1", 
                                              size = 12),
                           left = text_grob("Cumulative number of cases in populations 1 and 2",
                                            size = 12, 
                                            rot = 90))
figS1ann
figS1tot <- ggarrange(figS1ann, aon_legend, nrow = 2, heights = c(6.5,1.5))
figS1tot

dev.copy2pdf(file = "prophylactic_scenarios/figures/FigureS1.pdf", useDingbats=FALSE, 
              width = 8, height = 5)


