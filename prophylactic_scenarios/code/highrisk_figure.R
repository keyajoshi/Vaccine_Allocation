####################################
## This script generates Figure 5 
## and Figure S8 using the data 
## generated in highrisk_mortality.R
## and highrisk_transmission.R
#######################

###################################
## Import libraries
####################################
library(deSolve)
library(grid)
library(gridExtra)
library(viridis)
library(tidyverse)
options(scipen = 999)

#################################
## Load & manipulate data for 
## figures 
#################################
v_list <- seq(0, 1e6, 1e5)
pv_list <- seq(0,1, 0.1)
ph_list <- c(0.25, 0.5)
r_list <- c(2,4,8,16)

#####
## high risk mortality
#####

load("prophylactic_scenarios/data/highrisk_mortality.Rdata")

highrisk_mortality %>% 
  filter(time == max(time)) %>% 
  mutate(D1D2 = (D1H + D2H + D1L + D2L), 
         D1D2H = (D1H + D2H), 
         D1D2L = (D1L + D2L)) -> highrisk_mortality.sub

highrisk_mortality.sub%>%
  count(RHH1)

max_list_hd <- length(pv_list)*
  length(v_list) * 
  length(ph_list)

highrisk_mortality.sub$xaxis <- rep(1:max_list_hd, 
                                    each = length(r_list))

#####
## high risk transmission 
#####

load("prophylactic_scenarios/data/highrisk_transmission.Rdata")

highrisk_transmission %>% 
  data.frame() %>%
  dplyr::filter(time == max(time)) %>%
  mutate(C1C2 = C1H + C2H + C1L + C2L, 
         C1C2H = C1H + C2H, 
         C1C2L = C1L + C2L) -> highrisk_transmission.sub

highrisk_transmission.sub%>%
  count(R)

max_list_hr <- length(pv_list)*
  length(v_list) * 
  length(ph_list)

highrisk_transmission.sub$xaxis <- rep(1:max_list_hr, 
                     each = length(r_list))

##############
## plot options
##############
my.col.update <- c("#FF410DFF", "#6EE2FFFF", "#5B1A18", "#748AA6FF", "#F79D1EFF", "#D0DFE6FF",  "#02401B", "#FD6467", "#C27D38", "#C6CDF7" )

#############
## figure 5
## 25% of population high risk 
#############

hr_death <- highrisk_mortality.sub %>% 
  filter(v !=0,
         ph1 == 0.25) %>%
  group_by(v, pv, D1D2, RHH1, xaxis) %>% 
  summarize() %>%
  ggplot(aes(x = xaxis , 
             y = D1D2, 
             col = as.factor(v),
             linetype = as.factor(RHH1))) + 
  geom_line() + 
  theme_light() + 
  labs(x = "Combination of total number of vaccine doses & proportion given to population 1", 
       y = "Cumulative number of deaths \n in population 1 & 2", 
       color = "Vaccine Doses", 
       linetype = "Proportion high risk \n in population 1") + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none")  + 
  ggtitle("25% High Mortality in both Populations") +
  scale_color_manual(values = my.col.update) +
  scale_y_continuous(limits = c(0,4000), 
                     labels = seq(0,4000, 1000), 
                     breaks = seq(0,4000,1000)) + 
  guides(linetype = "none")

hr_transmission <-highrisk_transmission.sub %>% 
  filter(v !=0,
         ph1 == 0.25) %>%
  group_by(v, pv, C1C2, R, xaxis) %>% 
  summarize() %>%
  ggplot(aes(x = xaxis, 
             y = C1C2, col = as.factor(v))) + 
  geom_line(aes(linetype = as.factor(R))) +
  theme_light() +
  labs(x = "Combination of total number of vaccine doses & proportion given to population 1", 
       y = "Cumulative number of cases \n in population 1 & 2", 
       color = "Vaccine Doses",
       linetype = "R0 value") + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "bottom",
        legend.box = "vertical",
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold")) + 
  ggtitle("25% High Transmitters in both Populations") + 
  coord_cartesian(ylim = c(0, 20e5), clip = "off") + 
  scale_y_continuous(breaks= seq(0,20e5,2e5), labels = seq(0,20e5,2e5)) + 
  scale_color_manual(values = my.col.update) + 
  guides(linetype = guide_legend(order = 2),col = guide_legend(order = 1))

ggarrange(hr_transmission, hr_death, 
          nrow = 2, 
          common.legend = TRUE, 
          legend = "bottom")

dev.copy2pdf(file = "prophylactic_scenarios/figures/Figure5.pdf", useDingbats=FALSE, 
             width = 11, height = 8)

#############
## figure S8
## 50% of population high risk 
#############

hr_death_50 <- highrisk_mortality.sub %>% 
  filter(v !=0,
         ph1 == 0.5) %>%
  group_by(v, pv, D1D2, RHH1, xaxis) %>% 
  summarize() %>%
  ggplot(aes(x = xaxis , 
             y = D1D2, 
             col = as.factor(v),
             linetype = as.factor(RHH1))) + 
  geom_line() + 
  theme_light() + 
  labs(x = "Combination of total number of vaccine doses & proportion given to population 1", 
       y = "Cumulative number of deaths \n in population 1 & 2", 
       color = "Vaccine Doses", 
       linetype = "Proportion high risk \n in population 1") + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none")  + 
  ggtitle("50% High Mortality in both Populations") +
  scale_color_manual(values = my.col.update) +
  scale_y_continuous(limits = c(0,6000), 
                     labels = seq(0,6000, 1000), 
                     breaks = seq(0,6000,1000)) + 
  guides(linetype = "none")

hr_transmission_50 <-highrisk_transmission.sub %>% 
  filter(v !=0,
         ph1 == 0.5) %>%
  group_by(v, pv, C1C2, R, xaxis) %>% 
  summarize() %>%
  ggplot(aes(x = xaxis, 
             y = C1C2, col = as.factor(v))) + 
  geom_line(aes(linetype = as.factor(R))) +
  theme_light() +
  labs(x = "Combination of total number of vaccine doses & proportion given to population 1", 
       y = "Cumulative number of cases \n in population 1 & 2", 
       color = "Vaccine Doses",
       linetype = "R0 value") + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "bottom",
        legend.box = "vertical",
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold")) + 
  ggtitle("50% High Transmitters in both Populations") + 
  coord_cartesian(ylim = c(0, 20e5), clip = "off") + 
  scale_y_continuous(breaks= seq(0,20e5,2e5), labels = seq(0,20e5,2e5)) + 
  scale_color_manual(values = my.col.update) + 
  guides(linetype = guide_legend(order = 2),col = guide_legend(order = 1))

ggarrange(hr_transmission_50, hr_death_50, 
          nrow = 2, 
          common.legend = TRUE, 
          legend = "bottom")

dev.copy2pdf(file = "prophylactic_scenarios/figures/FigureS8.pdf", useDingbats=FALSE, 
             width = 11, height = 8)
