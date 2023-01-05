###################################################################
## This script generates Figure S2 
## This scenario shows vaccine allocation (incl optimal one)
## across two populations of *unequal* size where one 
## population is 10x the size of the other.
## The inits are more detailed here to allow for v > N.
###################################################################

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
# Figure 2 example is with two populations of size 1e6 and 2e6.
N1 <- 1e6 ## population 1 size
N2 <- 10e6 ## population 2 10 times size of population 1
ve <- 0.95 ## vaccine efficacy
tau <- 1-ve ## failure rate 
phi1 <- 0 ## underlying immunity population 1
phi2 <- 0 ## underlying immunity population 2 
years <- 3

dt <- seq(0, 365, 1)

R01 = 2 ## R0 is beta*N/gamma 
R02 = 2 ## R0 is beta*N/gamma

## use params from kissler et al. 
params <- c(sigma = 1/3, ## 1/latent period 
            nu = 1/5, ## 1/infectious period
            beta1 = R01/(5*N1), ## per contact probability of transmission
            beta2 = R02/(5*N2)
            )

##########################################
## Function for two populations 
##########################################

SEIR.twopop.nb <- function(t, x, params){
  
  with(as.list(c(params,x)),{
    
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

unequal10 <- data.frame()

v_list <- seq(0, 10e6, 2500) # seq(0, 2e6, 2500)
pv_list <- seq(0, 1, 0.01)
pct <- N1 * (1/1000)
pct2 <- N2 * (1/1000)

for(v in v_list){
  for(pv1 in pv_list){
    #for(phi1 in phi_list){
    
    #Following the run
    print(paste("v =",v, 
                "pv =",pv1, 
                sep =""))
    
    #Need to account for I in each population
    N1 = N1 - pct
    N2 = N2 - pct2
    
    # Only allow v to be up to total population size (N1+N2)
    # we can't vaccinate more than the total population size (in both pops)
    if(v > (N1+N2)){
      waste <- v - (N1+N2)
      v <- N1+N2 ## cap v to be the max number of doses. 
    }
    
    # If number of doses for one pop > total size that pop, 
    # we give the leftover doses to the other population.
    
    # First, calculate if there will be leftover and how much 
    if((v * pv1) > N1){
      left_for_v2 = (v*pv1) - N1
    } else{
      left_for_v2 = 0 ## added to sum later 
    }
    
    if((v * (1-pv1)) > N2){
      left_for_v1 = (v*(1-pv1)) - N2
    } else{
      left_for_v1 = 0 ## added to sum later 
    }
    
    # Then assign S, V and R accounting for the leftover doses that might exist
    if((v * pv1) > N1){
      S1 = (tau)*N1
      V1 = (1-tau)*N1
      R1 = 0
    } else {
      S1 = N1 - left_for_v1*(1-tau) - (v*pv1)*(1-tau)
      V1 = left_for_v1*(1-tau) + (v*pv1)*(1-tau)
      R1 = 0
    }
    
    if((v * (1-pv1)) > N2) {
      S2 = (tau)*N2
      V2 = (1-tau)*N2
      R2 = 0
    } else {
      S2 = N2 - left_for_v2*(1-tau) - (v*(1-pv1))*(1-tau)
      V2 = left_for_v2*(1-tau) + (v*(1-pv1))*(1-tau)
      R2 = 0
    }
    
    # Set E, I and C as they're always the same 
    # (doing it here and not in the if loops works much better)
    E1 = 0
    I1 = pct
    C1 = pct
    E2 = 0
    I2 = pct2
    C2 = pct2
    
    # Bring back N1 and N2 to their actual total size now that we have some individual infected (I)
    N1 <- N1+pct
    N2 <- N2+pct2
    
    # Set the initial conditions
    inits <- c(
      S1 = S1,
      V1 = V1,
      E1 = E1, 
      I1 = I1, 
      R1 = R1, 
      C1 = C1,
      S2 = S2, 
      V2 = V2, 
      E2 = E2, 
      I2 = I2, 
      R2 = R2, 
      C2 = C2
    )
    
    output2 <- lsoda(y = inits,
                     times = dt,
                     func = SEIR.twopop.nb,
                     parms = params)
    out2.df <- data.frame(output2)[length(dt),]
#    out2.df <- data.frame(output2)[366,]
    out2.df$pv = pv1
    out2.df$v = v
    out2.df$phi1 = phi1
    out2.df$phi2 = phi2
    out2.df$left_for_v1 = left_for_v1
    out2.df$left_for_v2 = left_for_v2
    unequal10 = rbind(unequal10,out2.df)
    #}
  }
}


######################################
## Creating Figure 2 plot pop 2 2 times as large
######################################

### Check df

twopop2.2 %>%
  mutate(vtau = v*(1-tau),
         vtot = V1 + V2) -> twopop2.2

#save(twopop2.2, file = "unequalallocation_20200528.RData")
#save(twopop2.2, file = "unequalallocation_pop1-1_pop2-2_20221116.RData")
#save(twopop2.2, file = "unequalallocation_pop1-1_pop2-10_20221116.RData")

### Load data
#load("unequalallocation_20200528.RData")
load("unequalallocation_pop1-1_pop2-2_20221116.RData")
#load("unequalallocation_pop1-1_pop2-10_20221116.RData")

### Subset data 
twopop <- twopop2.2 %>% 
  filter(time == max(time)) %>% 
  mutate(C1C2 = C1 + C2) %>%
  select(C1C2, pv, v, C1, C2) # %>% 
#  filter(phi1 == 0)

data_top <- twopop %>%
  filter(v < 1.8e6) %>%
  mutate(pop1_doses = pv * v) %>%
  mutate(pop2_doses = v - pop1_doses) %>%
  group_by(v) %>%
  slice(which.min(C1C2))

data_bottom <- twopop %>%
  filter(v %in% c(200000, 600000, 800000, 1000000,  1400000)) #c(2e5, 6e5, 8e5, 10e5, 14e5))
# filter(v %in% c(1e5, 3e5, 5e5, 7e5, 9e5 ))
# filter(v == 2100| v == 3500 | v == 4160 | v == 4500 | v == 6000)

### Labs
labs <- c("Regime 1\n(200000 doses)", "Regime 2\n(600000 doses)", "Regime 3\n(800000 doses)",
          "Regime 4\n(1000000 doses)", "Regime 5\n(1400000 doses)")
#names(labs) <- c("2100", "3500", "4160", "4500", "6000")
names(labs) <- c(200000, 600000, 800000, 1000000,  1400000)
labs

### Define my colors
my_col <- c("#B2182B", "grey65", "#4393C3")

### Top half of Figure 2
options(scipen = 99999)
plot_top <- ggplot(data_top, mapping = aes(x = v)) +
  geom_line(aes(y = pop1_doses, color = "Population 1 (size 1000000)"), size = 1.3) +
  geom_line(aes(y = pop2_doses, color = "Population 2 (size 2000000)"), size = 1.3) +
  geom_vline(xintercept= c(455000, 777500, 887500, 1225000), linetype="dotted", colour = "gray32") +
  geom_text(x=227500, y=1100000, label="1", color="gray32", size = 4) +
 # geom_text(x=227500, y=900000, label="pop 1 < pop 2", color="gray32", size = 4) +
  geom_text(x=616250, y=1100000, label="2", color="gray32", size = 4) +
  geom_text(x=832500, y=1100000, label="3", color="gray32", size = 4) +
  geom_text(x=1056250, y=1100000, label="4", color="gray32", size = 4) +
  geom_text(x=1500000, y=1100000, label="5", color="gray32", size = 4) +
  labs(#title = "Remaking Keeling/Duijzer plot",
    x = "Total number of vaccine doses",
    y = "Optimal vaccine allocation across \n populations 1 and 2") +
  theme_bw() + 
  theme(legend.title = element_blank(), 
        legend.position = 'bottom',
        #legend.spacing.x = unit(0.3, "cm"), 
       # panel.grid.major = element_blank(), 
      #  panel.grid.minor = element_blank(),
       # axis.text.x = element_blank(),
       # axis.ticks.x = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)) +
      #  legend.title = element_text(size = 10)) +
  coord_cartesian(ylim = c(0,1.2e6), clip = "off") +
  scale_y_continuous(breaks = seq(0, 1.2e6, 2e5), labels = seq(0, 1.2e6, 2e5)) +
# scale_x_continuous(breaks = seq(0, 2e6, 4e5), labels = seq(0, 2e6, 4e5)) +
  scale_x_continuous(breaks = seq(0, 2e6, 2e5), labels = seq(0, 2e6, 2e5)) +
  scale_color_manual(values = c("Population 1 (size 1000000)" = my_col[3], "Population 2 (size 2000000)" = my_col[1])) +
  geom_text(x = 250000, y=43500, label="population 2", size = 3, color="#B2182B") +
  geom_text(x = 220000, y=270000, label="population 1", size = 3, color="#4393C3", angle = 31)
plot_top

plot_top_legend <- get_legend(plot_top)

### Create df for plotting the minimal points
min_point <- data_bottom %>%
  group_by(v) %>%
  slice(which.min(C1C2)) %>%
  select(pv, v, C1C2) 
min_point <- as.data.frame(min_point)
min_point

### Bottom half of Figure 2
N1 <- 1e6
N2 <- 2e6
ratio <- N1/(N1+N2)
plot_bottom <- ggplot(data = data_bottom, aes(x = pv, y = C1C2, color = pv)) +
  geom_line(size = 2) +
  facet_wrap(~v, ncol = 5, labeller = labeller(v = labs)) + 
  theme_bw() +
  scale_color_gradientn(colours = my_col) + 
  labs(#title = "Title",
    x = "Proportion of vaccine doses given to population 1",
    y = "Cumulative number of cases in \n populations 1 and 2") +
  theme(legend.position = "none") +
  #scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
#  scale_x_continuous(labels = c(0, 0.2, 0.4, 0.8, 1)) +
  scale_x_continuous(labels = c(0, 0.25, 0.5, 0.75, 1)) + # c(0, 0.2, 0.4, 0.8, 1))
  theme(#panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()) +
  #geom_point(data = min_point, aes(x = pv[1], y = C1C2[1]), color = "black", inherit.aes=FALSE)
  geom_point(data = data.frame(x = min_point$pv[1], y = min_point$C1C2[1], v = 200000), aes(x=x, y=y), inherit.aes=FALSE, size = 1.5) +
  geom_point(data = data.frame(x = min_point$pv[2], y = min_point$C1C2[2], v = 600000), aes(x=x, y=y), inherit.aes=FALSE, size = 1.5) +
  geom_point(data = data.frame(x = min_point$pv[3], y = min_point$C1C2[3], v = 800000), aes(x=x, y=y), inherit.aes=FALSE, size = 1.5) +
  geom_point(data = data.frame(x = min_point$pv[4], y = min_point$C1C2[4], v = 1000000), aes(x=x, y=y), inherit.aes=FALSE, size = 1.5) +
  geom_point(data = data.frame(x = min_point$pv[5], y = min_point$C1C2[5], v = 1400000), aes(x=x, y=y), inherit.aes=FALSE, size = 1.5) +
  geom_vline(xintercept = ratio, linetype = "dashed", color = "grey")
plot_bottom

### Combine two plots
ggarrange(plot_top + rremove("legend"), 
          plot_bottom, 
          plot_top_legend, nrow = 3, heights = c(5,5,0.5))
dev.copy2pdf(file = "Fig2_unequal_1-2_20221202.pdf", useDingbats=FALSE,
            width = 8, height = 7)

#
###


######################################
## Creating Figure 2 plot pop 2 10 times as large
######################################

### Check df
twopop2.2 %>%
  mutate(vtau = v*(1-tau),
         vtot = V1 + V2) -> twopop2.2

#save(twopop2.2, file = "unequalallocation_20200528.RData")
#save(twopop2.2, file = "unequalallocation_pop1-1_pop2-2_20221116.RData")
#save(twopop2.2, file = "unequalallocation_pop1-1_pop2-10_20221116.RData")

### Load data
#load("unequalallocation_20200528.RData")
#load("unequalallocation_pop1-1_pop2-2_20221116.RData")
load("unequalallocation_pop1-1_pop2-10_20221116.RData")

N1 <- 1e6
N2 <- 10e6
v_list <- seq(0, 10e6, 2500)

### Subset data 
twopop <- twopop2.2 %>% 
  filter(time == max(time)) %>% 
  mutate(C1C2 = C1 + C2) %>%
  select(C1C2, pv, v, C1, C2) # %>% 
#  filter(phi1 == 0)

data_top <- twopop %>%
  filter(v < (N2*0.6)) %>%
  mutate(pop1_doses = pv * v) %>%
  mutate(pop2_doses = v - pop1_doses) %>%
  group_by(v) %>%
  slice(which.min(C1C2))

data_bottom <- twopop %>%
  filter(v %in% c(400000, 2000000, 3600000, 4300000, 5000000))

### Labs
labs <- c("Regime 1\n(400000 doses)", "Regime 2\n(2000000 doses)", "Regime 3\n(3600000 doses)",
          "Regime 4\n(4300000 doses)", "Regime 5\n(5000000 doses)")
#names(labs) <- c("2100", "3500", "4160", "4500", "6000")
names(labs) <- c(400000, 2000000, 3600000, 4300000, 5000000)
labs

### Define my colors
my_col <- c("#B2182B", "grey65", "#4393C3")

### Top half of Figure 2
options(scipen = 99999)
plot_top <- ggplot(data_top, mapping = aes(x = v)) +
  geom_line(aes(y = pop1_doses, color = "Population 1 (size 1,000,000)"), size = 1.3) +
  geom_line(aes(y = pop2_doses, color = "Population 2 (size 10,000,000)"), size = 1.3) +
  geom_vline(xintercept= c(455000, 3252500, 4200000, 4412500), linetype="dotted", colour = "gray32") +
   geom_text(x=150000, y=5000000, label="1", color="gray32", size = 4) +
   geom_text(x=1800000, y=5000000, label="2", color="gray32", size = 4) +
   geom_text(x=3750000, y=5000000, label="3", color="gray32", size = 4) +
   geom_text(x=4300000, y=5000000, label="4", color="gray32", size = 4) +
   geom_text(x=5100000, y=5000000, label="5", color="gray32", size = 4) +
  labs(#title = "Remaking Keeling/Duijzer plot",
    x = "Total number of vaccine doses",
    y = "Optimal vaccine allocation across \n populations 1 and 2") +
  theme_bw() + 
  theme(legend.title = element_blank(), 
        legend.position = 'bottom',
        #legend.spacing.x = unit(0.3, "cm"), 
        # panel.grid.major = element_blank(), 
        #  panel.grid.minor = element_blank(),
        # axis.text.x = element_blank(),
        # axis.ticks.x = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)) +
  #  legend.title = element_text(size = 10)) +
  #coord_cartesian(ylim = c(0,1.2e6), clip = "off") +
  scale_y_continuous(breaks = seq(0, 10e6, 10e5), labels = seq(0, 10e6, 10e5)) +
  scale_x_continuous(breaks = seq(0, 10e6, 8e5), labels = seq(0, 10e6, 8e5)) +
  scale_color_manual(values = c("Population 1 (size 1,000,000)" = my_col[3], "Population 2 (size 10,000,000)" = my_col[1])) +
  geom_text(x = 250000, y=43500, label="population 2", size = 3, color="#B2182B", angle = 31) +
  geom_text(x = 220000, y=270000, label="population 1", size = 3, color="#4393C3")
plot_top

plot_top_legend <- get_legend(plot_top)

### Create df for plotting the minimal points
min_point <- data_bottom %>%
  group_by(v) %>%
  slice(which.min(C1C2)) %>%
  select(pv, v, C1C2) 
min_point <- as.data.frame(min_point)
min_point

### Bottom half of Figure 2
ratio <- N1/(N1+N2)
plot_bottom <- ggplot(data = data_bottom, aes(x = pv, y = C1C2, color = pv)) +
  geom_line(size = 2) +
  facet_wrap(~v, ncol = 5, labeller = labeller(v = labs)) + 
  theme_bw() +
  scale_color_gradientn(colours = my_col) + 
  labs(#title = "Title",
    x = "Proportion of vaccine doses given to population 1",
    y = "Cumulative number of cases in \n populations 1 and 2") +
  theme(legend.position = "none") +
  #scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
  #  scale_x_continuous(labels = c(0, 0.2, 0.4, 0.8, 1)) +
  scale_x_continuous(labels = c(0, 0.25, 0.5, 0.75, 1)) + # c(0, 0.2, 0.4, 0.8, 1))
  theme(#panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()) +
  #geom_point(data = min_point, aes(x = pv[1], y = C1C2[1]), color = "black", inherit.aes=FALSE)
  geom_point(data = data.frame(x = min_point$pv[1], y = min_point$C1C2[1], v = 400000), aes(x=x, y=y), inherit.aes=FALSE, size = 1.5) +
  geom_point(data = data.frame(x = min_point$pv[2], y = min_point$C1C2[2], v = 2000000), aes(x=x, y=y), inherit.aes=FALSE, size = 1.5) +
  geom_point(data = data.frame(x = min_point$pv[3], y = min_point$C1C2[3], v = 3600000), aes(x=x, y=y), inherit.aes=FALSE, size = 1.5) +
  geom_point(data = data.frame(x = min_point$pv[4], y = min_point$C1C2[4], v = 4300000), aes(x=x, y=y), inherit.aes=FALSE, size = 1.5) +
  geom_point(data = data.frame(x = min_point$pv[5], y = min_point$C1C2[5], v = 5000000), aes(x=x, y=y), inherit.aes=FALSE, size = 1.5) +
  geom_vline(xintercept = ratio, linetype = "dashed", color = "grey")
plot_bottom

### Combine two plots
ggarrange(plot_top + rremove("legend"), 
          plot_bottom, 
          plot_top_legend, nrow = 3, heights = c(5,5,0.5))
dev.copy2pdf(file = "Fig2_unequal_1-10_20221202.pdf", useDingbats=FALSE,
             width = 8, height = 7)

#