rm(list=ls())
source("01_functions_cluster.R")
library(statmod)
require(tidyr)
library(cowplot)
library(dplyr)
library(readr)
library(beepr)
library(viridis)
library(ggdist)
library(dplyr)
library(ggplot2)
library(beepr)
library(viridis)





mydir = 'datasets_1'
myfiles = list.files(path=mydir, pattern="*.csv", full.names=TRUE)
#myfiles<-myfile
newfiles<-myfiles[1:154]
#load processed network data from the cluster
load("Network_all_data.RData")


#plot of temperature vs biomass figure 2A
(r0<-net_dat %>% ggplot(aes(x= Temperature, y =biomass, color=factor(h2)))+
    geom_point(size =2.25, alpha =0.25)+
    theme_bw()+
    scale_color_brewer(palette="Dark2")+
    facet_wrap(.~individual_variation)+
    ylab("Equilibrium biomass")+
    labs(color="heritability")+
    stat_smooth(geom="smooth", 
                method = "loess", 
                se=F , size =2.5, alpha =0.25))


## appendix figure S4 connectance and network size vs biomass
(app1<-net_dat %>% filter(Temperature < 30) %>% 
    ggplot(aes(x= network_size, y =biomass, color=factor(h2)))+
    geom_point(size =2.25, alpha =0.25)+
    theme_bw()+
    xlab("Network size")+
    scale_color_brewer(palette="Dark2")+
    facet_wrap(.~individual_variation)+
    ylab("Equilibrium biomass")+
    labs(color="heritability")+
    stat_smooth(geom="smooth", 
                method = "lm", 
                se=F , size =2.5, alpha =0.25))
(app2<-net_dat %>% filter(Temperature < 30) %>% 
    ggplot(aes(x= connectance, y =biomass, color=factor(h2)))+
    geom_point(size =2.25, alpha =0.25)+
    theme_bw()+
    xlab("Connectance")+
    scale_color_brewer(palette="Dark2")+
    facet_wrap(.~individual_variation)+
    ylab("Equilibrium biomass")+
    labs(color="heritability")+
    stat_smooth(geom="smooth", 
                method = "lm", 
                se=F , size =2.5, alpha =0.25))
## appendix figure S4 connectance and network size vs biomass
ggpubr::ggarrange(app1,app2, ncol =1,nrow = 2, labels = c("A","B"))



#temperature at collapse analysis 
webfiles<-unique(net_dat$web)

fact_1<- expand.grid(h2=c(0,0.2, 0.4),
                     individual_variation = c("high variation","low variation"),
                     webfiles = webfiles) %>% 
  mutate(temperature_collapse_richness=0,
         temperature_collapse_biomass=0,
         abrupt_collapse=0,
         nestedness=0,
         connectance=0,
         network_size=0,
         trait_lag=0,
         trait_matching=0)
abruptness<-numeric()
for(i in 1: nrow(fact_1)){
  
  temp_high<-net_dat %>% filter(individual_variation == fact_1$individual_variation[i],
                                web == fact_1$webfiles[i], h2 ==fact_1$h2[i] )
  
  temp_at_collapse_richness<-min(temp_high$Temperature[temp_high$richness <= 1])#min(temp_high$Temperature[temp_high$richness <= 0])
  temp_at_collapse_biomass<-min(temp_high$Temperature[temp_high$biomass < 0.1])
  
  index<-which(temp_high$Temperature == min(temp_high$Temperature[temp_high$richness <= 0]))
  fact_1$temperature_collapse_richness[i] <-temp_at_collapse_richness 
  fact_1$temperature_collapse_biomass[i]<-temp_at_collapse_biomass
  
  change_richness<- abs(temp_high$richness[index] - 
                          temp_high$richness[(index-2):index])
  
  #which((-diff(temp_high$richness)) > round(0.25*temp_high$network_size[1]))
  
  d<-min(which(abs(diff(temp_high$trait_lag)) > 5) )
  
  if((max(tail(change_richness)) >= 0.1*temp_high$network_size[1])){
    abruptness <- 1
  } else{abruptness <- 0}
  
  fact_1$abrupt_collapse[i]<- abruptness
  fact_1$nestedness[i]<-temp_high$nestedness[1]
  fact_1$connectance[i] <- temp_high$connectance[1]
  fact_1$network_size[i] <- temp_high$network_size[1]
 # fact_1$skewness_lag[i] <- skewness(temp_high$trait_lag)
  fact_1$trait_lag[i]<-  (temp_high$trait_lag[d])
  
  fact_1$trait_matching[i]<-  max(temp_high$trait_matching[1:index])
  print(i)
  
}

#figure 2B
(r1<-fact_1 %>% ggplot(aes(x =factor(h2), y = temperature_collapse_biomass, 
                           color=individual_variation))+
    geom_boxplot(outlier.colour  =NA)+
    theme_bw()+
    xlab("Heritability")+
    ylab("Temperature at collapse")+
    ylim(c(25,38))+
    scale_color_brewer(palette="Set1")+
    labs(color="Trait variation")+
    geom_jitter(aes(x =factor(h2), y = temperature_collapse_biomass,color=individual_variation),
                position = position_dodge2(0.8), size =3, alpha =0.3))

#grouping figure 2 A, 2B
ggpubr::ggarrange(r0, r1, labels = c("A","B"), nrow = 1, ncol = 2)


#fact_1$delta[fact_1$delta == 0.05] = "threshold 5~%"
#fact_1$delta[fact_1$delta == 0.15] = "threshold 15~%"
#some aggregated stats on collapse and trait variatoin and h2
fact_1 %>% 
  group_by( factor(h2), individual_variation ) %>% 
  summarise(mean_collapse = mean(temperature_collapse_biomass,na.rm = T),
            sd = sd(temperature_collapse_biomass,na.rm=T))


# Appendix figure 1 of trait matching
(rx<-fact_1 %>% ggplot(aes(x =factor(h2), y = trait_matching, 
                           color=individual_variation))+
    geom_boxplot(outlier.colour  =NA)+
    theme_bw()+
    xlab("Heritability")+
    ylab("Trait matchiing")+
    scale_color_brewer(palette="Set1")+
    labs(color="Trait variation")+
    geom_jitter(aes(x =factor(h2), y = trait_matching,color=individual_variation),
                position = position_dodge2(0.8), size =3, alpha =0.3))



#temperature at collapse for network size and connectance and abrupt collapses i.e. figure 3A-D

(r3<-fact_1 %>% ggplot(aes(y=temperature_collapse_biomass, x = network_size, color = factor(h2)))+
    geom_point(size=5)+
    theme_bw()+
    scale_color_brewer(palette="Dark2")+
    xlab("Network size")+
    ggtitle("A")+
    labs(col="heritability")+
    ylim(c(27,36))+
    ylab("Temperature at collapse")+
    geom_smooth(method = "lm", formula = y~x, se=T)+
    facet_wrap(.~individual_variation))  

(r4<-fact_1 %>% ggplot(aes(y=temperature_collapse_biomass, x = connectance, color = factor(h2)))+
    geom_point(size=5)+
    theme_bw()+
    scale_color_brewer(palette="Dark2")+
    xlab("Connectance")+
    labs(col="heritability")+
    ylim(c(27,36))+
    ggtitle("B")+
    ylab("Temperature at collapse")+
    geom_smooth(method = "lm", formula = y~x, se=T)+
    facet_wrap(.~individual_variation)  )


(r5<-fact_1 %>% ggplot(aes(x=network_size, y = abrupt_collapse, color =factor(h2)))+
    geom_point(size= 5, alpha=0.5)+
    geom_smooth(method = "glm", 
                method.args = list(family = "quasibinomial"), 
                se = T, size=2.5 , alpha=0.15) +
    scale_color_brewer(palette="Dark2")+
    ylab("Abrupt collapse")+
    xlab("Network size")+
    ggtitle("C")+
    labs(col="heritability")+
    theme_bw()+
    facet_wrap(.~individual_variation))

(r6<-fact_1 %>% ggplot(aes(x=connectance, y = abrupt_collapse, color =factor(h2)))+
    geom_point(size= 5, alpha=0.5)+
    geom_smooth(method = "glm", 
                method.args = list(family = "quasibinomial"), 
                se = T, size=2.5 , alpha=0.15) +
    scale_color_brewer(palette="Dark2")+
    #ggtitle("D")+
    xlab("Connectance")+
    ylab("Abrupt collapse")+
    labs(col="heritability")+
    theme_bw()+
    facet_grid(delta~individual_variation))

library(grid)
library(gridExtra)
#grouping the figures for 3A,B, 3C, 3D
grid_arrange_shared_legend(r3,r4,r5,r6,ncol = 2,nrow = 2)




#Figure S2 --- in the appendix- NODF nestedness vs collapse threshold and aburpt collapse
(r7<-fact_1 %>% ggplot(aes(y=temperature_collapse_biomass, x = nestedness, color = factor(h2)))+
    geom_point(size=5)+
    theme_bw()+
    scale_color_brewer(palette="Dark2")+
    xlab("Nestedness(NODF)")+
    ggtitle("A")+
    labs(col="heritability")+
    ylim(c(25,36))+
    ylab("Temperature at collapse")+
    geom_smooth(method = "lm", formula = y~x, se=T)+
    facet_wrap(.~individual_variation))  


(r8<-fact_1 %>% ggplot(aes(x=nestedness, y = abrupt_collapse, color =factor(h2)))+
    geom_point(size= 5, alpha=0.5)+
    geom_smooth(method = "glm", 
                method.args = list(family = "quasibinomial"), 
                se = T, size=2.5 , alpha=0.15) +
    scale_color_brewer(palette="Dark2")+
     ggtitle("B")+
    xlab("Nestedness(NODF)")+
    ylab("Abrupt collapse")+
    labs(col="heritability")+
    theme_bw()+
    facet_wrap(.~individual_variation))
#Figure S2 grouping--- in the appendix- NODF nestedness vs collapse threshold and aburpt collapse
ggpubr::ggarrange(r7,r8)


# figure 4B- trait lag and temperature 
(n1<-net_dat %>% ggplot(aes(x= Temperature, y =trait_lag, color=factor(h2)))+ 
    geom_point(size =3, alpha =0.5)+
    theme_bw()+
    ylim(c(-4,20))+
    ggtitle("B")+
    scale_color_brewer(palette="Dark2")+
    facet_wrap(.~individual_variation)+
    ylab("network trait lag")+
    labs(color="heritability"))



