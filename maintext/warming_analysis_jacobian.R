rm(list=ls())
source("01_functions_cluster.R")
library(statmod)
#require(deSolve) ## for integrating ordinary differential equations
require(tidyr) ## for efficient data manipulation & plotting
library(cowplot) ## for arranging plots in a grid
library(dplyr)
library(readr)
library(beepr)
library(viridis)
library(ggdist)
library(dplyr)
#library(readr)
library(beepr)
library(viridis)
library(moments)
library(ggplot2)
#theme_set(theme_classic()) 
load("Network_data.RData")

mydir = 'datasets_1'
myfiles = list.files(path=mydir, pattern="*.csv", full.names=TRUE)
#myfiles<-myfile
newfiles<-myfiles[1:154]

# load("cluster_warming_1.RData")
# 
#  net_dat<-sp_dat<-NULL
#  for(i in 1:10720){
#    print(i)
#  net_dat<-rbind(net_dat, output[[i]]$output_network_data)  
#  sp_dat<-rbind(sp_dat,output[[i]]$species_level_data)  
#  }
#   
# # 
 # load("cluster_warming_3.RData")
# # 
net_dat<-sp_dat1<-NULL
  for(i in 1:9348){
    print(i)
    net_dat<-rbind(net_dat, output[[i]]$output_network_data)  
   sp_dat1<-rbind(sp_dat1,output[[i]]$species_level_data)  
  }
# 
# str(net_dat)
# 
# net_dat<-rbind(net_dat,net_dat1)
# 
# spdat<-rbind(sp_dat,sp_dat1)
# save(spdat, file ="Species_data.RData")
# 
# save(net_dat, file ="Network_data.RData")
#temperature vs biomass and richness

(r0<-net_dat %>% ggplot(aes(x= Temperature, y =biomass, color=factor(h2)))+
  geom_point(size =4, alpha =0.1)+
  theme_bw()+
  scale_color_brewer(palette="Dark2")+
  facet_wrap(.~individual_variation)+
  ylab("Equilibrium biomass")+
  labs(color="heritability")+
  stat_smooth(geom="smooth", method = "loess", se=F , size =3, alpha =0.25))
                  

#temperature vs. richness 

net_dat %>% ggplot(aes(x= Temperature, y =trait_matching, color=factor(h2)))+
  geom_point(size =4, alpha =0.1)+
  theme_bw()+
  facet_wrap(.~individual_variation)+
  ylab("Richness")+
  labs(color="heritability")

net_dat$individual_variation<- plyr::revalue(net_dat$individual_variation ,c("high"= "high variation",
                                                                             "low" = "low variation"))
#temperature at collapse analysis 
webfiles<-unique(net_dat$web)

fact_1<- expand.grid(h2=c(0,0.4),
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
  temp_at_collapse_biomass<-min(temp_high$Temperature[temp_high$biomass <= 0.1])
  
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
  fact_1$skewness_lag[i] <- skewness(temp_high$trait_lag)
  fact_1$trait_lag[i]<-  (temp_high$trait_lag[d])
  
  fact_1$trait_matching[i]<-  (temp_high$trait_matching[temp_at_collapse_richness])
  print(i)
  
}

fact_1 %>% 
  group_by( factor(h2), individual_variation ) %>% 
  summarise(mean_collapse = mean(temperature_collapse_biomass,na.rm = T),
            sd = sd(temperature_collapse_biomass,na.rm=T))
#TEMPERATURE AT temp_at_collapse_biomass
(r1<-fact_1 %>% ggplot(aes(x =factor(h2), y = temperature_collapse_biomass, color=individual_variation))+
  geom_boxplot(outlier.colour  =NA)+
  theme_bw()+
  xlab("Heritability")+
  ylab("Temperature at collapse")+
 ylim(c(25,38))+
  scale_color_brewer(palette="Set1")+
  labs(color="Trait variation")+
  geom_jitter(aes(x =factor(h2), y = temperature_collapse_biomass,color=individual_variation),
              position = position_dodge2(0.8), size =3, alpha =0.5))
  

ggpubr::ggarrange(r0, r1, labels = c("A","B"), nrow = 1, ncol = 2)

#TEMPERATURE AT COLLAPSE
fact_1 %>% ggplot(aes(x =factor(h2), y = temperature_collapse_richness,
                      color=individual_variation))+
  geom_boxplot(outlier.colour  =NA)+
  theme_bw()+
  ylim(c(25,38))+
  xlab("Heritability")+
  ylab("Temperature at collapse")+
  scale_color_brewer(palette="Dark2")+
  labs(color="Trait variation")+
  geom_jitter(aes(x =factor(h2), y = temperature_collapse_richness,color=individual_variation),
              position = position_dodge2(0.8), size =3, alpha =0.25)


#temperature at collapse and network size

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
    ggtitle("D")+
    xlab("Connectance")+
  ylab("Abrupt collapse")+
    labs(col="heritability")+
  theme_bw()+
  facet_wrap(.~individual_variation))

library(grid)
library(gridExtra)

grid_arrange_shared_legend(r3,r4,r5,r6,ncol = 2,nrow = 2)

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

grid_arrange_shared_legend(r7,r8, nrow = 1,ncol = 2)

(n1<-net_dat %>% ggplot(aes(x= Temperature, y =trait_lag, color=factor(h2)))+ 
  geom_point(size =4, alpha =0.21)+
  theme_bw()+
  ylim(c(-4,20))+
    ggtitle("B")+
  scale_color_brewer(palette="Dark2")+
  facet_wrap(.~individual_variation)+
  ylab("network trait lag")+
  labs(color="heritability"))



(n2<-fact_1 %>% ggplot(aes(y= (trait_lag), x = network_size, color = factor(h2)))+
  geom_point(size=5, alpha=0.8)+
  theme_bw()+
  scale_color_brewer(palette="Dark2")+
  xlab("connectance")+
  ylab("Network trait lag at collapse")+
    labs(color="heritability")+
  geom_smooth(method = "lm", formula = y~x, se=F)+
  facet_wrap(.~individual_variation)  )


ggpubr::ggarrange(n1,n2, nrow=2,ncol = 1,labels = c("A","B"))
fact_1 %>% ggplot(aes(y=trait_lag, x = nestedness, color = factor(h2)))+
  geom_point(size=5)+
  theme_bw()+
  scale_color_brewer(palette="Dark2")+
  xlab("nestedness (NODF)")+
  ylab("Network trait lag")+
  geom_smooth(method = "lm", formula = y~x, se=F)+
  facet_wrap(.~individual_variation)  


fact_1 %>% ggplot(aes(x =factor(h2), y = (trait_lag), color=individual_variation))+
  geom_boxplot(outlier.colour  =NA)+
  theme_bw()+
  xlab("Heritability")+
  ylab("Mean network trait-lag at collapse")+
  ylim(c(-2,5))+
  scale_color_brewer(palette="Dark2")+
  labs(color="Trait variation")+
  geom_jitter(aes(x =factor(h2), y = trait_lag,color=individual_variation),
              position = position_dodge2(0.8), size =3, alpha =0.25)




(r9<-fact_1 %>% ggplot(aes(y=trait_matching, x = network_size, color = factor(h2)))+
    geom_point(size=5)+
    theme_bw()+
    scale_color_brewer(palette="Dark2")+
    xlab("Nestedness(NODF)")+
    ggtitle("A")+
    labs(col="heritability")+
    ylab("Temperature at collapse")+
    geom_smooth(method = "lm", formula = y~x, se=F)+
    facet_wrap(.~individual_variation))  





#######
load("Species_data.RData")
str(sp_dat1)

sp_dat1$Degree<-as.numeric(as.character(sp_dat1$Degree))
sp_dat1$Species<-as.numeric(as.character(sp_dat1$Species))
sp_dat1$Nestedness<-as.numeric(as.character(sp_dat1$Nestedness))
sp_dat1$Connectance<-as.numeric(as.character(sp_dat1$Connectance))
sp_dat1$Individual_variation<-as.factor(as.character(sp_dat1$Individual_variation))
sp_dat1$Network_size<-as.numeric(as.character(sp_dat1$Network_size))
sp_dat1$mutualism_strength<-as.numeric(as.character(sp_dat1$mutualism_strength))
sp_dat1$density<-as.numeric(as.character(sp_dat1$density))
sp_dat1$trait_values<-as.numeric(as.character(sp_dat1$trait_values))
sp_dat1$Temperature<-as.numeric(as.character(sp_dat1$Temperature))
str(sp_dat1)


spdat$Degree<-as.numeric(as.character(spdat$Degree))
spdat$Species<-as.numeric(as.character(spdat$Species))
spdat$Nestedness<-as.numeric(as.character(spdat$Nestedness))
spdat$Connectance<-as.numeric(as.character(spdat$Connectance))
spdat$Individual_variation<-as.factor(as.character(spdat$Individual_variation))
spdat$Network_size<-as.numeric(as.character(spdat$Network_size))
spdat$mutualism_strength<-as.numeric(as.character(spdat$mutualism_strength))
spdat$density<-as.numeric(as.character(spdat$density))
spdat$trait_values<-as.numeric(as.character(spdat$trait_values))
spdat$Temperature<-as.numeric(as.character(spdat$Temperature))
str(spdat)

# for(i in 421868:nrow(spdat)){
#   
#   
#   spdat$web[i]<- as.character(fact_1$webfiles[which(fact_1$connectance==spdat$Connectance[i])][1])
#   
#   print(i)
#   
#   
#   
# }
# 
#spdat<- sp_dat1
mydir = 'datasets_1'
myfiles = list.files(path=mydir, pattern="*.csv", full.names=TRUE)
#myfiles<-myfile
webfiles<-myfiles


#example dynamics of three species plant-pollinator system 
fact_2<- expand.grid(`web` =webfiles[1:153]) %>%
  as_tibble %>%
  mutate(`Nestedness`= 0,
         `Connectance` = 0,
         `network_size`=0)
# a loop to go over the rows of the dataframe to calculate nestedness and network size for each network
for (r in 1:nrow(fact_2)){
  g<-adj.mat(myfiles[which(myfiles == fact_2$web[r])]) #network web names
  Aspecies<- dim(g)[2] # no of animal species
  Plantspecies<- dim(g)[1]
  
  fact_2$network_size[r]<-Aspecies+Plantspecies
  fact_2$Nestedness[r]<-  round(nestedness_NODF(g),2)
  fact_2$Connectance[r] <- Connectance(g)
}

fact_2 <- fact_2 %>% filter(network_size < 105 )

webfiles<-as.character(fact_2$web)

web<-webfiles[40:86]
fact<-expand.grid(Temperature=seq(12,45,0.5),
                  h2=c(0,0.4),
                  gamma=1.5,
                  var=c("high","low"),
                  `web` =seq(1,47,1),
                  `replicates`=1+(1:1)*100) %>% mutate(Dominant_eigenvalue=0,
                                                       Average_robustness=0,
                                                       network_size=0,
                                                       connectance=0,
                                                       nestedness=0)

for(i in 1:nrow(fact)){
  
 temp<-spdat %>% filter(Temperature == fact$Temperature[i],heritability == fact$h2[i],
                  Individual_variation == fact$var[i], webname == fact$web[i])
  
density<-temp$density 
traitvals<-temp$trait_values 

variation <- fact$var[i]


g<-adj.mat(myfiles[which(myfiles ==  web[fact$web[i]])]) #network web names
# g<-g[-1,-1] 


Aspecies<- dim(g)[2] # no of animal species
Plantspecies<- dim(g)[1] # no of plant species
degree.plants<-degree.animals<-numeric()
for(r in 1:Plantspecies){
  degree.plants[r]<-sum(g[r,])} # degree of plants
for(x in 1:Aspecies){
  degree.animals[x]<-sum(g[,x]) # degree of animals
}

density_N <-density[1:Aspecies]
density_P <-density[1:Plantspecies]

trait_val_N<-traitvals[1:Aspecies]
trait_val_P<-traitvals[1:Plantspecies]     

if(variation == "high"){
  s<- runif((Aspecies+Plantspecies), 0.05,0.09) 
  
}else{
  s<- runif((Aspecies+Plantspecies), 0.005,0.009) 
}


Jij<- matrix(0, nrow=Aspecies,ncol=Plantspecies)
Jij_2<-matrix(0, nrow=Plantspecies,ncol=Aspecies)
for(r in 1:Aspecies){
  for(l in 1:Plantspecies){
    #
    m.temp<-list(ma=trait_val_N[r],mp=trait_val_P[l])
    sigma1<-list(sa=sqrt(s[r]),sp=sqrt(s[(Aspecies)+l]))
    temp1<-gausquad.animals(m=m.temp,sigma=sigma1,w=7,h=0.25,np=density_P[l],na=density_N[r],
                            mut.strength=1.5, points=7,
                            mat=g[l,r],
                            degree.animal = degree.animals[r], vec_na=density_N, vec_np=density_P)
     Jij[r,l]<- temp1$J
    
  }

}

for(k in 1:Plantspecies){
  for(o in 1:Aspecies){
    m2.temp<-list(ma=trait_val_N[o],mp=trait_val_P[k])
    sigma2<-list(sa=sqrt(s[o]),sp=sqrt(s[(Aspecies+k)]))
    temp2<-gausquad.plants(m=m2.temp,sigma=sigma2,w=7,h=0.25,np=density_P[k],na=density_N[o],
                           mut.strength=1.5,
                           points=7,mat=g[k,o], 
                           degree.plant =degree.plants[k], vec_na=density_N, vec_np=density_P)
       Jij_2[k,o]<- temp2$J
  }

}

Amatrix <- mat.comp(g)$Amatrix  #competition matrix , aij, for animals
Pmatrix <- mat.comp(g)$Pmatrix  #competition matrix , aij, for plants
Jc_p<-matrix(0, nrow=Plantspecies, ncol= Plantspecies)
Jc_a <- matrix(0, nrow=Aspecies, ncol= Aspecies)
for(t in 1:Aspecies){
  for(v in 1:Aspecies){
    
    Jc_a[t,v]<- -Amatrix[t,v]*density_N[t]
  }
}

for(u in 1:Plantspecies){
  for(w in 1:Plantspecies){
    
    Jc_p[u,w]<- -Pmatrix[u,w]*density_P[u]
  }
}

#THE JACOBIAN MATRIX
Jacobian<-cbind(rbind( Jc_a, Jij_2),rbind(Jij,Jc_p))
    
eJ<-Re(eigen(Jacobian, only.values =T)$values)

average_robustness<- exp(mean(log(abs(eJ)))) 
dominant.eigenvalue<-eJ[1]

fact$Average_robustness[i]<- average_robustness
fact$Dominant_eigenvalue[i]<-dominant.eigenvalue
fact$network_size[i]<- Plantspecies+Aspecies
fact$connectance[i]<- Connectance(g)
fact$nestedness[i]<-nestedness_NODF(g)
print(i)
}

fact$var<-plyr::revalue(fact$var, c("high"= "high variation", "low" = "low variation"))

(r10<-fact %>% ggplot(aes(x=Temperature,y=Dominant_eigenvalue,color=factor(h2)))+
  geom_point(size=4,alpha =0.5)+
  theme_bw()+
  labs(color="heritability")+
  scale_color_brewer(palette="Dark2")+
  facet_wrap(.~var)+
    ylab("Re(dominant eigenvalue)")+
  xlab("Temperature"))




(r11<-fact %>% filter(Temperature <25) %>%  
  ggplot(aes(x=network_size,y=Dominant_eigenvalue,color=factor(h2)))+
  geom_point(size=4,alpha =0.5)+
  scale_color_brewer(palette="Dark2")+
  theme_bw()+
  labs(color="heritability")+
  ylab("Re(dominant eigenvalue)")+
  xlab("Network size")+
  scale_color_brewer(palette="Dark2")+
  geom_hline(yintercept =0, linetype = "dashed",alpha=0.5, lwd =2)+
  geom_smooth(method = "lm",formula = y~x, se=F, lwd=2)+
  facet_wrap(.~var))

(r12<-fact %>% filter(Temperature <25) %>%  
  ggplot(aes(x=nestedness,y=Dominant_eigenvalue,color=factor(h2)))+
  geom_point(size=4,alpha =0.5)+
  scale_color_brewer(palette="Dark2")+
  theme_bw()+
  labs(color="heritability")+
  ylab("Re(dominant eigenvalue)")+
  xlab("Nestedness (NODF)")+
  ggtitle("A")+
  scale_color_brewer(palette="Dark2")+
  geom_hline(yintercept =0, linetype = "dashed",alpha=0.5, lwd =2)+
  geom_smooth(method = "lm",formula = y~x, se=F, lwd=2)+
  facet_wrap(.~var))


grid_arrange_shared_legend(r12,n1, nrow = 2,ncol = 1)


fact %>% filter(Temperature <30) %>%  
  ggplot(aes(x=connectance,y=Dominant_eigenvalue,color=factor(h2)))+
  geom_point(size=4,alpha =0.5)+
  theme_bw()+
  geom_smooth(method = "lm",formula = y~x, se=F)+
  facet_wrap(.~var)


fact %>% filter(Temperature <30) %>%  
  ggplot(aes(x=nestedness,y=Dominant_eigenvalue,color=factor(h2)))+
  geom_point(size=4,alpha =0.5)+
  theme_bw()+
  geom_smooth(method = "lm",formula = y~x, se=F)+
  facet_wrap(.~var)





















#gaussian quadrature approximation used to numerically determine the integrals.
gausquad.animals_jac<-function(m,sigma,w,h,np,na,mut.strength,points,mat,degree.animal, vec_na,vec_np){
  
  
  temp2<-dat2<-jac1<-x2<-x3<-array(dim=c(points))
  if(mat == 0){
    return(list(G= 0, B = 0, J=0))
  }
  else if(mat == 1){
    #nodes oir points in the abscissa where the integral will be evaluated numerically
    
    z1<-gauss.quad.prob(points, dist = "normal", mu=m$ma, sigma =sigma$sa)$nodes #z'
    z2<-gauss.quad.prob(points, dist = "normal", mu=m$mp, sigma =sigma$sp)$nodes #z''
    
    #weights of the gaussian distribution given by mean trait value mu_i and its variance \sigma_i
    w1<-gauss.quad.prob(points, dist = "normal", 
                        mu=m$ma,sigma =sigma$sa)$weights #pi(z')
    w2<-gauss.quad.prob(points, dist = "normal", 
                        mu=m$mp,sigma =sigma$sp)$weights #pj(z'')
    
    
    #for the pairwise model however there are only two species interacting and hence i and j
    #or in other words the integral goes over z and z'
    for (i in 1: points){
      
      
      f <-  exp(-(z1[i]- z2)^2/w^2) # + 2*alpha*(sign(z1[i] - z2))*(1- exp(-(z1[i]-z2)^2/w^2)) + sign(alpha))
      
      temp2[i]<- sum(np*(mut.strength/degree.animal)*f/(1+h*np*(mut.strength/degree.animal)*f)*w2*w1[i])
      
      jac1[i] <- sum(na*(mut.strength/degree.animal)*f/(1+h*sum( np)*(mut.strength/degree.animal)*f)^2*w2*w1[i])
      
      
      dat2[i]<- sum( ((z1[i]-m$ma)*f*(mut.strength/degree.animal)*np/(1+h*np*(mut.strength/degree.animal)*f) )*w2*w1[i])
      
      
      
    }
    G = sum(temp2)
    B = sum(dat2) 
    J=sum(jac1)
    
    
    
    return(list(G= G, B = B, J=J))
  }
}

#gaussian quadrature approximation used to numerically determine the integrals.
gausquad.plants_jac<-function(m,sigma,w,h,np,na,mut.strength,points,mat,degree.plant, vec_na,vec_np){
  
  temp2<-dat2<-x3<-x4<-jac2<-array(dim=c(points))
  
  if(mat == 0){
    
    return(list(G= 0, B = 0, J=0))
    
  }
  else if (mat==1){
    #nodes oir points in the abscissa where the integral will be evaluated numerically
    z1<-gauss.quad.prob(points, dist = "normal", mu=m$mp, sigma =sigma$sp)$nodes #z'
    z2<-gauss.quad.prob(points, dist = "normal", mu=m$ma, sigma =sigma$sa)$nodes #z''
    
    #weights of the gaussian distribution given by mean trait value mu_i and its variance \sigma_i
    
    w1<-gauss.quad.prob(points, dist = "normal", 
                        mu=m$mp,sigma =sigma$sp)$weights #pi(z')
    w2<-gauss.quad.prob(points, dist = "normal", 
                        mu=m$ma,sigma =sigma$sa)$weights #pj(z'')
    
    
    #for the pairwise model however there are only two species interacting and hence i and j
    #or in other words the integral goes over z and z'
    for (i in 1: points){
      
      f <-  exp(-(z1[i]- z2)^2/w^2) # + 2*alpha*(sign(z1[i] - z2))*(1- exp(-(z1[i]-z2)^2/w^2)) + sign(alpha))
      
      temp2[i]<- sum(na*(mut.strength/degree.plant)*f/(1+h*na*(mut.strength/degree.plant)*f)*w2*w1[i])
      
      jac2[i]<- sum(np*(mut.strength/degree.plant)*f/(1+h*sum(na)*(mut.strength/degree.plant)*f)^2*w2*w1[i])
      
      
      dat2[i]<- sum( ((z1[i]-m$mp)*f*(mut.strength/degree.plant)*na/(1+h*na*(mut.strength/degree.plant)*f) )*w2*w1[i])
      
      # x4[i]<-sum(h*w2*exp(-(z1[i]-z2)^2/w^2))
    }
    
    G= sum(temp2)
    B = sum(dat2)
    J=sum(jac2)
  } 
  
  
  
  return(list(G= G, 
              B = B, J=J))
  
  
}
  