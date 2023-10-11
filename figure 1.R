rm(list=ls())
source("01_functions_cluster.R")
library(statmod) ## for integrating ordinary differential equations
require(tidyverse) ## for efficient data manipulation & plotting
library(cowplot) ## for arranging plots in a grid
library(dplyr)
library(readr)
library(beepr)
library(deSolve)
library(network)
library(sna)
library(ggplot2)
library(ggnet)



##### example simulations

mydir = 'datasets_1'
myfiles = list.files(path=mydir, pattern="*.csv", full.names=TRUE)
#myfiles<-myfile
webfiles<-myfiles



fact<-expand.grid(Temperature=28,
                  h2=c(0,0.4),
                  gamma=1,
                  var="high",
                  `web` = c(webfiles[67]),
                  `random_seed`=4327+(1:1)*100)



  
  g<-adj.mat(myfiles[which(myfiles == fact$web[1])]) #network web names
  # g<-g[-1,-1] 
  Aspecies<- dim(g)[2] # no of animal species
  Plantspecies<- dim(g)[1] # no of plant species
  nestedness<-nestedness_NODF(g)
  C<-Connectance(g)
  
  
  #parameters for modelling: intialisation of the model
  if(fact$var[1] == "high"){  sigma<-runif((Aspecies+Plantspecies), 0.02,0.09) 
  }else {sigma <- sigma<-runif((Aspecies+Plantspecies), 0.002,0.009)
  }#three species trait variance 
  Na <- runif( (Aspecies) , 1,1)
  Np<- runif( (Plantspecies) , 1,1)
  muA<- runif(Aspecies, 15, 30)  #initial mean phenotypic optimum trait values
  muP<- runif(Plantspecies, 15, 30)   #intial mean phenotypic optimum trait values
  
  
  
  
  #parameters below are taken from another paper Akesson et al 2021 Nat Comms.
  bw  <- 2
  aw  <- 0.1
  gi <- 1
  ki <-0.1 #mortality rate
  w<- 7#mutualism interaction width
  Temp<-fact$Temperature[1]<-28 #seq(Tmin, Tmax, by=(Tmax-Tmin)/50)
  h2<-fact$h2[1] #heritability of trait variance 
  Amatrix <- mat.comp(g)$Amatrix  #competition matrix , aij, for animals
  Pmatrix <- mat.comp(g)$Pmatrix  #competition matrix , aij, for plants
  mut.strength=2 #average mutualistic strength
  
  params<-list(matrix=g,bw=bw,aw=aw,h2=h2,w=w,Amatrix=Amatrix,Pmatrix=Pmatrix,
               gi=gi,ki=ki,Temp=Temp, sigma=sigma,
               mut.strength=1.5)
  
  
  ic<-c(Na,Np,muA,muP) ## initial conditions coerced into a vector
  
  t1<-plot_snapshot(Na=Na,Np=Np, m=c(muA,muP), sigma = sigma,
                    limits = c(14,31), res = 1001)
  
  
  
  tmax <- 1e4## time to integrate equations fors
  sol<-ode(func=eqs, y=ic, parms=params, 
           times=seq(0, tmax, by=tmax/100)) %>% ## solve ODEs
    organize_results(params) 
  
  h0<-sol %>% plot_density()
  #for time point 1
h0
  sol_1<- sol %>% filter(time == 1e4)
  Na<- (sol_1 %>% filter(type=="N"))$v
  Np<-(sol_1 %>% filter(type =="P"))$v
  ma<-(sol_1 %>% filter(type == "ma"))$v
  mp<-(sol_1 %>% filter(type == "mp"))$v
  t2<- plot_snapshot(Na=Na, Np=Np, m = c(ma,mp), sigma = sigma, limits = c(15,31),res = 1001)
  t2
  #fa<-functional_diversity(n = Na$v, m = ma$v,s = sigma[1:Aspecies], q=2, delta=0.001)
  #fp<-functional_diversity(n = Np$v, m = mp$v,s = sigma[1:Plantspecies], q=2, delta=0.001)
  
  
 
 
 
 #### heritability -0.4 
 
 
 
 h2<-fact$h2[2]<-0.25 #heritability of trait variance 
 Amatrix <- mat.comp(g)$Amatrix  #competition matrix , aij, for animals
 Pmatrix <- mat.comp(g)$Pmatrix  #competition matrix , aij, for plants
 mut.strength=1 #average mutualistic strength
 
 
 Na <- runif( (Aspecies) , 1,1)
 Np<- runif( (Plantspecies) , 1,1)
 muA<- runif(Aspecies, 15, 30)  #initial mean phenotypic optimum trait values
 muP<- runif(Plantspecies, 15, 30)   #intial mean phenotypic optimum trait values
 
 
 ic<-c(Na,Np,muA,muP) ## initial conditions coerced into a vector
 
 t3<-plot_snapshot(Na=Na,Np=Np, m=c(muA,muP), sigma = sigma,limits = c(14,31), res = 1001)
 
 bw  <- 2
 aw  <- 0.1
 gi <- 1
 ki <-0.1 #mortality rate
 w<- 7#mutualism interaction width
 Temp<-fact$Temperature[1]<-28  #seq(Tmin, Tmax, by=(Tmax-Tmin)/50)
 h2<-fact$h2[2] #heritability of trait variance 
 Amatrix <- mat.comp(g)$Amatrix  #competition matrix , aij, for animals
 Pmatrix <- mat.comp(g)$Pmatrix  #competition matrix , aij, for plants
 mut.strength=1.5 #average mutualistic strength
 
 params<-list(matrix=g,bw=bw,aw=aw,h2=h2,w=w,Amatrix=Amatrix,Pmatrix=Pmatrix,
              gi=gi,ki=ki,Temp=Temp, sigma=sigma,
              mut.strength=1.5)
 
 
 
 tmax <- 1e4## time to integrate equations fors
 sol<-ode(func=eqs, y=ic, parms=params, times=seq(0, tmax, 
                                                  by=tmax/100)) %>% ## solve ODEs
   organize_results(params) 
 
 h4<-sol %>% plot_density()
 #for time point 1
 
 sol_1<-  sol %>% filter(time == 1e3)
 Na<- (sol_1 %>% filter(type=="N"))$v
 Np<-(sol_1 %>% filter(type =="P"))$v
 ma<-(sol_1 %>% filter(type == "ma"))$v
 mp<-(sol_1 %>% filter(type == "mp"))$v
 
 t4<- plot_snapshot(Na=Na, Np=Np, m = c(ma,mp), 
                    sigma = sigma, limits = c(15,30),res = 1001)
 
 t4
 
 net = network(g, directed = FALSE)
  #colorblind palette
 cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
 
 #vertex names 
 names<-network.vertex.names(net)
 net %v% "groups" = ifelse( names[1:sum(dim(g))] %in% c( as.character(seq(1:dim(g)[1])) ), "plants", "animals")
 net %v% "color" = ifelse(net %v% "groups" == "plants", "#0072B2", "#E69F00" )
 #ggnet2(net,  mode="circle",  color ="groups", edge.size = 1,edge.alpha = 1, edge.color = "black", edge.lty = 1)
 
 webg<-ggnet2(net, mode="circle", size=5, 
              edge.size = 1.1,max_size =12, 
              color ="color",edge.alpha = 1, legend.position = "")
 
 
 ggpubr::ggarrange(webg,t1,t2,h0,
                   t1,t4,h4,
                   ncol=4,nrow=2,
                   labels = c("A", "B", "C",
                              "D",
                              "E", "F","G"))
 