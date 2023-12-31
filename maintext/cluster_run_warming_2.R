rm(list=ls())
source("functions.R")
library(foreach)
library(doMC)
library(statmod)
require(deSolve) ## for integrating ordinary differential equations
require(tidyr) ## for efficient data manipulation & plotting
library(dplyr)
#library(readr)
library(beepr)
#library(viridis)
numcores<- 20
registerDoMC(numcores)



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
##### example simulations


fact<-expand.grid(Temperature=seq(12,45,0.5),
                  h2=c(0,0.4),
                  gamma=1.5,
                  var=c("high","low"),
                  `web` = webfiles[1:86],
                  `replicates`=1+(1:1)*100) 


output<-foreach(r = 1:nrow(fact))%dopar%{
  
  g<-adj.mat(myfiles[which(myfiles == fact$web[r])]) #network web names
  # g<-g[-1,-1] 
  
  
  Aspecies<- dim(g)[2] # no of animal species
  Plantspecies<- dim(g)[1] # no of plant species
  
  #degree of plants and anichmals
  degree.plants<-degree.animals<-numeric()
  for(i in 1:Plantspecies){
    degree.plants[i]<-sum(g[i,])} # degree of plants
  for(j in 1:Aspecies){
    degree.animals[j]<-sum(g[,j]) # degree of animals
  }
  
  degree<-c(degree.animals, degree.plants)
  nestedness<-nestedness_NODF(g)
  C<-Connectance(g)
  network_size<-Aspecies+Plantspecies
  
  #parameters for modelling: intialisation of the model
  if(fact$var[r] == "high"){  sigma<-runif((Aspecies+Plantspecies), 0.05,0.09) 
  }else {sigma <- runif((Aspecies+Plantspecies), 0.005,0.009)
  }#three species trait variance 
  na <- runif( (Aspecies) , 1,1)
  np<- runif( (Plantspecies) , 1,1)
  muA<- runif(Aspecies, 10, 30)  #initial mean phenotypic optimum trait values
  muP<- runif(Plantspecies, 10, 30)   #intial mean phenotypic optimum trait values
  
  
  
  
  #parameters below are taken from another paper Akesson et al 2021 Nat Comms.
  bw  <- 2
  aw  <- 0
  gi <- 1
  ki <-0.1 #mortality rate
  w<- 7#mutualism interaction width
  Temp<-fact$Temperature[r]  #seq(Tmin, Tmax, by=(Tmax-Tmin)/50)
  h2<-fact$h2[r] #heritability of trait variance 
  Amatrix <- mat.comp(g)$Amatrix  #competition matrix , aij, for animals
  Pmatrix <- mat.comp(g)$Pmatrix  #competition matrix , aij, for plants
  mut.strength=1.5 #average mutualistic strength
  
  params<-list(matrix=g,bw=bw,aw=aw,h2=h2,w=w,Amatrix=Amatrix,Pmatrix=Pmatrix,web=fact$web[r],
               gi=gi,ki=ki,Temp=Temp, sigma=sigma,A=Aspecies,P=Plantspecies,degree=degree,
               var=fact$var[r],g=g,
               mut.strength=1.5,nestedness=nestedness,C=C,network_size=network_size)
  
  ic<-c(na,np,muA,muP) ## initial conditions coerced into a vector
  
  tmax <- 1e4## time to integrate equations fors
  
  out<-cluster_run( params = params,ic = ic,tmax = tmax)
  
}

#save(output, file ="cluster_warming.RData")
