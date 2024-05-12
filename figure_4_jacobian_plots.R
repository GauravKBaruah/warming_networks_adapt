rm(list=ls())
source("01_functions_cluster.R")
library(statmod)
require(tidyr)
library(cowplot)
library(readr)
library(beepr)
library(viridis)
library(ggdist)
library(dplyr)
library(readr)
library(beepr)
library(viridis)
library(ggplot2)

#functions used for analysis of the results for the dominant eigenvalue
#gaussian quadrature approximation used to numerically determine the double integrals in the equations 1 and 2.
gausquad.animals_jac<-function(m,sigma,w,h,np,na,mut.strength,points,mat,degree.animal, vec_na,vec_np){
  
  
  temp2<-dat2<-jac1<-x2<-x3<-array(dim=c(points))
  if(mat == 0){ #if there is no existing species interaction from the adjacency matrix
    return(list(G= 0, B = 0, J=0))
  }
  else if(mat == 1){ # if there is an existing species interaction from the ajdacency matrix
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
#returns an approximate value of the double integrals in the main-text
gausquad.plants_jac<-function(m,sigma,w,h,np,na,mut.strength,points,mat,degree.plant, vec_na,vec_np){
  
  temp2<-dat2<-x3<-x4<-jac2<-array(dim=c(points))
  
  if(mat == 0){ #if there is no existing species interaction from the ajdacency matrix
    
    return(list(G= 0, B = 0, J=0))
    
  }
  else if (mat==1){ #if there is an existing species interaction from the ajdacency matrix
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
      
    }
    
    G= sum(temp2)
    B = sum(dat2)
    J=sum(jac2)
  } 
  
  
  
  return(list(G= G, 
              B = B, J=J))
  
  
}


####### FOR jacobian stability analysis#################
load("all_data_species.RData") # load the species level data from simulations
spdat<-sp_dat
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


#spdat<- sp_dat1
mydir = 'datasets_1'
myfiles = list.files(path=mydir, pattern="*.csv", full.names=TRUE)
webfiles<-myfiles


#data frame for all the species interaction networks
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

#isolating only those matrices that has number of species less than 105
fact_2 <- fact_2 %>% filter(network_size < 105 )
webfiles<-as.character(fact_2$web)


#empty dataframe to organise further results and analysis for each treatment parameter:
#temperature,
#heritability, h2
#var:individual variation
#mutualistic strenght = 1.5
fact<-expand.grid(Temperature=seq(12,40,0.5),
                  h2=c(0, 0.2, 0.4),
                  gamma=1.5,
                  var=c("high","low"),
                  `web` =unique(spdat$webname),
                  `replicates`=1+(1:1)*100) %>% mutate(Dominant_eigenvalue=0,
                                                       Average_robustness=0,
                                                       network_size=0,
                                                       connectance=0,
                                                       nestedness=0)

#loop that goes over all these treatments and calculates the dominant eigenvalues at equilibrium
for(i in 1:nrow(fact)){
  
  temp<-spdat %>% filter(Temperature == fact$Temperature[i],heritability == fact$h2[i],
                         Individual_variation == fact$var[i], webname == fact$web[i])
  
  #extracts the equilibrium density and trait values of all species for a particular web, for a particular treatment 
  #of temperature, h2, variation
  density<-temp$density 
  traitvals<-temp$trait_values 
  variation <- fact$var[i]
  
  #adjacency matrix of the web in question
  g<-adj.mat(myfiles[which(myfiles ==  fact$web[i])]) #network web names

  
  
  Aspecies<- dim(g)[2] # no of animal species
  Plantspecies<- dim(g)[1] # no of plant species
  degree.plants<-degree.animals<-numeric()
  for(r in 1:Plantspecies){
    degree.plants[r]<-sum(g[r,])} # degree of plants
  for(x in 1:Aspecies){
    degree.animals[x]<-sum(g[,x]) # degree of animals
  }
  
  #density and trait values at equilibrium for animal and plants
  density_N <-density[1:Aspecies]
  density_P <-density[1:Plantspecies]
  trait_val_N<-traitvals[1:Aspecies]
  trait_val_P<-traitvals[1:Plantspecies]     
  
  #trait variation used in the main-text
  if(variation == "high"){
    s<- runif((Aspecies+Plantspecies), 0.05,0.09) 
    }else{
    s<- runif((Aspecies+Plantspecies), 0.005,0.009) 
  }
  
  #calculation of the jacobian matrices using the gaussian quadrature method.
  Jij<- matrix(0, nrow=Aspecies,ncol=Plantspecies)
  Jij_2<-matrix(0, nrow=Plantspecies,ncol=Aspecies)
  for(r in 1:Aspecies){
    for(l in 1:Plantspecies){
      #
      m.temp<-list(ma=trait_val_N[r],mp=trait_val_P[l])
      sigma1<-list(sa=sqrt(s[r]),sp=sqrt(s[(Aspecies)+l]))
      temp1<-gausquad.animals_jac(m=m.temp,sigma=sigma1,w=7,h=0.25,np=density_P[l],na=density_N[r],
                                  mut.strength=1.5, points=7,
                                  mat=g[l,r],
                                  degree.animal = degree.animals[r], vec_na=density_N, vec_np=density_P)
      Jij[r,l]<- temp1$J
      
    }
    
  }
  #calculation of the jacobian matrices using the gaussian quadrature method.
  for(k in 1:Plantspecies){
    for(o in 1:Aspecies){
      m2.temp<-list(ma=trait_val_N[o],mp=trait_val_P[k])
      sigma2<-list(sa=sqrt(s[o]),sp=sqrt(s[(Aspecies+k)]))
      temp2<-gausquad.plants_jac(m=m2.temp,sigma=sigma2,w=7,h=0.25,np=density_P[k],na=density_N[o],
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
      
      Jc_a[t,v]<- -Amatrix[t,v]*density_N[t] #jacobian of competition matrix of animals
    }
  }
  
  for(u in 1:Plantspecies){
    for(w in 1:Plantspecies){
      
      Jc_p[u,w]<- -Pmatrix[u,w]*density_P[u] #jacobian of competition matrix of plants (See supplementary)
    }
  }
  
  #THE FINAL JACOBIAN MATRIX that has competition as well as mutualistic matrix part.
  Jacobian<-cbind(rbind( Jc_a, Jij_2),rbind(Jij,Jc_p))
  
  #eigenvalues of the jacobian matrix
  eJ<-Re(eigen(Jacobian, only.values =T)$values)
  
  average_robustness<- exp(mean(log(abs(eJ)))) #robustness of the network
  dominant.eigenvalue<-eJ[1] #dominant eigenvalue of the network
  
  fact$Average_robustness[i]<- average_robustness
  fact$Dominant_eigenvalue[i]<-dominant.eigenvalue
  fact$network_size[i]<- Plantspecies+Aspecies
  fact$connectance[i]<- Connectance(g)
  fact$nestedness[i]<-nestedness_NODF(g)
  print(i)
}


#renaming variables and factors
fact$var<-plyr::revalue(fact$var, c("high"= "high variation", "low" = "low variation"))

#for figure 4 A of the main-text nestedness vs. dominant eigenvalue
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

#for figure S3 in the appendix of temperature and dominant eigenvalue
(r10<-fact %>% ggplot(aes(x=Temperature,y=Dominant_eigenvalue,color=factor(h2)))+
    geom_point(size=4,alpha =0.5)+
    theme_bw()+
    labs(color="heritability")+
    scale_color_brewer(palette="Dark2")+
    facet_wrap(.~var)+
    ylab("Re(dominant eigenvalue)")+
    xlab("Temperature"))

#
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

grid_arrange_shared_legend(r12,n1, nrow = 2,ncol = 1)



#connectance vs. dominant eigenvalue
fact %>% filter(Temperature <30) %>%  
  ggplot(aes(x=connectance,y=Dominant_eigenvalue,color=factor(h2)))+
  geom_point(size=4,alpha =0.5)+
  theme_bw()+
  geom_smooth(method = "lm",formula = y~x, se=F)+
  facet_wrap(.~var)









