#R Program to Read in Greg's value data and deflate to common year
rm(list=ls())
################################################################################
library(data.table)
library(tidyverse)
library(dplyr)
library(plyr)
library(readr)
library(readxl)
library(bea.R)
library(ggplot2)
library(gridExtra)
library(copula)
library(here)
source('beakey.R')  #BEA Key for building deflator
################################################################################
#FUNCTIONS
################################################################################
source('Misc_functions.R')  #External Functions Provided by JOE
###############################################################################################################
#Function to Generate Samples
dens_sample<-function(X,samples){
  data1=density(X,from=min(X),to=max(X))
  tmpsim=sample(data1$x,samples,prob=data1$y,replace=T)
  
  return(tmpsim)
}
#####################################################################################
#Function for Simulated Values Based on the Density Function (dens_sample)
###########################################################################################################
sims<-function(X,N){
  valsim=dens_sample(X$val2022,N)  #Value
  agesim=dens_sample(X$age2,N)             #Age
  lensim=dens_sample(X$len,N)              #Length
  vhpsim=dens_sample(X$vhp,N)            #VHP
  gtsim<-dens_sample(X$gtons,N)           #GTONS
  YI<-cbind(valsim,lensim,vhpsim,agesim,gtsim)
  
  colnames(YI)<-c("valsim","agesim","lensim","vhpsim","gtsim")
  
  return(YI)
}
#############################################################################################
#Function
hists<-function (X0, X1){
  
  valplot<-ggplot(as.data.frame(X0$val2022),aes(x=X0$val2022))+
    geom_histogram(bins=50, color="black", fill="gray")+
    labs(title="Historical Value",
         x="Vessel Value ($2022)",
         y="count")+
    theme(plot.title=element_text(hjust=0.5))
  
  valplot2<-ggplot(as.data.frame(X1),aes(x=valsim))+
    geom_histogram(bins=50, color="black", fill="gray")+
    labs(title="Simulated Value",
         x="Vessel Value ($2022)",
         y="count")+
    theme(plot.title=element_text(hjust=0.5))
  
  
  grid1<-grid.arrange(valplot, valplot2)
  
  return(grid1)
  
}
#Define function which can be used to bind vessel attributes
Xval<-function(data){
  X1<-data$val2022
  X2<-data$len
  X3<-data$vhp
  X4<-data$age2
  X5<-data$gtons
  XFINAL<-cbind(X1,X2,X3,X4,X5)
  
  names<-c("value","len","vhp","age2","tons")
  colnames(XFINAL)<-names
  return(XFINAL)
}
######################################################################################################
histcopula<-function(X,XSIM,s_scale){
  ###############################################################################
  # historical copula  binding approach
  #X=Historical Data
  #XSIM=Simulated DATA from call to function sims
  ###############################################################################
  XTOC1<-as.data.frame(Xval(X)) #Call other function
  X<-XTOC1
  
  (n=length(X$value))
  (N=length(XSIM$valsim))
  (nreps=floor(N/n)+1)
  
  Y=X
  
  HCOP=cbind(val1=rep(X$value,nreps),len1=rep(X$len,nreps), vhp1=rep(X$vhp,nreps),age1=rep(X$age2,nreps),gtons=rep(X$tons,nreps))
  HCOP=HCOP[1:N,]
  
  YCE = XSIM
  
  j=1
  for(j in 1:ncol(XSIM)){
    y=sort(XSIM[,j])
    YCE[,j]=y[rank(HCOP[,j], ties.method = "first")]
  }

  YCE<-as.data.frame(YCE)
  colnames(YCE)<-c("value","len","vhp","age2","tons")
  

  ###############################################################################
  # "smooth" the historical copula
  ###############################################################################
  (sdevs=apply(HCOP,2,sd))
  
  #s_scale = 0.10
  s_scale = s_scale
  #s_scale = 0.50
  
  # Potential correlation between z values used in Copula smoothing
  
  rho
  tmp=FactorYcZc(XSIM,rhos=0)
  #tmp=FactorYcZc(YI,rhos=sqrt(0.75*rho))
  #tmp=FactorYcZc(YI,rhos=sqrt(rho))
  
  ZC=tmp$Zc
  # cor(ZC)
  # 
  # 
  # cor(HCOP)
  ########################################################
  HCOP0=HCOP
  hpick=sample(1:nrow(HCOP),nrow(HCOP))
  HCOP2=HCOP0[hpick,]
  ########################################################
  
  #for(j in 1:2) HCOP[,j] = HCOP[,j] + s_scale*rnorm(N,0,sdevs[j])
  
  for(j in 1:ncol(HCOP)) HCOP[,j] = HCOP0[,j] + s_scale*ZC[,j]*sdevs[j]
  
  YCE = XSIM
  
  j=1
  for(j in 1:ncol(XSIM)){
    y=sort(XSIM[,j])
    YCE[,j]=y[rank(HCOP[,j], ties.method = "first")]
  }
  cor(YCE)
  cor(HCOP)
  
  YCE<-as.data.frame(YCE)
  colnames(YCE)<-c("value","len","vhp","age2","tons")  
  # cor(HCOP0)
  # cor(HCOP2)
  # cor(HCOP)
  # 
  # 
  # cor(YI)
  # cor(HCOP)
  
  
  A<-list(X,YCE)
  return(A)
}
######################################################################################################
#Function to Graph Copula Points
simpoints<-function(Y,YSIM){
  
  x11()
  plot(YSIM$value~YSIM$len,main='Bound with HCOP' ,xlim=c(min(Y1$len),max(Y1$len)),ylim=c(min(Y1$value),max(Y1$value)))
  points(Y1$value~Y1$len,pch=20,col=2,cex=3)
  
  x11()
  plot(YSIM$value~YSIM$vhp,main='Bound with HCOP',xlim=c(min(Y1$vhp),max(Y1$vhp)),ylim=c(min(Y1$value),max(Y1$value)))
  points(Y1$value~Y1$vhp,pch=20,col=2,cex=3)
  
  
  x11()
  plot(YSIM$value~YSIM$age2,main='Bound with HCOP',xlim=c(min(Y1$age2),max(Y1$age2)),ylim=c(min(Y1$value),max(Y1$value)))
  points(Y1$value~Y1$age2,pch=20,col=2,cex=3)
  
  x11()
  plot(YSIM$value~YSIM$tons,main='Bound with HCOP',xlim=c(min(Y1$tons),max(Y1$tons)),ylim=c(min(Y1$value),max(Y1$value)))
  points(Y1$value~Y1$tons,pch=20,col=2,cex=3)
}
######################################################################################################
#End of Functions
###############################################################################################################
set.seed(1001)
(N=2000)
###########################################################################################################
#Extract Data and Place in Dataframes
ves_val<- read_excel("V:/Ardini_Cost_Survey/SAS_Cost_Survey/Profitability_Profiles/Calculate_Profit/Depreciation/survey_response_2015_2022_char_gear.xlsx",sheet="Data")
colnames(ves_val)[3]<-"YEAR"
ves_val<-subset(ves_val, YEAR==2015 | YEAR==2022)
################################################################################
#GDP Deflator
#calculate from BEA Data
SpecList<-list('UserID'= beaKey,
               'Method' = 'GetData',
               'datasetname' = 'NIPA',
               'TableName' = 'T10109',
               'RowNumber' = '10',
               'Frequency'='A',
               'Year'= '2015,2022'
)

GDP<-beaGet(SpecList, asWide=FALSE)
GDP<-GDP[(LineNumber==1),c("TimePeriod","DataValue")]
colnames(GDP)<-c("YEAR","GDP")
GDP$GDPD=GDP$GDP/100
GDP[,]<-lapply(GDP, function(x) type.convert(as.character(x), as.is=TRUE))
gdpd2022=subset(GDP, YEAR==2022)
GDP$GDPD2022=GDP$GDPD/gdpd2022$GDPD
GDP2<-GDP[,c(1,4)]
###############################################################################
#Join Vessel Value Data to GDP data
ves_val<-join(ves_val,GDP2, by="YEAR", type="inner")
ves_val$vessel_val=as.numeric(ves_val$vessel_val)
ves_val$val2022=ves_val$vessel_val/ves_val$GDPD2022
ves_val<-subset(ves_val, val2022>0) #for 2015 and 2022 data, this step drops 240 obs. (576 down to 336)
################################################################################
ves_val$age2=100-(ves_val$vessel_age)
ves_val<-subset(ves_val, toc==1 | toc==2 | toc==3) #for 2015 and 2022 data, this step drops 12 obs. (336 down to 324)
ves_val$primary_gear=ifelse(is.na(ves_val$primary_comm_gear),"Other",ves_val$primary_comm_gear)
#############################################################################################################
#Subset vessel value into separate sets. TOC3 = Steel Hulled Vessels
TOC3<-subset(ves_val, toc==3)

TOC3_Dredge<-subset(TOC3, primary_gear=="Dredge_Scallop")
hist(TOC3_Dredge$val2022,breaks=10)
TOC3_Dredge<-subset(TOC3_Dredge, val2022 < 2000000) #for 2015 and 2022, removes 1 outlier obs.(13 down to 12)

TOC3_Trawl<-subset(TOC3, primary_gear=="Trawl")
hist(TOC3_Trawl$val2022, breaks=10)
TOC3_Trawl<-subset(TOC3_Trawl, val2022<7000000) #for 2015 and 2022, removes 0 outlier obs.(14)

TOC3_Other<-subset(TOC3, primary_gear != "Trawl" & primary_gear != "Dredge_Scallop")
hist(TOC3_Other$val2022, breaks=10)
TOC3_Other<-subset(TOC3_Other, val2022<4000000) #for 2015 and 2022, removes 1 outlier obs.(13 down to 12)
##################################################################################################
#TOC2=Fiberglass Hulled Vessels
TOC2<-subset(ves_val, toc==2)

TOC2_HG<-subset(TOC2, primary_gear=="Handgear")
hist(TOC2_HG$val2022, breaks=20)
TOC2_HG<-subset(TOC2_HG, val2022<600000) #for 2015 and 2022, removes 2 outlier obs. (58 down to 56)

TOC2_POT<-subset(TOC2, primary_gear=="Pot/Trap")
hist(TOC2_POT$val2022, breaks=30)
TOC2_POT<-subset(TOC2_POT, val2022<1050000) #for 2015 and 2022, removes 2 outlier obs. (155 down to 153)

TOC2_TRAWL<-subset(TOC2, primary_gear=="Trawl")
hist(TOC2_TRAWL$val2022, breaks=10) #for 2015 and 2022, removes 0 outlier obs. (11)

TOC2_OTHER<-subset(TOC2, primary_gear!="Trawl" & primary_gear!="Pot/Trap"
                   & primary_gear!="Handgear") 
hist(TOC2_OTHER$val2022, breaks=10) #for 2015 and 2022, removes 0 outlier obs. (41)
################################################################################
#TOC1 = Wood Vessels
TOC1<-subset(ves_val, toc==1)
hist(TOC1$val2022, breaks=10) #for 2015 and 2022, removes 0 outlier obs. (19)
########################################################################################
#Simulate Data from Distributions for each hull type.
Y_TOC1<-as.data.frame(sims(TOC1,N))
Y_TOC2_POT<-as.data.frame(sims(TOC2_POT,N))
Y_TOC2_TRAWL<-as.data.frame(sims(TOC2_TRAWL,N))
Y_TOC2_OTHER<-as.data.frame(sims(TOC2_OTHER,N))
Y_TOC2_HG<-as.data.frame(sims(TOC2_HG,N))
Y_TOC3_TRAWL<-as.data.frame(sims(TOC3_Trawl,N))
Y_TOC3_DREDGE<-as.data.frame(sims(TOC3_Dredge,N))
Y_TOC3_OTHER<-as.data.frame(sims(TOC3_Other,N))
################################################################################################################
#Plot histograms comparison of actual data and simulated data
TOC1_grid<-hists(TOC1,Y_TOC1)
TOC2_POT_grid<-hists(TOC2_POT,Y_TOC2_POT)
TOC2_TRAWL_grid<-hists(TOC2_TRAWL,Y_TOC2_TRAWL)
TOC2_OTHER_grid<-hists(TOC2_OTHER,Y_TOC2_OTHER)
TOC2_HG_grid<-hists(TOC2_HG,Y_TOC2_HG)
TOC3_TRAWL_grid<-hists(TOC3_Trawl,Y_TOC3_TRAWL)
TOC3_DREDGE_grid<-hists(TOC3_Dredge,Y_TOC3_DREDGE)
TOC3_OTHER_grid<-hists(TOC3_Other,Y_TOC3_OTHER)
#####################################################################################################################
#Generate Copulas. The 0.5 below can be changed to get different copula
#TOC1
Y0<-histcopula(TOC1,Y_TOC1,0.5)
Y1<-Y0[[1]]
YSIM<-Y0[[2]]
#Note: The function simpoints is designed to immediately graph the plots
#comparing the actual points and simulated points values.
#They should be cleared from the workspace when done
simpoints(Y1,YSIM)
#The simulated points are saved in an R Data file for use in the Value model
save(YSIM, file="./data/TOC1_simulated_data.RData")
###################################################################################
#TOC2 Pot
Y0<-histcopula(TOC2_POT,Y_TOC2_POT,0.5)
Y1<-Y0[[1]]
YSIM<-Y0[[2]]
simpoints(Y1,YSIM)
save(YSIM, file="./data/TOC2_POT_simulated_data.RData")
########################################################################################
#TOC2 HG
Y0<-histcopula(TOC2_HG,Y_TOC2_HG,0.5)
Y1<-Y0[[1]]
YSIM<-Y0[[2]]
simpoints(Y1,YSIM)
save(YSIM, file="./data/TOC2_HG_simulated_data.RData")
#######################################################################################
#TOC2 TRAWL
Y0<-histcopula(TOC2_TRAWL,Y_TOC2_TRAWL, 0.5)
Y1<-Y0[[1]]
YSIM<-Y0[[2]]
simpoints(Y1,YSIM)
save(YSIM, file="./data/TOC2_TRAWL_simulated_data.RData")
#######################################################################################
#TOC2 OTher
Y0<-histcopula(TOC2_OTHER,Y_TOC2_OTHER, 0.5)
Y1<-Y0[[1]]
YSIM<-Y0[[2]]
simpoints(Y1,YSIM)
save(YSIM, file="./data/TOC2_OTHER_simulated_data.RData")
#######################################################################################
#TOC3 TRAWL
Y0<-histcopula(TOC3_Trawl,Y_TOC3_TRAWL, 0.5)
Y1<-Y0[[1]]
YSIM<-Y0[[2]]
simpoints(Y1,YSIM)
save(YSIM, file="./data/TOC3_TRAWL_simulated_data.RData")
#######################################################################################
#TOC3 Dredge
Y0<-histcopula(TOC3_Dredge,Y_TOC3_DREDGE, 0.5)
Y1<-Y0[[1]]
YSIM<-Y0[[2]]
simpoints(Y1,YSIM)
save(YSIM, file="./data/TOC3_DREDGE_simulated_data.RData")
#######################################################################################
#TOC3 OTHER
Y0<-histcopula(TOC3_Other,Y_TOC3_OTHER, 0.5)
Y1<-Y0[[1]]
YSIM<-Y0[[2]]
simpoints(Y1,YSIM)
save(YSIM, file="./data/TOC3_OTHER_simulated_data.RData")
