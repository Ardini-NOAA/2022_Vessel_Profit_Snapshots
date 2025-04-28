###############################################################################
rm(list=ls())
###############################################################################
#library(dplyr)
#library(clpAPI)
library(glpkAPI)
library(lpSolveAPI)
library(Rglpk)
library(Rsymphony)
library(Matrix)
library(readr)
library(here)
library(psych)
library(asbio)
library(AER)
library(bea.R)
library(readxl)
library(plyr)
source('V:/Ardini_Cost_Survey/SAS_Cost_Survey/Profitability_Profiles/Calculate_Profit/Depreciation/Vessel_Values/beakey.R')  #BEA Key for building deflator
###############################################################################
source('V:/Ardini_Cost_Survey/SAS_Cost_Survey/Profitability_Profiles/Calculate_Profit/Depreciation/Vessel_Values/LP-Modeling-Joe.R')
source('V:/Ardini_Cost_Survey/SAS_Cost_Survey/Profitability_Profiles/Calculate_Profit/Depreciation/Vessel_Values/LP_CODES.R')
###############################################################################
###############################################################################
#Functions to solve LP model with various solvers
###############################################################################
# clpAPI solver
# status code=0 means the LP model solved
LP_clpAPI<-function(LP){
  SOLC=LPSOLVE(LP,solver='clp',printobj=F)
  status<-SOLC$status
  objret<-SOLC$objval
  sol<-SOLC$solution
  
  solList <- list("status" = status, "obj" = objret, "sol"=sol)
  
  return(solList)
}
###############################################################################
LP_Symphony<-function(LP){
  # Rsymphony solver
  #status code=0 means the LP model solved
  
  SOLS=LPSOLVE(LP,solver='Rsym',printobj=F)
  status<-SOLS$status
  objret<- SOLS$objval
  sol=SOLS$solution
  
  solList <- list("status" = status, "obj" = objret, "sol"=sol)
  
  return(solList)
}
##############################################################################
LP_High<-function(LP){
  # highs solver
  #status code=7 means the LP Model Solved
  
  SOLH<-LPSOLVE(LP,solver='highs',printobj=F)
  status<-SOLH$status
  objret<-SOLH$objval
  sol<-SOLH$solution
  
  solList <- list("status" = status, "obj" = objret, "sol"=sol)
  
  return(solList)
}
#############################################################################
LP_RGLPK<-function(LP){
  # Rglpk solver method 4
  # status code=0 means the model solved
  
  SOLRG4=LPSOLVE(LP,solver='Rglpk',printobj=F)
  status<-SOLRG4$status
  objret<-SOLRG4$objval
  sol<-   SOLRG4$solution
  
  solList <- list("status" = status, "obj" = objret, "sol"=sol)
  
  return(solList)
}
###############################################################################
LP_LPSOLVE<-function(LP){
  # lpsolveAPI solver
  # status code=0 means the model solved
  # NOTE: Ignore any lpSolveAPI "lprec" error messages
  
  SOL_LS=LPSOLVE(LP,solver='lpsolve',printobj=F)
  status<-SOL_LS$status
  objret<-SOL_LS$objval
  sol<-SOL_LS$solution
  
  solList <- list("status" = status, "obj" = objret, "sol"=sol)
  
  return(solList)
}
###############################################################################
#End of Solver Functions
#Can find additional solver calls commented out at bottom of program
############################################################################
#Read in data 
load("V:/Ardini_Cost_Survey/SAS_Cost_Survey/Profitability_Profiles/Calculate_Profit/Depreciation/Vessel_Values/data/TOC1_simulated_data_v2.RData")
df1<-YSIM
###############################################################################
#Transform data for use in models
###############################################################################
V1 = as.matrix(df1[,"value"])      # 
Z1 = as.matrix(df1[,c("age2","vhp","len","tons")])  # inputs
V=log(V1)
Z=log(Z1)
(N=ncol(Z))
(K=nrow(df1))                            #Number of observations
################################################################################
###############################################################################
#Set-up the DEA model This step sets up the A 
#matrix, rhs values and objective function to use in DEA Model
###############################################################################
LP=createLP(sense='min')
# define variables. Variable labels for LP model
DL=addvarLP(LP,v='DL',vindex=expand.grid(k=1:K),vtype='C')
a0=addvarLP(LP,v='a0',vindex=expand.grid(1),vtype='C',vlower=-Inf)
a1=addvarLP(LP,v='a1',vindex=expand.grid(1),vtype='C',vlower=-Inf)
a11=addvarLP(LP,v='a11',vindex=expand.grid(1),vtype='C',vlower=-Inf)
B=addvarLP(LP,v='B',vindex=expand.grid(n=1:N),vtype='C',vlower=-Inf)
BB=addvarLP(LP,v='BB',vindex=expand.grid(n=1:N,np=1:N),vtype='C',vlower=-Inf)
G=addvarLP(LP,v='G',vindex=expand.grid(n=1:N),vtype='C',vlower=-Inf)
##########################################################################
# set obj
obj=LP$obj
for(k in 1:K) obj[DL[k]] = 1.0

# add constraints
# number of rows
(nr=(K+N+2)+(N^2-N)/2 + N*K)
A=matrix(0,nr,LP$nc); 
colnames(A)=LP$cnames
dirs=LP$dirs; 
rhs=LP$rhs
row=0
# restriction (2)
row=row+1; 
for(n in 1:N) A[row,B[n]]=1 ;dirs[row]='=='; rhs[row]=1
# restriction (3)
for(n in 1:N){
	row=row+1;for(np in 1:N) A[row,BB[n,np]]=1; dirs[row]='=='; rhs[row]=0
}  
# restriction (4)
row=row+1; for(n in 1:N) A[row,G[n]]=1; dirs[row]='=='; rhs[row]=0
# restrictions (5)
for(n in 1:N){for(np in 1:N){
 if(n<np){  
  row=row+1 
  A[row,BB[n,np]]=1.0; A[row,BB[np,n]]= -1.0; dirs[row]='==';rhs[row]=0
 } # end if  
}} # end double loops
# restrictions (6)
for(n in 1:N){for(k in 1:K){
 row=row+1; 
 A[row,B[n]]=1; for(np in 1:N) A[row,BB[n,np]]=Z[k,np]; A[row,G[n]]=V[k]
 dirs[row]='>='; rhs[row]=0.0
}}
# restrictions to compute DL
for(k in 1:K){
 row=row+1
 A[row,a0]=1;A[row,a1]=V[k];A[row,a11]=0.5*V[k]^2
 for(n in 1:N) A[row,B[n]]=Z[k,n]
 for(n in 1:N) for(np in 1:N) A[row,BB[n,np]]=0.5*Z[k,n]*Z[k,np]
 for(n in 1:N) A[row,G[n]]=Z[k,n]*V[k]
 A[row,DL[k]]=(-1.0)
 dirs[row]='=='; rhs[row]=0
} 
###############################################################################
row;
nr;
LP$nc
###############################################################################
# put values in LP list object
LP$nr=nr;LP$obj=obj;LP$A=A; LP$dirs=dirs; LP$rhs=rhs
###############################################################################
length(grep('LP',search()))
if(length(grep('LP',search()))>0) detach(LP)
nrow(df1)
if(nrow(df1)<=10) Write.LP(LP)
###############################################################################
#Solve model and get results
#Use Symphony solver
###############################################################################
model1<-LP_Symphony(LP)

(status<-model1$status) #Status of zero means the model solved
(obj<-model1$obj)
sol<-model1$sol         #Put model solution into sol
##############################################################################
#Set up shadow prices for each attribute and observation
###############################################################################
Beta<-matrix(0,1,N)
for(n in 1:N){
  Beta[1,n]=sol[(K+3+n)]
}

BB<-matrix(0,1,(N^2))

for(n in 1:(N^2)){
  BB[1,n]=sol[K+3+N+n]
}

BB2<-matrix(BB, N, byrow=TRUE)
Gamma<-matrix(0,1,N)

for(n in 1:N){
 Gamma[1,n]=sol[(K+3+N+(N^2)+n)]
}

BB2
Gamma
Beta
val1<-matrix(0,K,N)
#Set up matrix to hold vessel value divided by input values
for (n in 1:N){
  val1[,n]<-V1/Z1[,n]
}
 
sprice_inp<-matrix(0,K,N)
#####################################################################################################
#CAUTION: The equation in the loop below needs to be changed if the number of inputs is different than 4
for(k in 1:K){
 for(n in 1:N){
  sprice_inp[k,n]=Beta[n]+(BB2[n,1]*Z[k,1])+(BB2[n,2]*Z[k,2])+(BB2[n,3]*Z[k,3])+ (BB2[n,4]*Z[k,4])+Gamma[n]*V[k]  
 }
}
#######################################################################################################
#Final shadow prices for each vessel and input
val2<-round(val1*sprice_inp,3) 
val2[val2==0]<-NA

# R program to find the confidence interval

inp_names<-colnames(Z1)
colnames(val2)<-inp_names

summary(val2)


evalue<-t(as.matrix(winsor.mean(val2, trim=0.025, na.rm=TRUE)))
evalue<-as.data.frame(evalue)
evalue[is.na(evalue)]<-0

#######################################################################################
#These next steps are to plot out observed values and estimated values 
#Extract Data and Place in Dataframes
df1$est_value=(df1$age2*evalue$age2)+(df1$vhp*evalue$vhp)+(df1$len*evalue$len)
                   +(df1$tons*evalue$tons)
#Now read in survey Values
ves_val<- read_excel("V:/Ardini_Cost_Survey/SAS_Cost_Survey/Profitability_Profiles/Calculate_Profit/Depreciation/survey_response_2015_2022_char_gear.xlsx",sheet="Data")
colnames(ves_val)[3]<-"YEAR"
ves_val<-subset(ves_val, YEAR==2015 | YEAR==2022)
################################################################################
#Calculate GDP Deflator from BEA Data
SpecList<-list('UserID'= beaKey,
               'Method' = 'GetData',
               'datasetname' = 'NIPA',
               'TableName' = 'T10109',
               'RowNumber' = '10',
               'Frequency'='A',
               'Year'= '2015,2022'
)
GDP<-beaGet(SpecList, asWide=FALSE)  #Pull from BEA using API key
GDP<-GDP[(LineNumber==1),c("TimePeriod","DataValue")]
colnames(GDP)<-c("YEAR","GDP")
GDP$GDPD=GDP$GDP/100
GDP[,]<-lapply(GDP, function(x) type.convert(as.character(x), as.is=TRUE))
gdpd2022=subset(GDP, YEAR==2022)
GDP$GDPD2022=GDP$GDPD/gdpd2022$GDPD
GDP2<-GDP[,c(1,4)]
###############################################################################
#Convert vessel value survey data to 2022 values
ves_val<-join(ves_val,GDP2, by="YEAR", type="inner")
ves_val<-subset(ves_val, vessel_value>0)
ves_val$vessel_value=as.numeric(ves_val$vessel_value)
ves_val$val2022=ves_val$vessel_val/ves_val$GDPD2022
################################################################################
ves_val$age2=100-(ves_val$vessel_age)
ves_val<-subset(ves_val, toc==1 | toc==2 | toc==3)
ves_val$primary_gear=ifelse(is.na(ves_val$primary_comm_gear),"Other",ves_val$primary_comm_gear)
#############################################################################################################
#Subset into different TOC (Hull type 1=Wood, 2=Fiberglass, 3=Steel) and Gear Types
TOC1<-subset(ves_val, toc==1)

TOC3<-subset(ves_val, toc==3)

TOC3_Dredge<-subset(TOC3, primary_gear=="Dredge_Scallop")
TOC3_Dredge<-subset(TOC3_Dredge, val2022<2000000)  #Eliminate outliers

TOC3_Trawl<-subset(TOC3, primary_gear=="Trawl")
TOC3_Trawl<-subset(TOC3_Trawl, val2022<7000000)

TOC3_Other<-subset(TOC3, primary_gear != "Trawl" & primary_gear != "Dredge_Scallop")
TOC3_Other<-subset(TOC3_Other, val2022<4000000)
##################################################################################################
TOC2<-subset(ves_val, toc==2)

TOC2_HG<-subset(TOC2, primary_gear=="Handgear")
TOC2_HG<-subset(TOC2_HG, val2022<600000)

TOC2_POT<-subset(TOC2, primary_gear=="Pot/Trap")
TOC2_POT<-subset(TOC2_POT, val2022<1050000)

TOC2_TRAWL<-subset(TOC2, primary_gear=="Trawl")

TOC2_OTHER<-subset(TOC2, primary_gear!="Trawl" & primary_gear!="Pot/Trap"
                   & primary_gear!="Handgear")


#The code below is to compare estimated values with actual values and simulated Values
TOC_0<-TOC2_POT  #This has to be changed for each data set. Set to TOC1, TOC_TRAWL for example
 
dens_o <- density(df1$value)
dens_e <- density(df1$est_value)
dens_s <-density(TOC_0$val2022)
TOC_0$est_value=(TOC_0$age2*evalue$age2)+(TOC_0$vhp*evalue$vhp)+(TOC_0$len*evalue$len)
     +(TOC_0$gtons*evalue$tons)
# 
dens_se<-density(TOC_0$est_value)
 
xr <- range(c(dens_o$x, dens_e$x))
yr <- range(c(dens_o$y, dens_e$y))
# 
#Execute all four lines below at once
X11()
plot(dens_s, xlim=xr, ylim=yr,main="TOC1 Values") #Black Line original values
lines(dens_se, lty=4,col=5)  #estimated value original data, dashed line light blue
lines(dens_e, lty=3,col=2)  #estimated value simulated data, dashed line red
lines(dens_o, lty=2, col=1) #simulated values, black dashed line
