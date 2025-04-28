###############################################################################
#debug on run off
#LP=LPM;solver=solver
###############################################################################
LPSOLVE=function(LP,solver='clp',printobj=F){
  if(solver=='highs'|solver=='HiGHS')     tmpLPM=LP_highs(LP)  
  if(solver=='clp'|solver=='CLP')         tmpLPM=LP_clp(LP)
  if(solver=='glpk'|solver=='GLPK')       tmpLPM=LP_glpk(LP)
  if(solver=='lpsolve'|solver=='lpSolve') tmpLPM=LP_lpSolve(LP)
  if(solver=='rglpk'|solver=='Rglpk')     tmpLPM=LP_Rglpk(LP)  
  if(solver=='Rsym'|solver=='rsym')       tmpLPM=LP_Rsymphony(LP)
  if(solver=='lindo'|solver=='LINDO')     tmpLPM=LP_Lindo(LP)
  if(solver=='gurobi'|solver=='GUROBI')   tmpLPM=LP_gurobi(LP)
  if(printobj==T) print(paste(solver,'objval',tmpLPM$objval))
  tmpLPM
} # end function LPSOLVE
###############################################################################
my_seconds=function(){
    time1=proc.time()
    as.numeric(time1[3])
}
###############################################################################
source('A2SM-scripts.R')
###############################################################################



###############################################################################
# HiGHS LP solver function
LP_highs=function(LP){
###############################################################################
require(highs);require(Matrix)
#######################################
# Note: This code is set up to allow
# the user to measure the time used in
# different procedures. The following
#lines pull initial times from the system
#######################################
 time000=my_seconds()
 time00=my_seconds()
#######################################
#Convert LP matrix format (if necessary)
#######################################
 if(LP$SMM!='CRI') {
   LP=LP2LP(LP,SMM2='CRI',CINDEX2=F)
 }

 if(LP$CINDEX!=F)     LP=LP2LP(LP,SMM2='CRI',CINDEX2=F)
# Compute sparse matrix conversion time used
 SMtime=my_seconds()-time00
#######################################
# Construct lower and upper bounds on x
# If xbounds are not specified in LP
# we assume 0 <= x <= Inf
#######################################
 if('xlower' %in% names(LP) == F) LP$xlower = rep(0,length(obj))
 if('xupper' %in% names(LP) == F) LP$xupper = rep(Inf,length(obj))
#######################################
# Construct lower and upper bounds on Ax
# If rbounds are not specified in LP
# we construct them from the LP list's
# 'rhs' and 'dirs' vectors
#######################################
 vartypes=LP$vartypes
 if(length(which(vartypes=='B'))>0){
  pick=which(vartypes=='B')
  xlower=LP$xlower
  xupper=LP$xupper
  vartype[pick]='I'
  xlower=ifelse(xlower>0.5,1,0)
  xupper=ifelse(xlower>0.5,1,0)
  LP$vartypes=vartypes
  LP$xlower=xlower
  LP$xupper=xupper 
 } # end if(length(which(vartypes=='B'))>0)
 
 if('rlower' %in% names(LP) == F) {
  rlower     = LP$rhs
  rlower     = ifelse(LP$dirs=='<='|LP$dirs=='<',-Inf,rlower)
  rlower     = ifelse(LP$dirs=='>='|LP$dirs=='>',LP$rhs,rlower)
  LP$rlower  = rlower
 }

 if('rupper' %in% names(LP) == F) {
  rupper     = LP$rhs
  rupper     = ifelse(LP$dirs=='<='|LP$dirs=='<',LP$rhs,rupper)
  rupper     = ifelse(LP$dirs=='>='|LP$dirs=='>',Inf,rupper)
  LP$rupper  = rupper
 }

 (LPsense = ifelse(LP$LPsense!='max'|LP$LPsense!='MAX',F,T))

 ia=LP$ia
 ja=LP$ja
 ra=LP$ra

 ASM<- sparseMatrix(ia, ja, x = ra)

if(!exists('threads')) threads=1 
#print(paste('highs threads =',threads))
time0=my_seconds()

tmp <- highs_solve(L = LP$obj, lower = LP$xlower, upper = LP$xupper,
        A = ASM, lhs = LP$rlower, rhs = LP$rupper,types=LP$vartypes,
        maximum=LPsense,control=list(threads=threads))

proctime=my_seconds()-time0

 # retrieve the results
 (status=tmp$status)
 (objval=tmp$objective_value)
 (solution=tmp$primal_solution)



 TIME=my_seconds()-time000

 list(status=status,objval=objval,solution=solution,TIME=TIME,SMtime=SMtime,
  proctime=proctime,LP=LP)

} # end function LP_highs
###########################################################
###############################################################################
# clpAPI LP solver function
LP_clp=function(LP){
###############################################################################
require(clpAPI)
#######################################
# Note: This code is set up to allow
# the user to measure the time used in
# different procedures. The following
#lines pull initial times from the system
#######################################
 time000=my_seconds()
 time00=my_seconds()
#######################################
#Convert LP matrix format (if necessary)
#######################################
LP=LP2LP(LP,SMM2='CMO',CINDEX2=T)
#######################################
typecheck=c(grep('I',LP$vartypes),grep('B',LP$vartypes))
if(length(typecheck)>0) stop('clp cannot solve integer or binary LP problems') 
# Compute sparse matrix conversion time used
 SMtime=my_seconds()-time00
#######################################
# Construct lower and upper bounds on x
# If xbounds are not specified in LP
# we assume 0 <= x <= Inf
#######################################
 if('xlower' %in% names(LP) == F) LP$xlower = rep(0,length(obj))
 if('xupper' %in% names(LP) == F) LP$xupper = rep(Inf,length(obj))
#######################################
# Construct lower and upper bounds on Ax
# If rbounds are not specified in LP
# we construct them from the LP list's
# 'rhs' and 'dirs' vectors
#######################################

 if('rlower' %in% names(LP) == F) {
  rlower     = LP$rhs
  rlower     = ifelse(LP$dirs=='<='|LP$dirs=='<',-Inf,rlower)
  rlower     = ifelse(LP$dirs=='>='|LP$dirs=='>',LP$rhs,rlower)
  LP$rlower  = rlower
 }

 if('rupper' %in% names(LP) == F) {
  rupper     = LP$rhs
  rupper     = ifelse(LP$dirs=='<='|LP$dirs=='<',LP$rhs,rupper)
  rupper     = ifelse(LP$dirs=='>='|LP$dirs=='>',Inf,rupper)
  LP$rupper  = rupper
 }

 nr = LP$nr
 nc = LP$nc
 obj = LP$obj
 (sense = ifelse(LP$LPsense!='max'&LP$LPsense!='MAX',1,-1))

 ia=LP$ia
 ja=LP$ja
 ra=LP$ra

 xlower     = LP$xlower
 xupper     = LP$xupper
 rlower     = LP$rlower
 rupper     = LP$rupper
 dirs       = LP$dirs

lp <- initProbCLP()
 # supress output
 setLogLevelCLP(lp,0)
 # direction of optimization
 setObjDirCLP(lp, sense)  # 1 = min, -1 = max

 # load problem data
 loadProblemCLP(lp, ncols=nc, nrows=nr, ia=ia, ja=ja, ra=ra,
  lb=xlower, ub=xupper,obj_coef=obj,rlb=rlower,rub=rupper)

   scaleModelCLP(lp,3)

  # solve lp problem
 time0=my_seconds()
  solveInitialCLP(lp)
  #primalCLP(lp,ifValP=0)
  #dualCLP(lp,ifValP=0)
  proctime=my_seconds()-time0

 # retrieve the results
 (status=getSolStatusCLP(lp))
 (objval=getObjValCLP(lp))
 (solution=getColPrimCLP(lp))
 (dual=getRowDualCLP(lp))
 (rcost=getColDualCLP(lp))
 (lhs=getRowPrimCLP(lp))


 # remove problem
 delProbCLP(lp)

 TIME=my_seconds()-time000

 list(status=status,objval=objval,solution=solution,dual=dual,rcost=rcost,
      lhs=lhs,TIME=TIME,SMtime=SMtime,proctime=proctime,LP=LP)

} # end function LP_CLP
###########################################################


###########################################################
LP_glpk=function(LP){
###########################################################
 require(glpkAPI)
 time000=my_seconds()
 time00=my_seconds()

# Check and convert LP matrix format if necessary
 LP=LP2LP(LP,SMM2='CRI',CINDEX2=F)

 SMtime=my_seconds()-time00

 obj=LP$obj
 
 
 if('vartypes'  %in% names(LP) == F) LP$vartypes  = rep('C',length(obj))

 if('xlower' %in% names(LP) == F) LP$xlower = rep(0,length(obj))
 if('xupper' %in% names(LP) == F) LP$xupper = rep(Inf,length(obj))

 if('rlower' %in% names(LP) == F) {
 #if(is.null(LP$rlower)){
  rlower     = LP$rhs
  rlower     = ifelse(LP$dirs=='<='|LP$dirs=='<',-Inf,rlower)
  rlower     = ifelse(LP$dirs=='>='|LP$dirs=='>',LP$rhs,rlower)
  LP$rlower  = rlower
 }

 if('rupper' %in% names(LP) == F) {
 #if(is.null(LP$rupper)){
  rupper     = LP$rhs
  rupper     = ifelse(LP$dirs=='<='|LP$dirs=='<',LP$rhs,rupper)
  rupper     = ifelse(LP$dirs=='>='|LP$dirs=='>',Inf,rupper)
  LP$rupper  = rupper
 }

 (nr = LP$nr)
 (nc = LP$nc)
 obj = LP$obj
 (LPsense = LP$LPsense)

 ia=LP$ia
 ja=LP$ja
 ra=LP$ra
 (nnz=length(ra))

 vartypes   = LP$vartypes
 xlower     = LP$xlower
 xupper     = LP$xupper
 rlower     = LP$rlower
 rupper     = LP$rupper

 lp <- initProbGLPK()
 # direction of optimization
 if(LPsense=='min'|LPsense=='MIN') setObjDirGLPK(lp,GLP_MIN)
 if(LPsense=='max'|LPsense=='MAX') setObjDirGLPK(lp,GLP_MAX)

 # add rows and columns
 addRowsGLPK(lp, nr)
 addColsGLPK(lp, nc)

 setObjCoefsGLPK(lp, j=1:nc, obj_coef=obj)
# setColsBndsGLPK(lp, j=1:nc, lb=xlower, ub=xupper, type=vartypes)
 setColsBndsGLPK(lp, j=1:nc, lb=xlower, ub=xupper)
 setRowsBndsGLPK(lp, c(1:nr), rlower, rupper)

 # load constraint matrix
 loadMatrixGLPK(lp, nnz, ia, ja, ra)
 # suppress terminal output
 termOutGLPK(GLP_OFF)

 time0=my_seconds()
 # solve lp problem
 solveSimplexGLPK(lp)
 proctime=my_seconds()-time0


 # retrieve the results
 (status = getSolStatGLPK(lp))
 (objval = getObjValGLPK(lp))
 solution = getColsPrimGLPK(lp)
 dual=getRowsDualGLPK(lp)
 rcost=getColsDualGLPK(lp)
 lhs=getRowsPrimGLPK(lp) 
 
 # remove problem object
 delProbGLPK(lp)

 TIME=my_seconds()-time000

 list(status=status,objval=objval,solution=solution,
            dual=dual,rcost=rcost,lhs=lhs,
      TIME=TIME,SMtime=SMtime,proctime=proctime,LP=LP)

} # end function LP_glpk
###########################################################

###########################################################
LP_gurobi=function(LP,Threads=1,Cutoff=NULL,TimeLimit=36000,
  params=list()){
###########################################################
 require(gurobi)
 time000=my_seconds()
 model <- list()

if(LP$SMM!='A') {
 LP=LP2LP(LP,SMM2='A',CINDEX2=F)
}

 model$A          = LP$A
 model$obj        = LP$obj
 (model$modelsense = LP$LPsense)
 model$rhs        = LP$rhs
 model$lb         = LP$xlower
 model$ub         = LP$xupper
 GBsense          = LP$dirs
 GBsense=ifelse(GBsense=='<='|GBsense=='<','<',GBsense)
 GBsense=ifelse(GBsense=='>='|GBsense=='>','>',GBsense)
 GBsense=ifelse(GBsense=='=='|GBsense=='=','=',GBsense)
 #GBsense

 (model$sense      = GBsense)
 (model$vtype      = LP$vartypes)

 if(is.null(Cutoff)){
   Cutoff=ifelse(LP$LPsense=='min'|LP$LPsense=='MIN',1E30,-1E30)
 }


# params <- list(OutputFlag=0)
 params$OutputFlag=0
 params$Threads=Threads
 params$TimeLimit=TimeLimit

 time0=my_seconds()
 result <- gurobi(model, params)
 proctime=my_seconds()-time0

 (objval=result$objval)
 solution=result$x
 # Clear space
 rm(model, result, params)

 TIME=my_seconds()-time000


 list(objval=objval,solution=solution,
       TIME=TIME,proctime=proctime)

} # end function LP_Gurobi
###########################################################


###########################################################
LP_Lindo=function(LP){
###########################################################
  require(rLindo)
  time000=my_seconds()
  #Create LINDO enviroment object
  rEnv <- rLScreateEnv()

  #Create LINDO model object
  rModel <- rLScreateModel(rEnv)

  #Disable printing log
  rLSsetPrintLogNull(rModel)

  time00=my_seconds()
#  tmp=Lindo_A2SM(LP$A)
#  tmp=A2SM(LP$A,SMM='CMO',CINDEX=T)

  if(LP$SMM!='CMO'|LP$CINDEX!='T') {
   LP=LP2LP(LP,SMM2='CMO',CINDEX2=T)
  }

  (SMtime=my_seconds()-time00 )
  #Define the model data

  (nr=as.integer(LP$nr))
  (nc=as.integer(LP$nc))
  ia=as.integer(LP$ia)
  ja=as.integer(LP$ja)
  ra=LP$ra

  (nNZ=length(ra))


  (LPsense <- ifelse(LP$LPsense=='max'|LP$LPsense=='MAX',LS_MAX,LS_MIN))

  obj=LP$obj

  Ltypes=LP$dirs
  Ltypes=ifelse(Ltypes=='<='|Ltypes=='<','L',Ltypes)
  Ltypes=ifelse(Ltypes=='>='|Ltypes=='>','G',Ltypes)
  Ltypes=ifelse(Ltypes=='=='|Ltypes=='=','E',Ltypes)
  (Ltype=paste(Ltypes,sep='',collapse=''))


  rhs=LP$rhs

  pdLower =  LP$xlower
  pdUpper =  LP$xupper

  (vartypes=LP$vartypes)
  (Lvartypes=paste(vartypes,sep='',collapse=''))


  rLSloadLPData(model=rModel,nCons=nr,nVars=nc,nObjSense=LPsense,
    dObjConst=0,padC=obj,padB=rhs,pszConTypes=Ltype,nAnnz=nNZ,
    paiAcols=ja,panAcols = NULL,padAcoef=ra,paiArows=ia,
    padL=pdLower,padU=pdUpper)

  #rLSgetErrorMessage(rEnv,2053)

  time0=my_seconds()
  #Solve the model
  rLSoptimize(rModel,LS_METHOD_FREE)
  proctime=my_seconds()-time0


  #Get solution information
  (status=rLSgetDInfo(rModel,LS_DINFO_POBJ)$ErrorCode)
  (objval=rLSgetDInfo(rModel,LS_DINFO_POBJ)$pdResult)
  solution=rLSgetPrimalSolution(rModel)$padPrimal


  #Delete the model and environment
  rLSdeleteModel(rModel)
  rLSdeleteEnv(rEnv)

  time1=my_seconds()
  TIME=my_seconds()-time000

 list(status=status,objval=objval,solution=solution,
      TIME=TIME,SMtime=SMtime,proctime=proctime,LP=LP)
} # end function LP_LINDO
###########################################################


###########################################################
LP_lpSolve=function(LP){
###########################################################
require(lpSolveAPI)
time000=my_seconds()

if(LP$SMM!='CMO') LP=LP2LP(LP,SMM2='CMO',CINDEX2=F)
if(LP$CINDEX!=F)  LP=LP2LP(LP,SMM2='CMO',CINDEX2=F)

(LPsense = LP$LPsense)

nr=LP$nr;nc=LP$nc;obj=LP$obj;ia=LP$ia;ja=LP$ja;ra=LP$ra
dirs=LP$dirs; rhs=LP$rhs
####################################################################
# Construct lower and upper bounds on x
# If xbounds are not specified in LP
# we assume 0 <= x <= Inf
####################################################################
if('xlower' %in% names(LP) == F) LP$xlower = rep(0,length(obj))
if('xupper' %in% names(LP) == F) LP$xupper = rep(Inf,length(obj))
xlower=LP$xlower;xupper=LP$xupper
####################################################################
 if('vartypes'  %in% names(LP) == F) LP$vartypes  = rep('C',length(obj))
 vartypes=LP$vartypes
 vartypes[vartypes=='C']='real'
 vartypes[vartypes=='B']='binary'
 vartypes[vartypes=='I']='integer'
####################################################################
 dirs[dirs=='==']='='
 dirs[dirs=='<']='<='
 dirs[dirs=='>']='>='
 #################################################################### 
 lp=make.lp(nrow=nr,ncol=nc)
 lp.control(lp,sense=LPsense)

 j=1
 for(j in 1:nc){
  (jdxs=(ja[j]):(ja[j+1]-1))
  set.column(lp,j,ra[jdxs],indices=ia[jdxs])  
 }
 
 set.objfn(lp,obj)
 set.constr.type(lp,dirs)
 set.rhs(lp,rhs)
 set.bounds(lp, lower = xlower, upper = xupper)

 for(j in 1:nc) set.type(lp,j,vartypes[j])

 # solve lp problem
 time0=my_seconds()
 (status=solve(lp))
 (proctime=my_seconds() - time0)

(objval=get.objective(lp))
(solution=get.variables(lp))
(duals=get.dual.solution(lp))
(LHS=get.constraints(lp))

# write.lp(lp,file='LP_MPS_DATA.mps',type='freemps',use.names=c(F,F))
 
 delete.lp(lp)
 
 TIME=my_seconds()-time000

 list(status=status,objval=objval,solution=solution,
   duals=duals,LHS=LHS,TIME=TIME,proctime=proctime)

} # end function LP_lpSolve
###########################################################


###########################################################
LP_Rglpk=function(LP){
###########################################################
 require(Rglpk); require(Matrix)
 time000=my_seconds()
 

 LP=LP2LP(LP,SMM2='CRI',CINDEX2=F)
 ASM = sparseMatrix(LP$ia, LP$ja, x = LP$ra)

 (LPsense = ifelse(LP$LPsense!='max'|LP$LPsense!='MAX',F,T)) 
 
 nc=ncol(LP$A)

 if('xlower' %in% names(LP) == F) LP$xlower = rep(0,length(LP$obj))
 if('xupper' %in% names(LP) == F) LP$xupper = rep(Inf,length(LP$obj))
 

 bounds=list(lower=list(ind=1:LP$nc,val=LP$xlower),
             upper=list(ind=1:LP$nc,val=LP$xupper))

 
 time0=my_seconds()

 tmpLP= Rglpk_solve_LP(LP$obj, ASM, LP$dirs, LP$rhs, bounds, 
 	         LP$vartypes, max=LPsense) 
 

 proctime=my_seconds()-time0


 (status=as.numeric(tmpLP$status))
 (objval=tmpLP$optimum)
 solution=tmpLP$solution
 rcost=tmpLP$solution_dual
 dual=tmpLP$auxiliary$dual
 lhs=tmpLP$auxiliary$primal
 TIME=my_seconds()-time0
 list(status=status,objval=objval,solution=solution,dual=dual,rcost=rcost,
      lhs=lhs,TIME=TIME,proctime=proctime)
} # end function LP_Rglpk
###########################################################

###########################################################
LP_Rsymphony=function(LP){
###########################################################
 require(Rsymphony)
 time000=my_seconds()
 (LPsense = ifelse(LP$LPsense!='max'|LP$LPsense!='MAX',F,T))
 
 if(LP$SMM!='A'){
  LP=LP2LP(LP,SMM2='A',CINDEX2=F)
 }

 nc=ncol(LP$A)

 if('xlower' %in% names(LP) == F) LP$xlower = rep(0,length(obj))
 if('xupper' %in% names(LP) == F) LP$xupper = rep(Inf,length(obj))
 
 lower=list(ind=1:nc,val=LP$xlower)
 upper=list(ind=1:nc,val=LP$xupper)
 bounds=list(lower=lower,upper=upper)


 
 time0=my_seconds()

 tmpLP=Rsymphony_solve_LP(obj=LP$obj, mat=LP$A, dir=LP$dirs,
       rhs=LP$rhs,types=LP$vartypes,max = LPsense,bounds = bounds) 
 
 proctime=my_seconds()-time0


 (status=as.numeric(tmpLP$status))
 (objval=tmpLP$objval)
 solution=tmpLP$solution
 TIME=my_seconds()-time0
 list(status=status,objval=objval,solution=solution,
      TIME=TIME,proctime=proctime)
} # end function LP_Rsymphony
###########################################################


###########################################################
MILP_Gurobi=function(LPB,Threads=1,Cutoff=NULL,
    TimeLimit=36000,params=list()){
###########################################################
 require(gurobi)
 time000=my_seconds()
 model <- list()

 if(LP$SMM!='A') {
  LP=LP2LP(LP,SMM2='A',CINDEX2=F)
 }

 model$A          = LPB$A
 model$obj        = LPB$obj
 (model$modelsense = LPB$LPsense)
 model$rhs        = LPB$rhs
 model$lb         = LPB$xlower
 model$ub         = LPB$xupper
 GBsense          = LPB$dirs
 GBsense=ifelse(GBsense=='<='|GBsense=='<','<',GBsense)
 GBsense=ifelse(GBsense=='>='|GBsense=='>','>',GBsense)
 GBsense=ifelse(GBsense=='=='|GBsense=='=','=',GBsense)
 #GBsense

 (model$sense      = GBsense)
 (model$vtype      = LPB$vartypes)

 if(is.null(Cutoff)){
   Cutoff=ifelse(LPB$LPsense=='min'|LPB$LPsense=='MIN',1E30,-1E30)
 }


# params <- list(OutputFlag=0)
 params$OutputFlag=0
 params$Threads=Threads
 params$TimeLimit=TimeLimit

 time0=my_seconds()
 result <- gurobi(model, params)
 proctime=my_seconds()-time0


 (objval=result$objval)
 solution=result$x
 # Clear space
 rm(model, result, params)
 TIME=my_seconds()-time000
 list(objval=objval,solution=solution,
      TIME=TIME,proctime=proctime)

} # end function MILP_Gurobi
###########################################################


###########################################################
MILP_Lindo=function(LPB,TimeLimit=1200){
###########################################################
 require(rLindo)
 time000=my_seconds()
# tmp=Lindo_A2SM(LPB$A)
# tmp=A2SM(LPB$A,SMM='CMO',CINDEX = T)

 if(LPB$SMM!='CMO'|LPB$CINDEX!='T'){
  LPB=LP2LP(LPB,SMM2='CMO',CINDEX2=T)
 }

 SMtime=my_seconds()-time000

 #Create LINDO environment object
 rEnv <- rLScreateEnv()

 #Create LINDO model object
 rModel <- rLScreateModel(rEnv)

 #Disable printing log
 rLSsetPrintLogNull(rModel)


  #Define the model data

  (nr=as.integer(LPB$nr))
  (nc=as.integer(LPB$nc))
  ia=as.integer(LPB$ia)
  ja=as.integer(LPB$ja)
  ra=LPB$ra

 (nNZ=length(ra))
  (LPsense=LPB$LPsense)

 (LPsense <- ifelse(LPB$LPsense=='max'|LPB$LPsense=='MAX',LS_MAX,LS_MIN))

 obj=LPB$obj

 Ltypes=LPB$dirs
 Ltypes=ifelse(Ltypes=='<='|Ltypes=='<','L',Ltypes)
 Ltypes=ifelse(Ltypes=='>='|Ltypes=='>','G',Ltypes)
 Ltypes=ifelse(Ltypes=='=='|Ltypes=='=','E',Ltypes)
 (Ltype=paste(Ltypes,sep='',collapse=''))

 rhs=LPB$rhs

 pdLower =  LPB$xlower
 pdUpper =  LPB$xupper

 (vartypes=LPB$vartypes)
 (Lvartypes=paste(vartypes,sep='',collapse=''))

 rLSloadLPData(model=rModel,nCons=nr,nVars=nc,nObjSense=LPsense,
     dObjConst=0,padC=obj,padB=rhs,pszConTypes=Ltype,nAnnz=nNZ,
     paiAcols=ja,panAcols = NULL,padAcoef=ra,paiArows=ia,
     padL=pdLower,padU=pdUpper)

 #Load data to the model
 rLSloadVarType(rModel,Lvartypes)

 # load time limits
 rLSsetModelDouParameter(model=rModel,nParameter=LS_DPARAM_MIP_TIMLIM,dValue=TimeLimit)

 #rLSgetErrorMessage(rEnv,2053)


 #Solve the model
 time0=my_seconds()
 rLSsolveMIP(rModel)
 proctime=my_seconds()-time0


 #Get solution information
 (status=rLSgetIInfo(rModel,LS_IINFO_MIP_STATUS)$pnResult) # integer type'd info
 (objval=rLSgetDInfo(rModel,LS_DINFO_MIP_OBJ)$pdResult) # double type'd info

  xvals=rLSgetMIPPrimalSolution(rModel)$padPrimal


 #Delete the model and environment
 rLSdeleteModel(rModel)
 rLSdeleteEnv(rEnv)

TIME=my_seconds()-time000
list(status=status,objval=objval,xvals=xvals,
     TIME=TIME,SMtime=SMtime,proctime=proctime,LPB=LPB)

} # end function MILP_LINDO
###########################################################




###############################################################################
# debug on
#LP1=LP1;SMM2='CMO';CINDEX2='T'
#LP1=LP2;SMM2='A';CINDEX2='F'
#LP1=LP2;SMM2='CRI';CINDEX2='F'
###############################################################################
LP2LP=function(LP1,SMM2='CMO',CINDEX2=F){
###############################################################################
  LP2=LP1
  (SMM1=LP1$SMM)
  (CINDEX1=LP1$CINDEX)


  if(SMM1=='A'& SMM2!='A') {
    A2 = A2SM(LP1$A,SMM=SMM2,CINDEX=CINDEX2)
  }


  if(SMM1!='A'& SMM2!='A') {
    A1=list(nnz=LP1$nnz,nr = LP1$nr,nc = LP1$nc,ia = LP1$ia,ja = LP1$ja,
        ra = LP1$ra,SMM = SMM1,CINDEX = CINDEX1)
    A2 = SM2SM(A1,SMM2=SMM2,CINDEX2=CINDEX2)
  }


  if(SMM1!='A' & SMM2=='A'){
   A = SM2A(LP1$A)
  }


 if(SMM2 != 'A'){
  LP2$nnz=A2$nnz
  LP2$nr=A2$nr
  LP2$nc=A2$nc
  LP2$ia=A2$ia
  LP2$ja=A2$ja
  LP2$ra=A2$ra
  LP2$SMM=SMM2
  LP2$CINDEX=A2$CINDEX
  LP2=LP2[names(LP2) %in% c('A','ASM')==FALSE]
 }

 if(SMM2 == 'A'){
  LP2$A=A
  LP2$SMM=SMM2
  LP2=LP2[names(LP2)%in%c('nnz','nr','nc','ia','ja','ra','ASM','CINDEX')==FALSE]
 }

  LP2
} # end function LP2LP
###############################################################################

###############################################################################
#' Write.LPC: Writes LP to .csv file format
#'
#' @param LP      A sparse LP object
#' @param fname   .csv file name (include path if desired)
#' @return        = .csv file (can be opened/solved with spreadsheet solver)
###############################################################################
# debugging
#fname='LPC.csv'
###############################################################################
Write.LP=function(LP,fname='LP.csv'){
###############################################################################
  if(length(grep('LP',search()))>0) detach(LP)
  if(exists('A')) rm(A,nr,nc)
 attach(LP)
  A = SM2A(LP)
  (nc=ncol(A))
  (nr=nrow(A))
  M=A

  if(exists('rhs'))    M=cbind(M,'|','lhsA','|',LP$dirs,'|',LP$rhs,'|')
  if(exists('rlower')) M=cbind(M,LP$rlower,'|',LP$rupper,'|')

  (ncolM=ncol(M))

  (tmp=c(paste('x',1:nc,sep=''),'|','LHS','|'))
  if(exists('xnames')){
    (tmp=c(xnames,'|','LHS','|'))
  } 
  

  if(exists('rhs'))    (tmp=c(tmp,'rest','|','rhs','|'))
  if(exists('rlower')) (tmp=c(tmp,'rlower','|','rupper','|'))

  M=rbind(rep('========',ncol(M)),M)
  M=rbind(c(tmp,rep('',ncolM-length(tmp))),M)
  M=rbind(rep('========',ncol(M)),M)

  M=rbind(M,rep('========',ncol(M)))

  (tmp=c(LP$obj,'|','lhsOBJ','=',paste(LP$LPsense)))
  M=rbind(M,c(tmp,rep('',ncolM-length(tmp))))

  M=rbind(M,rep('========',ncol(M)))

  (tmp=c('nrows =',nr,'ncols =',nc))
  M=rbind(c(tmp,rep('',ncolM-length(tmp))),M)


  (tmp=c(LP$xlower,'|','xlower'))
  M=rbind(M,c(tmp,rep('',ncolM-length(tmp))))
  tmp=c(LP$xupper,'|','xupper')
  M=rbind(M,c(tmp,rep('',ncolM-length(tmp))))

  M=rbind(M,rep('========',ncol(M)))

  (tmp=c(0*LP$obj,'|','xvalues'))
  M=rbind(M,c(tmp,rep('',ncolM-length(tmp))))

  M=rbind(M,rep('========',ncol(M)))
  M=rbind(rep('========',ncol(M)),M)

  (tmp=c(vartypes,'|','vartypes'))
  M=rbind(M,c(tmp,rep('',ncolM-length(tmp))))
  M=rbind(rep('========',ncol(M)),M)

 detach(LP)

 rm(A,nr,nc)

 write.table(M,file=fname,row.names=F,col.names=F,sep=',')
###############################################################################
} # end function Write.LP
###############################################################################

