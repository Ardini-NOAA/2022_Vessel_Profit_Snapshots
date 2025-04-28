###############################################################################
#rm(list=ls()); graphics.off()
###############################################################################
# cc=function(nm='X',...){
#  arglist=list(...)
#  tmp=nm
#  if(length(arglist>1)){
#   print(paste('length dots = ',length(arglist)))
#   (tmp=paste(tmp,'[',arglist[[1]],sep=''))
#   for(j in 2:length(arglist)){
#    tmp=paste(tmp,arglist[[j]],sep=',')
#   }
#  tmp=paste(tmp,']',sep='')
#  } # end if(length(arglist>1))
#  tmp
# } # end function cc
###############################################################################
ccc=function(nm='X',index=NULL){
 tmp=nm
 if(!is.null(index)){
  #print(paste('length index = ',length(index)))
  (tmp=paste(tmp,'[',index[1],sep=''))
  if(length(index)>=2){
 	 for(j in 2:length(index)){
   tmp=paste(tmp,index[j],sep=',')
   } # end for
  } # end if	
  tmp=paste(tmp,']',sep='')
 } # end if(!is.null(index))
 tmp
} # end function ccc
###############################################################################
ulong=function(x,na.rm=F){
  if(na.rm==T) x=x[!is.na(x)]
  length(unique(x))
}
###############################################################################
# 
# ###############################################################################
# cc('X',1,2,3,4)
# cc('X')
# 
# ccc('Y')
# ccc('Y',1)
# ccc('Y',c(1,3,5))
# ###############################################################################




###############################################################################
createLP=function(sense='min'){
 list(LPsense=sense,nr=0,nc=0,obj = NULL,A=NULL,dirs=NULL,rhs=NULL,
      cnames=NULL,vartypes=NULL,xlower=NULL,xupper=NULL, 
      nvblocks=0,vblocks=list(NULL),SMM='A',CINDEX=F)
} # end function createLP
###############################################################################
# LP=createLP()
# LP
###############################################################################




###############################################################################
# debug on
#v='X';vindex=expand.grid(i=1);vtype='C';vlower=0.0;vupper=Inf
#
#v='X';vindex=expand.grid(i=1:4);vtype='C';vlower=0.0;vupper=Inf
#
#v='X';vindex=expand.grid(i=1:4,j=1:3);vtype='C';vlower=0.0;vupper=Inf
#vindex=doBy::orderBy(~i+j,data=vindex)
#
#v='X';vindex=expand.grid(i=1:4,j=1:3,k=1:2);vtype='C';vlower=0.0;vupper=Inf
#vindex=doBy::orderBy(~i+j+k,data=vindex)
#
###############################################################################
addvarLP=function(LP,v='X',vindex=expand.grid(i=1:5,j=1:3,k=1:2),
                   vtype='C',vlower=0.0,vupper=Inf){
###############################################################################
    (dims=apply(vindex,2,ulong))
    vnames=''
    VAR=array(' ',dim=apply(vindex,2,ulong))
    for(icell in 1:nrow(vindex)) vnames[icell]=ccc(v,vindex[icell,])
    VAR[as.matrix(cbind(vindex))] = vnames                                               
    (nc=length(vnames))
    LP$nc=LP$nc+nc
    LP$cnames=c(LP$cnames,vnames)
    LP$vartypes=c(LP$vartypes,rep(vtype,nc))
    LP$xlower=c(LP$xlower,rep(vlower,nc))
    LP$xupper=c(LP$xupper,rep(vupper,nc))
    LP$nvblocks=LP$nvblocks+1
    LP$vblocks[[LP$nvblocks]]=VAR
   	LP$obj=rep(0,LP$nc)
   	names(LP$obj)=LP$cnames
   	names(LP$xlower)=LP$cnames
    names(LP$xupper)=LP$cnames    
    LP <<- LP
    VAR
###############################################################################  
} # end function addvarLP
###############################################################################
#X=addvarLP(LP)
###############################################################################

