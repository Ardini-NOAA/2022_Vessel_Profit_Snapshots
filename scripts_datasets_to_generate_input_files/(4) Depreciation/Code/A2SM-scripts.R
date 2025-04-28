###############################################################################
###############################################################################
#' A2SM: Convert a matrix A to sparse matrix form
#'
#' @param A        The matrix A.
#' @param SMM      Sparse matrix method with SMM='A',CRI','RCI','CMO', or 'RMO'
#' @param CINDEX   T = Use zero indexing or F = Use one indexing.
#' @param eps      Value to use in non-zero test.
#' @param NAflag   Number to uses as NA flag
#'
#' @return A list object containing the matrix components:
#' @return nnz     = the number of non-zero elements in ra.
#' @return nr      = the number of rows in matrix A.
#' @return nc      = the number of columns in matrix A.
#' @return ia      = the row index.
#' @return ja      = the column index.
#' @return ra      = the non-zero coefficients in A
#' @return rnames  = matrix row names -- may be  ""
#' @return cnames  = matrix column names -- may be ""
#' @return SMM     = the sparse matrix type.
#' @return CINDEX with T = Use zero indexing or F = Use one indexing.
#' @examples
#' \dontrun{
#' (A = matrix(c(1,0,0,2,0,3,0,0,0,5,0,6),3,4,byrow=T))
#' (SM1 = A2SM(A,SMM='CMO',CINDEX=T))
#' (SM2=SM2SM(SM1,SMM='CRI',CINDEX=F))
#' SM2A(SM1)
#' SM2A(SM2)
#'
#' (A = matrix(c(1,0,0,2,0,3,0,0,0,5,0,6),3,4,byrow=T)); A[2,3]=NA; A
#' (ASM = A2SM(A,SMM='CMO',CINDEX=T)); SM2A(ASM)
#' (ASM2 = SM2SM(ASM,SMM2='CRI',CINDEX2=T)); A ; (A2 = SM2A(SM2))
#'
#' #clpAPI  documentation example
#' nr=5
#' nc=8
#' ra=c(3.0,5.6,1.0,2.0,1.1,1.0,-2.0,2.8,-1.0,1.0,1.0,-1.2,-1.0,1.9)
#' ia=c(0,4,0,1,1,2,0,3,0,4,2,3,0,4)
#' ja=c(0,2,4,6,8,10,11,12,14)
#' SMM='CMO'
#' CINDEX=T
#' ASM=list(nr=nr,nc=nc,ra=ra,ia=ia,ja=ja,SMM=SMM,CINDEX=CINDEX)
#' (A=SM2A(ASM))
#' }
###############################################################################
A2SM = function(A, SMM = 'CRI',CINDEX = F,
                eps = .Machine$double.eps,NA_flag=(-1E6)){
##############################################################################
  A[is.na(A)] = NA_flag
##############################################################################
  (nr = nrow(A))
  (nc = ncol(A))
  if(is.null(rownames(A))) rownames(A)=rep('',nr)
  if(is.null(colnames(A))) colnames(A)=rep('',nc)
  (rnames=rownames(A))
  (cnames=colnames(A))
##############################################################################
# jatwood algorithm 10/29/2022
# initially uses SMM RCI
  tmp=list(NULL); eps=1E-9
  j=1
  for(j in 1:nc){
   a=as.vector(A[,j])
   ipick=which(abs(a)>eps)
   ra=a[ipick]
   ia=ipick
   ja=rep(j,length(ipick))
   as.data.frame(cbind(ia,ja,ra))
   tmp[[j]]=as.data.frame(cbind(ia,ja,ra))
  } # end loop for(j in 1:nc)
  
  tmp=bind_rows(tmp)
  
  (nnz=nrow(tmp))

  ASM=list(nnz=nnz,nr=nr,nc=nc,ia=tmp$ia,ja=tmp$ja,ra=tmp$ra,
         rnames=rnames,cnames=cnames,SMM='CRI',CINDEX=F)
  
  rm(tmp)
  
# convert to user specified sparse matrix format     
  ASM=SM2SM(ASM,SMM2=SMM,CINDEX2=CINDEX)

  return(ASM)
} # end function A2SM
###############################################################################


###############################################################################
#' SM2A: Convert a sparse matrix object ASM into a matrix A
#'
#' ASM is a sparse matrix object containing:
#' nr = number of rows, nc = number of columns, ra = nonzero coeff,
#' ia = row indices, ja = col indices,
#' rnames = row names,cnames = column names # Note these may be missing or ""
#' SMM = space matrix method, and CINDEX = use 0-indexing
#' @param ASM  A sparse matrix object
#' @param NAflag   Number to uses as NA flag
#'
#' @return The matrix A.
#' @examples
#' \dontrun{
#' (A = matrix(c(1,0,0,2,0,3,0,0,0,5,0,6),3,4,byrow=T))
#' (SM1 = A2SM(A,SMM='CMO',CINDEX=T))
#' (SM2=SM2SM(SM1,SMM='CRI',CINDEX=F))
#' SM2A(SM1)
#' SM2A(SM2)
#'
#' (A = matrix(c(1,0,0,2,0,3,0,0,0,5,0,6),3,4,byrow=T)); A[2,3]=NA; A
#' (ASM = A2SM(A,SMM='CMO',CINDEX=T)); SM2A(ASM)
#' (ASM2 = SM2SM(ASM,SMM2='CRI',CINDEX2=T)); A ; (A2 = SM2A(SM2))
#'
#' #clpAPI  documentation example
#' nr=5
#' nc=8
#' ra=c(3.0,5.6,1.0,2.0,1.1,1.0,-2.0,2.8,-1.0,1.0,1.0,-1.2,-1.0,1.9)
#' ia=c(0,4,0,1,1,2,0,3,0,4,2,3,0,4)
#' ja=c(0,2,4,6,8,10,11,12,14)
#' SMM='CMO'
#' CINDEX=T
#' ASM=list(nr=nr,nc=nc,ra=ra,ia=ia,ja=ja,SMM=SMM,CINDEX=CINDEX)
#' (A=SM2A(ASM))
#' }
###############################################################################
SM2A = function(ASM,NA_flag=(-1E6)){
###############################################################################
 ASM$ra[is.na(ASM$ra)] = NA_flag
##
  nr=ASM$nr; nc=ASM$nc
  if(!exists('rnames',where=ASM)) ASM$rnames=rep('',nr)
  if(!exists('cnames',where=ASM)) ASM$cnames=rep('',nc)
##
 (SMM0=ASM$SMM)

 if(SMM0=='A') A = ASM$A

 if(SMM0!='A'){

  ASM2=SM2SM(ASM,SMM2='CRI',CINDEX2=F)
  nr=ASM2$nr;nc=ASM2$nc;
  ia=ASM2$ia;ja=ASM2$ja;ra=ASM2$ra


  (ra=as.numeric(ra))
  (ia=as.integer(ia))
  (ja=as.integer(ja))

  (nr=as.integer(nr))
  (nc=as.integer(nc))

  A=matrix(0,nr,nc)
  A[cbind(ia,ja)]=ra

 } # end if(SMM0!='A')
##
  A[A <= NA_flag] = NA
##
  if(!is.null(ASM$rnames)&length(ASM$rnames)==nrow(A)) rownames(A)=ASM$rnames
  if(!is.null(ASM$cnames)&length(ASM$cnames)==ncol(A)) colnames(A)=ASM$cnames
##
  return(A)
} # end function SM2A
###############################################################################


###############################################################################
#' SM2SM: Convert a sparse matrix form into a different sparse matrix form.
#'
#' ASM is a sparse matrix object containing:
#' nr = number of rows, nc = number of columns, ra = nonzero coeff,
#'  ia = row indices, ja = col indices, SMM = space matrix method,
#'  CINDEX=use 0-indexing
#' @param ASM   A sparse matrix object
#' @param SMM2  sparse matrix method with SMM2 ='A','CRI','RCI','CMO',or'RMO'
#' @param CINDEX2   T = Use zero indexing or F = Use one indexing.
#' @param NAflag   Number to uses as NA flag
#'
#' @return A list object containing the following sparse matrix components:
#' @return nnz    = the number of non-zero elements in ra.
#' @return nr     = the number of rows in matrix.
#' @return nc     = the number of columns in matrix.
#' @return ia     = the revised row index.
#' @return ja     = the revised column index.
#' @return rnames = row names  -- may be missing or ""
#' @return cnames = column names -- may be missing or ""
#' @return ra     = the revised non-zero coefficients in A
#' @return SMM    = the revised sparse matrix type.
#' @return CINDEX = T or F for the revised sparse form.
#' @examples
#' \dontrun{
#' (A = matrix(c(1,0,0,2,0,3,0,0,0,5,0,6),3,4,byrow=T))
#' (SM1 = A2SM(A,SMM='CMO',CINDEX=T))
#' (SM2=SM2SM(SM1,SMM='CRI',CINDEX=F))
#' SM2A(SM1)
#' SM2A(SM2)
#'
#' (A = matrix(c(1,0,0,2,0,3,0,0,0,5,0,6),3,4,byrow=T)); A[2,3]=NA; A
#' (ASM = A2SM(A,SMM='CMO',CINDEX=T)); SM2A(ASM)
#' (ASM2 = SM2SM(ASM,SMM2='CRI',CINDEX2=T)); A ; (A2 = SM2A(SM2))
#'
#' #clpAPI  documentation example
#' nr=5
#' nc=8
#' ra=c(3.0,5.6,1.0,2.0,1.1,1.0,-2.0,2.8,-1.0,1.0,1.0,-1.2,-1.0,1.9)
#' ia=c(0,4,0,1,1,2,0,3,0,4,2,3,0,4)
#' ja=c(0,2,4,6,8,10,11,12,14)
#' SMM='CMO'
#' CINDEX=T
#' ASM=list(nr=nr,nc=nc,ra=ra,ia=ia,ja=ja,SMM=SMM,CINDEX=CINDEX)
#' (A=SM2A(ASM))
#' }
###############################################################################
SM2SM = function(ASM, SMM2='CRI', CINDEX2=F, NA_flag=(-1E6)){
###############################################################################
  ASM$ra[is.na(ASM$ra)] = NA_flag
##
  nnz=ASM$nnz;nr=ASM$nr; nc=ASM$nc
  if(!exists('rnames',where=ASM)) ASM$rnames=rep('',nr)
  if(!exists('cnames',where=ASM)) ASM$cnames=rep('',nc)
##
  (rnames=ASM$rnames)
  (cnames=ASM$cnames)
## original sparse form
  (SMM0=ASM$SMM)
  
## if SMM0 indicates a matrix convert ASM to SMM='CRI' sparse form
  if(SMM0=='A') {
    (tmp=A2SM(ASM$A, SMM = 'CRI', CINDEX = F))
    ASM=merge_lists(tmp,ASM)
  } # end if(SMM0=='A')

## Pull sparse matrix elements from ASM 
  (nr=ASM$nr);(nc=ASM$nc);(ia=ASM$ia);(ja=ASM$ja);(ra=ASM$ra)
  (SMM=ASM$SMM);(CINDEX=ASM$CINDEX)

## convert "zero" to "one" indexing if needed
  if(CINDEX==T){
    ia=ia+1
    ja=ja+1
  }# end if (CINDEX == T)

  ia;ja

## convert CMO ja vector or RMO ia vector to complete index vectors     
  if(SMM=='CMO') ja=rep(1:nc,diff(ja))
  if(SMM=='RMO') ia=rep(1:nr,diff(ia))
  
  tmp=as.data.frame(cbind(ia,ja,ra))
   
## convert ASM to ASM2 form 
  if(SMM2=='CRI'|SMM2=='CMO'){
    # resort to be sure in "complete and contiguous CRI form
    (tmp=tmp[with(tmp,order(ja,ia)),])
    # pull resorted data
    ra=tmp$ra; ia=tmp$ia; ja=tmp$ja
     
    if(SMM2=='CMO'){ 
     (ja=c(match(1:nc,ja),nnz+1))
    } # end if(SMM2=='CMO')     
         
  } # end if(SMM2=='CRI'|SMM2=='CMO')
  
  if(SMM2=='RCI'|SMM2=='RMO'){
    # resort to be sure in "complete and contiguous RCI form
    (tmp=tmp[with(tmp,order(ia,ja)),])
    # pull resorted data
    ra=tmp$ra; ia=tmp$ia; ja=tmp$ja
     
    if(SMM2=='RMO'){ 
     ia=c(match(1:nr,ia),nnz+1)   
    } # end if(SMM2=='RMO')     
    
  } # end if(SMM2=='RCI'|SMM2=='RMO)

  (ra=as.numeric(ra))
  (ia=as.integer(ia))
  (ja=as.integer(ja))
  (nr=as.integer(nr))
  (nc=as.integer(nc))

  if(CINDEX2 == T) {ia = ia-1; ja = ja-1}
##
##
ASM2=ASM # inherits all objects in ASM
# put revised data into ASM2
ASM2$nnz=nnz;ASM2$nr=nr;ASM2$nc=nc;ASM2$ia=ia;ASM2$ja=ja;ASM2$ra=ra
ASM2$SMM=SMM2;ASM2$CINDEX=CINDEX2


if(SMM2=='A'){
 ASM2=ASM  # inherits all objects in  original ASM
 ASM2$A = SM2A(ASM)
 ASM2$SMM = SMM2
}
#(ASM2 = list(nnz=nnz,nr = nr,nc = nc,ia = ia2,ja = ja2,ra = ra2,
#            rnames=rnames,cnames=cnames,SMM = SMM2,CINDEX = CINDEX2))
###############################################################################
return(ASM2)
###############################################################################
} # end function SM2SM
###############################################################################


# ###############################################################################
# # Examples:
# ###############################################################################
# # create matrix with missing data 
# (A = matrix(c(1,0,0,2,0,3,0,0,0,5,0,6),3,4,byrow=T)); A[2,3]=NA; A
# # no names for in matrix A
# c(rownames(A),colnames(A))
# # Convert to sparse matrix with NA's replaced by NA_flag
# (SM1=A2SM(A,NA_flag=-1000))
# # add names to A matrix
# rownames(A)=paste('r',1:nrow(A),sep='');colnames(A)=paste('c',1:ncol(A),sep='')
# A; SM2=A2SM(A)
# # Note below: SM2A puts " NA's"  back in matrix if original matrix had them
# # of if user indicates levels in matrix that designate missing data
# SM2A(SM1);SM2A(SM2)
# #####################
# ## test heritability of SM objects
# SM4$M = matrix(1:6,2,3)
# SM4
# (SM5 = SM2SM(SM4,SMM2='CMO'))
# #####################
# (ASM = A2SM(A,SMM='CMO',CINDEX=T)); SM2A(ASM)
# (ASM2 = SM2SM(ASM,SMM2='CRI',CINDEX2=T)); A ; (A2 = SM2A(SM2))
# ##############################################################################
# 
# ##############################################################################
# A=t(matrix(1:15,5,3))
# A[cbind(c(3,1,2,3,1,3),c(1,2,3,3,4,5))]=0
# A
# #
# # (A=matrix(c(1,0,2,0,0,3,0,4,0,5,6,0),3,4,byrow=T))
# # (A=matrix(c(1,0,0,2,0,3,4,0,0,5,0,6),3,4,byrow=T))
# # (A=matrix(c(1,0,0,2,0,3,0,0,0,5,0,6),3,4,byrow=T))
# # (A=matrix(c(2,3,1,4,1,2,3,4,2),3,3,byrow=T))
# #
# A2SM(A,SMM='CMO',CINDEX=T)
# 
# SMMlist=c('RCI','CRI','RMO','CMO')
# CINDEXlist=c(F,T)
# 
# SMM1 = 'RCI'     # 'RCI', 'CRI', 'RMO', 'CMO'
# CINDEX1 = F      # T or F
# SMM2='CRI'       # 'RCI', 'CRI', 'RMO', 'CMO'
# CINDEX2 = F      # T or F
# 
# for(SMM1 in SMMlist){
#  for(CINDEX1 in CINDEXlist) {
#   for(SMM2 in SMMlist) {
#    for(CINDEX2 in CINDEXlist) {
#     print(paste(SMM1,CINDEX1,SMM2,CINDEX2))
#     (ASM=A2SM(A,SMM=SMM1,CINDEX=CINDEX1))
#     (ASM1=SM2SM(ASM,SMM2=SMM2,CINDEX2=CINDEX2))
#     (ASM2=A2SM(A,SMM=SMM2,CINDEX=CINDEX2))
#     print(summary(unlist(ASM1)==unlist(ASM2)))
#    } # end for(CINDEX2 in CINDEXlist)
#   } # end for(SMM2 in SMMlist)
#  } # end for(CINDEX1 in CINDEXlist)
# } # end for(SMM1 in SMMlist)
# 
# 
# summary(SM2A(ASM1)==SM2A(ASM2))
# #############################################################
# (A=matrix(c(1,0,2,0,0,3,0,4,0,5,6,0),3,4,byrow=T))
# 
# (ASM=A2SM(A,SMM='A'))
# ASM2=SM2SM(ASM,SMM2=SMM2,CINDEX2=CINDEX2)
# (A2=SM2A(ASM2))
# summary(as.vector(A)==as.vector(A2))
# 
#  for(SMM2 in SMMlist) {
#   for(CINDEX2 in CINDEXlist) {
#    print(paste(SMM2,CINDEX2))
#    (ASM=A2SM(A=A,SMM='A'))
#    (ASM2=SM2SM(ASM,SMM2=SMM2,CINDEX2=CINDEX2))
#    (A2=SM2A(ASM2))
#    print(summary(as.vector(A)==as.vector(A2)))
#   } # end for(CINDEX2 in CINDEXlist)
#  } # end for(SMM2 in SMMlist)
# ###########################################################
# (A=matrix(c(1,0,2,0,0,3,0,4,0,5,6,0),3,4,byrow=T))
# SMM2='CRI'      #  'RCI', 'CRI', 'RMO', 'CMO'
# CINDEX2 = F     # T or F
# 
# (ASM=A2SM(A,SMM='A'))
# ASM2=SM2SM(ASM,SMM2=SMM2,CINDEX2=CINDEX2)
# (A2=SM2A(ASM2))
# summary(as.vector(A)==as.vector(A2))
# (A3=SM2A(ASM))
# summary(as.vector(A)==as.vector(A3))
# ##########################################################
# # clpAPI  documentation example
# nr=5
# nc=8
# ra=c(3.0,5.6,1.0,2.0,1.1,1.0,-2.0,2.8,-1.0,1.0,1.0,-1.2,-1.0,1.9)
# ia=c(0,4,0,1,1,2,0,3,0,4,2,3,0,4)
# ja=c(0,2,4,6,8,10,11,12,14)
# SMM='CMO'
# CINDEX=T
# ASM=list(nr=nr,nc=nc,ra=ra,ia=ia,ja=ja,SMM=SMM,CINDEX=CINDEX)
# (A=SM2A(ASM))
# #############################################################################





###############################################################################
#' cbindSM: "column bind" two sparse matrices
#'
#' @param   SM1     First sparse matrix object
#' @param   SM2     Second sparse matrix object
#' @param   SMM     Sparse matrix method with 'CRI','RCI','CMO',or'RMO'
#' @param   CINDEX  T = Use zero indexing or F = Use one indexing.
#'
#' @return A 'column bound" sparse matrix object with components"
#' @return nnz     = the number of non-zero elements in ra.
#' @return nr      = the number of rows in matrix A.
#' @return nc      = the number of columns in matrix A.
#' @return ia      = the row index.
#' @return ja      = the column index.
#' @return ra      = the non-zero coefficients in A
#' @return rnames  = matrix row names -- may be ''
#' @return cnames  = matrix column names -- may be ''
#' @return SMM     = the sparse matrix type.
#' @return CINDEX with T = Use zero indexing or F = Use one indexing.
###############################################################################
cbindSM=function(SM1,SM2,SMM='CRI',CINDEX=F){
 nr1=SM1$nr; nc1=SM1$nc; nr2=SM2$nr; nc2=SM2$nc
 if(nr1!=nr2) stop('STOP: row numbers not compatible')
 tmp1=SM2SM(SM1,SMM='CRI',CINDEX=F)
 tmp2=SM2SM(SM2,SMM='CRI',CINDEX=F)
 nr=nr1; nc = nc1+nc2
 ia=c(tmp1$ia,tmp2$ia)
 ja=c(tmp1$ja,tmp2$ja+nc1)
 ra=c(tmp1$ra,tmp2$ra)
 tmp1$nnz=length(ra);tmp1$nr=nr;tmp1$nc=nc
 tmp1$ia=ia;tmp1$ja=ja;tmp1$ra=ra;tmp1$SMM='CRI';tmp1$CINDEX=F

 if(length(tmp1$cnames)>1) {
   if(length(tmp2$cnames)>1)  tmp1$cnames=c(tmp1$cnames,tmp2$cnames)
   if(length(tmp2$cnames)<=1) tmp1$cnames=c(tmp1$cnames,rep('',nc2))
 } # end

 SM2SM(tmp1,SMM=SMM,CINDEX=CINDEX)
} # end function cbindSPM
###############################################################################

###############################################################################
#' rbindSM: "row bind" two sparse matrices
#'
#' @param   SM1     First sparse matrix object
#' @param   SM2     Second sparse matrix object
#' @param   SMM     Sparse matrix method with 'CRI','RCI','CMO',or'RMO'
#' @param   CINDEX  T = Use zero indexing or F = Use one indexing.
#'
#' @return  A 'row bound" sparse matrix object with components"
#' @return nnz     = the number of non-zero elements in ra.
#' @return nr      = the number of rows in matrix A.
#' @return nc      = the number of columns in matrix A.
#' @return ia      = the row index.
#' @return ja      = the column index.
#' @return ra      = the non-zero coefficients in A
#' @return rnames  = matrix row names -- may be ''
#' @return cnames  = matrix column names -- may be ''
#' @return SMM     = the sparse matrix type.
#' @return CINDEX with T = Use zero indexing or F = Use one indexing.
###############################################################################
rbindSM=function(SM1,SM2,SMM='CRI',CINDEX=F){
 nr1=SM1$nr; nc1=SM1$nc; nr2=SM2$nr; nc2=SM2$nc
 if(nc1!=nc2) stop('STOP: column number not compatible')
 tmp1=SM2SM(SM1,SMM='RCI',CINDEX=F)
 tmp2=SM2SM(SM2,SMM='RCI',CINDEX=F)
 nr=nr1+nr2; nc = nc1
 ia=c(tmp1$ia,tmp2$ia+nr1)
 ja=c(tmp1$ja,tmp2$ja)
 ra=c(tmp1$ra,tmp2$ra)
 tmp1$nnz=length(ra);tmp1$nr=nr;tmp1$nc=nc
 tmp1$ia=ia;tmp1$ja=ja;tmp1$ra=ra;tmp1$SMM='RCI';tmp1$CINDEX=F

 if(length(tmp1$rnames)>1) {
   if(length(tmp2$rnames)>1)  tmp1$rnames=c(tmp1$rnames,tmp2$rnames)
   if(length(tmp2$rnames)<=1) tmp1$rnames=c(tmp1$rnames,rep('',nr2))
 } # end

 SM2SM(tmp1,SMM=SMM,CINDEX=CINDEX)
} # end function rbindSPM
###############################################################################
#
#
#
################################################################################
# # Examples of cbindSM and rbindSM
#  set.seed(1001)
#  nr=7;nc=7
#  (A=matrix(round(runif(nr*nc),3),nr,nc))
#  A[A<0.5]=0
#  colnames(A)=paste('c',1:nc,sep='')
#  rownames(A)=paste('r',1:nr,sep='')
#  A
#
#  (A11=A[1:5,1:5]); (A12=A[1:5,6:7])
#  (A21=A[6:7,1:5]); (A22=A[6:7,6:7])
#
#
#  (SM11=A2SM(A11))
#  (SM12=A2SM(A12))
#  (SM21=A2SM(A21))
#  (SM22=A2SM(A22))
#
#  cbind(A11,A12)
#  SM2A(cbindSM(SM11,SM12))
#
#
#  rbind(A11,A21)
#  SM2A(rbindSM(SM11,SM21))
#
#  A
#  SM2A(cbindSM(rbindSM(SM11,SM21),rbindSM(SM12,SM22)))
###############################################################################
#  A2SM(A,SMM='CMO')
#  cbindSM(rbindSM(SM11,SM21),rbindSM(SM12,SM22),SMM='CMO')
###############################################################################











###############################################################################
#' merge_lists: Merges two list objects.
#' This function appends any objects in list2 and not in list 1
#' to list1 with priority given to list 1 components.
#'
#' @param  list1    A list object
#' @param  list2    A list object
#' @return list3    A merged object
#' @examples
#' \dontrun{
#' L1=list(a=1:3,b=4:6,c=7:9,e=10:12)
#' L2=list(b=13:15,d=16:18,e=19:21,f=22:24)
#' (L3=merge_lists(L1,L2))
#' (L4=merge_lists(L2,L1))
#' }
###############################################################################
merge_lists = function(list1, list2){
#####################################################################
# This function appends any objects in list2 and not in list 1
# to list1 with priority given to list 1 components.
#####################################################################
   (names1 = names(list1))
   (names2 = names(list2))

   (list3 = list1)
   (n1=length(list1))

   (pick2 = which(is.na(match(names2, names1))))

   if (length(pick2) > 0) {
     for(j in 1:length(pick2)) list3 = c(list3,list2[pick2[j]])
   }

   sort_list(list3)
} # end function merge_lists
###############################################################################

###############################################################################
#' sort_list: Sorts objects in list by object name
#'
#' @param  list1    A list object
#' @return list2    A resorted list object
###############################################################################
sort_list=function(list1){
   (n=length(list1))
   df1=data.frame(Lnames=names(list1))
   df1$cell=1:n;  df1$cell=as.numeric(df1$cell)
   df1=df1[order(df1$Lnames),]
   list2=list1[df1$cell[1]]
   for(j in 2:n) list2=c(list2,list1[df1$cell[j]])
   list2
} # end function sort_list
################################################################################


