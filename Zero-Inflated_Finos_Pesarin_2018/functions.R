compute_orbit_cardinality <- function(tested_variable=1,data=dati[-1,]){
  Strata=apply(data[,-tested_variable],1,paste,collapse="_")
  Strata=as.numeric(factor(Strata))
  Strata_freq=table(Strata)
  
  Strata_size=sapply(Strata_freq,factorial)
  prod(Strata_size)
}


# compute_stat_lm <- function(id,keep_coef=2,Y01,Yno0){
#   X=cbind(dataX$X[,keep_coef],dummies)
#   X=X[id,]
#   X=X[dataX$Y!=0,]
#   n=nrow(X)
#   sv=svd(X[,-(1:length(keep_coef))])
#   sv$u=sv$u[,sv$d>.00001]
#   sv$v=sv$v[,sv$d>.00001]
#   sv$d=sv$d[sv$d>.00001]
#   
#   H=diag(n)-sv$u%*%diag(sv$d^-1)%*%t(sv$u)
#   # dim(solve(t(X)%*%X)%*%t(X))
#   # length(Yno0)
#   (t(X[,(1:length(keep_coef))])%*%H%*%Yno0)/n
# }
# 
# compute_stat_01glm <- function(id,keep_coef=2,Y01,Yno0){
#   md=glm(Y01[id]~dataX$X[,keep_coef]+dummies+0)
#   smd=summary(md)
#   smd$coefficients[1,"t value"]
# }


compute_stat_quant_naive <- function(id,keep_coef=2,Yno0,X,STRATA){
  Xno0=X[id,,drop=FALSE]
  STRATA=STRATA[id]
  id_Ynot0=which(!is.na(Yno0))
  n=length(id_Ynot0)
  for(istr in unique(STRATA)){
    selezionati=which(STRATA==istr)
    selezionati_not0=intersect(id_Ynot0,selezionati)
    if(length(selezionati_not0)>0){
      Xno0[selezionati_not0,]=scale(Xno0[selezionati_not0,,drop=FALSE])
      if(any(is.na(Xno0[selezionati_not0,]))){
        Xno0[selezionati_not0,]=ifelse(is.na(Xno0[selezionati_not0,]),0,Xno0[selezionati_not0,])
      }
    }
  }
  # browser()
  stats=colMeans(Xno0*matrix(Yno0,ncol=ncol(Xno0),nrow = nrow(Yno0)),na.rm = TRUE)/
    sqrt(apply(Xno0[id_Ynot0,,drop=FALSE],2,var)*var(Yno0,na.rm = TRUE))
  stats=stats*sqrt(n)
  return(stats)
}

compute_stat_01_naive <- function(id,Y01,X){
  t(X[id,,drop=FALSE])%*%Y01/nrow(Y01)
}



compute_test_stats_pvalues <- function(tested_var=2){
  STRATA <- factor(aaply(dataX$X[,-tested_var],1,paste,collapse="_"))
  Yno0=dataX$Y
  Yno0[Yno0==0]=NA
  # Yno0=log(Yno0)
  Y01=Y01_orig
  X=dataX$X[,tested_var,drop=FALSE]
  X01=X
  for(istr in unique(STRATA)){
    selezionati=which(STRATA==istr)
    media=mean(Y01[selezionati])
    Y01[selezionati]=Y01[selezionati]-media
    media=mean(X[selezionati])
    X01[selezionati,drop=FALSE]=X01[selezionati,drop=FALSE]-media
    if(length(intersect(id_Ynot0,selezionati))>0){
      media=mean(Yno0[intersect(id_Ynot0,selezionati)])
      Yno0[intersect(id_Ynot0,selezionati)]=Yno0[intersect(id_Ynot0,selezionati)]-media
    }
  }
  # Yno0=Yno0[dataX$Y!=0,]
  # dummies=model.matrix(~.,data=data.frame(STRATA))
  # Create permutation space
  permMtrx=matrix(NA,length(STRATA),perms)
  tab=table(STRATA)
  singleton=names(which(tab==1))
  singleton_id=which(STRATA%in%singleton)
  permMtrx[singleton_id,]=singleton_id
  strata=setdiff(names(tab),singleton)
  for(strt in strata){
    idobs=which(STRATA==strt)
    for(i in 1:perms){
      permMtrx[idobs,i]=sample(idobs)
    }
  }
  permMtrx=cbind(1:nrow(permMtrx),permMtrx)
  # dim(permMtrx)
  # debug(compute_stat_quant_naive)
  # compute_stat_quant_naive(permMtrx[,1],keep_coef=tested_var,Yno0,X)
  
  permT_quant=aaply(permMtrx,2,compute_stat_quant_naive,keep_coef=tested_var,Yno0,X,STRATA)
  permT_01=aaply(permMtrx,2,compute_stat_01_naive,Y01,X01)
  permT_01=permT_01/sqrt(apply(X,2,var)*rep(var(Y01),ncol(X))/nrow(X))
  permT=cbind(permT_01,permT_quant)
  colnames(permT)=paste0(colnames(permT),rep(c("Zeros","Quant"),ncol(X)))
  permT#[1,],p=t2p(permT,obs.only = 1,tail = 0))
}



######### utility for flip
cFlip <- function(...) {
  res=list(...)[[1]]
  if(length(list(...))>1){
    nperms=sapply(list(...),function(xx) nrow(xx@permT))
    if(length(unique(nperms))>1) {
      warning("The flip-objects have different number of permutations, the minimum number will be retained for each test.")
      nperms=min(nperms)
    } else nperms=nperms[1]
    
    for(i in 2:length(list(...)))  res@permT=cbind(res@permT[1:nperms,,drop=FALSE],list(...)[[i]]@permT[1:nperms,,drop=FALSE])
    
    res@tail = as.vector(unlist(sapply(1:length(list(...)), function(i)  rep(if(is.null(list(...)[[i]]@tail)) 0 else list(...)[[i]]@tail,length.out=ncol(list(...)[[i]]@permT))
    )))
    # migliore questo output, ammettere la presenza di altri elementi in extraInfoPre
    resNames=unique(unlist(sapply(list(...),function(xx) colnames(xx@res))))
    resNames=c(setdiff(resNames,c("Stat","p-value")),c("Stat","p-value"))
    res@res[,setdiff(resNames,colnames(res@res))]=NA
    res@res=res@res[,resNames]
    for(i in 2:length(list(...)))  {
      res@res[nrow(res@res)+(1:nrow(list(...)[[i]]@res)),colnames(list(...)[[i]]@res)]=list(...)[[i]]@res
    }
  }
  res
}
