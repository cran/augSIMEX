fastapproach<-function(family){
  if (missing(family)) return(list(scorefun=score.modifieduser,fast=TRUE))
  if (family$family=="gaussian"){
    if (family$link=="identity") {return(list(scorefun=scoregaussiancpp,fast=TRUE))}
  }
  if (family$family=="binomial"){
    if (family$link=="logit") {return(list(scorefun=scorelogitcpp,fast=TRUE))}
    if (family$link=="probit") {return(list(scorefun=scoreprobitcpp,fast=TRUE))}
    if (family$link=="cloglog") {return(list(scorefun=scorecloglogcpp,fast=TRUE))}
  }
  if (family$family=="poisson"){
    if (family$link=="log") {return(list(scorefun=scorelogitcpp,fast=TRUE))}
  }
  return(return(list(scorefun=score.modifiedglm,fast=FALSE)))
}

Getalpha<-function(ValidationData,pimodel=NULL,qimodel=NULL,mis.var,mis.true,err.var){
  ValidationZ1<-ValidationData[ValidationData[,mis.true]==1,]
  ValidationZ1.data<-data.frame(pi=1-(ValidationZ1[,mis.var]==ValidationZ1[,mis.true])*1,ValidationZ1)
  ValidationZ0<-ValidationData[ValidationData$Z==0,]
  ValidationZ0.data<-data.frame(qi=1-(ValidationZ0[,mis.var]==ValidationZ0[,mis.true])*1,ValidationZ0)
  
  if (is.null(pimodel)) {pimodel<-formula(ValidationZ1.data[!colnames(ValidationZ1.data) %in% c(err.var,mis.true,mis.var)])}
  if (is.null(qimodel)) {qimodel<-formula(ValidationZ0.data[!colnames(ValidationZ0.data) %in% c(err.var,mis.true,mis.var)])}
  
  alphahat1<-glm(pimodel,family=binomial,data=ValidationZ1.data)
  alphahat2<-glm(qimodel,family=binomial,data=ValidationZ0.data)
  return(list(alphahat1=alphahat1,alphahat2=alphahat2))
}

imputeX<-function(maindata,ValidationData,err.true,err.var,lambda,Sigma_ehat,nsize,repeated,repind){
  if (!repeated){
  e_ib<-mvrnorm(n=nsize,mu=rep(0,length=dim(Sigma_ehat)[1]),Sigma=Sigma_ehat)
    return(as.matrix(maindata[,err.true])+e_ib*sqrt(lambda))} else {
      Wbi<-lapply(1:length(repind),FUN=function(i){
        if (is.null(names(repind))) repvar<-repind[[i]] else repvar<-repind[[err.var[i]]]
        Wij<-as.matrix(maindata[,repvar])
        mi_i<-rowSums(1-is.na(maindata[,repvar]),na.rm = T)
        ci<-sqrt(lambda/mi_i/(mi_i-1))
        
        nrepvar<-length(repvar)
        e_ib<-matrix(rnorm(n=nsize*nrepvar),ncol=nrepvar)
        mean_eib<-rowMeans(e_ib,na.rm=T)
        sd_eib<-apply(e_ib,MARGIN = 1,function(x) {sd(x,na.rm=T)})
        mean_eib_M<-matrix(rep(mean_eib,each=nsize),ncol=nrepvar)
        sd_eib_M<-matrix(rep(sd_eib,each=nsize),ncol=nrepvar)
        T_bij<-(e_ib-mean_eib)/sd_eib
        return(rowMeans(Wij)+ci*rowSums(T_bij*Wij,na.rm = T))
      })
      return(matrix(unlist(Wbi),ncol=length(repind),byrow=T))
    }
}


score.glm<-function(beta,Y,DataM,weight,offset,linkinv,var,mueta){
  results<-lapply(1:dim(DataM)[2],FUN=function(i){
    S<-apply(cbind(Y,DataM),MARGIN=1,function(x){
      eta<- matrix(beta,nrow=1) %*% matrix(as.numeric(x[2:length(x)]),ncol=1)
      mu<-linkinv(eta)
      return(weight[i]*mueta(eta)/(var(mu))*(x[1]-mu)* x[i+1])})
    return(S)}
  )
  return(matrix(unlist(results),ncol=dim(DataM)[2]))
}

score.modifiedglm<-function(beta,Y,DataM,DataM0,DataM1,phat,qhat,weight,linkinv,var,mueta){
  ScoreZ0<-score.glm(beta,Y,DataM0,weight,offset,linkinv,var,mueta)
  ScoreZ1<-score.glm(beta,Y,DataM1,weight,offset,linkinv,var,mueta)
  Gere<-lapply(1:(dim(DataM)[2]),FUN = function(i){
    value<-((1-DataM[,dim(DataM)[2]])*(ScoreZ0[,i]*(1-phat)-ScoreZ1[,i]*qhat)-DataM[,dim(DataM)[2]]*(ScoreZ0[,i]*phat-ScoreZ1[,i]*(1-qhat)))/(1-phat-qhat)
    return(sum(value))
  })
  return(unlist(Gere))
}

score.modifieduser<-function(beta,Y,DataM,DataM0,DataM1,phat,qhat,weight,offset,sfun){
  ScoreZ0<-sfun(beta,Y,DataM0,weight,offset)
  ScoreZ1<-sfun(beta,Y,DataM1,weight,offset)
  Gere<-lapply(1:(dim(DataM)[2]),FUN = function(i){
    value<-((1-DataM[,dim(DataM)[2]])*(ScoreZ0[,i]*(1-phat)-ScoreZ1[,i]*qhat)-DataM[,dim(DataM)[2]]*(ScoreZ0[,i]*phat-ScoreZ1[,i]*(1-qhat)))/(1-phat-qhat)
    return(sum(value))
  })
  return(unlist(Gere))
}
