augSIMEX<-function(mainformula = formula(data), pimodel = NULL, qimodel = NULL, meformula = NULL, family = gaussian,
                   data,validationdata,
                   err.var, mis.var, err.true, mis.true, err.mat = NULL, cppmethod = TRUE,
                   repeated = FALSE, repind = list(),
                   subset,  offset, weights, na.action, scorefunction=NULL,
                   lambda = NULL, M = 50, B = 200, nBoot = 50, extrapolation = c("quadratic","linear"),...)

{ call <- match.call()

  if (!missing(family)){
    ### makesure family is matched
    if (is.character(family))
      family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family))
      family <- family()
    if (is.null(family$family)) {
      print(family)
      stop("'family' not recognized")
    }
    
    ### Obtain the function in the family object
    variance <- family$variance
    linkinv <- family$linkinv
    mu.eta <- family$mu.eta
    if (!is.function(variance) || !is.function(linkinv))
      stop("'family' argument seems not to be a valid family object",
           call. = FALSE)
    dev.resids <- family$dev.resids
    aic <- family$aic
    
    unless.null <- function(x, if.null) if (is.null(x))
      if.null
    else x
    valideta <- unless.null(family$valideta, function(eta) TRUE)
    validmu <- unless.null(family$validmu, function(mu) TRUE)
    
    ### check if fast approach available
    cppmethod<-fastapproach(family)$fast
    scorefun<-fastapproach(family)$scorefun
  } else {
    if (is.null(scorefunction)) stop("Family is missing and scorefunction is not specified. ")
    sfun <-scorefunction
    scorefun<-score.modifieduser
  }

 

  ### setup the weights and offset
  if (missing(weights)) weights<-NULL
  if (missing(offset))  offset<-NULL


  ### make data into design matrix
  if (missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("mainformula", "data", "subset","weights", "na.action", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  names(mf)[2]<-"formula"
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  Response <- model.response(mf, "any")
  nsize<-length(Response)
  if (length(dim(nsize)) == 1L) {
    nm <- rownames(nsize)
    dim(nsize) <- NULL
    if (!is.null(nm))
      names(nsize) <- nm
  }
  mainCovariates <- if (!is.empty.model(mt))
    model.matrix(mt, mf, contrasts) else matrix(, nsize, 0L)
  if (all(mainCovariates[,1]==1)) intercept<-TRUE
  weights <- as.vector(model.weights(mf))
  if (!is.null(weights) && !is.numeric(weights))
    stop("'weights' must be a numeric vector")
  if (!is.null(weights) && any(weights < 0))
    stop("negative weights not allowed")
  offset <- as.vector(model.offset(mf))
  if (!is.null(offset)) {
    if (length(offset) != nsize)
      stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                    length(offset), nsize), domain = NA)
  }else {offset<- rep.int(0, nsize)}
  if (is.null(weights)) weights <- rep.int(1, nsize)
   else { if (length(weights) != nsize)
     stop(gettextf("number of weights is %d should equal %d (number of observations)",
                       length(weights), nsize), domain = NA)}


  ### rearrange the validationData
  val.covariates<-c(mis.var,err.var,err.true, mis.true)
  if (missing(validationdata)){ validationdata<-mget(unique(val.covariates),envir=parent.frame())}
  if (!repeated){
    if (!all(val.covariates %in% names(validationdata)))
      {stop(paste("Variable",
                   val.covariates[!val.covariates %in% names(validationdata)],
                   "is not included in the validation data."))}}
   else {if (!all(c(mis.var,mis.true) %in% names(validationdata))) {
     stop(paste("Variable",
                mis.var, "or", mis.true,
                "is not included in the validation data."))
   }}
  colname<-colnames(data)

  ## SIMEXvariable
  err.var = unique(err.var)
  nerr.var = length(err.var)
  mis.var = unique(mis.var)
  nmis.var = length(mis.var)
  if (!is.character(err.var) | nerr.var > length(colname)) {
    stop("Invalid err.var object")
  }
  if (!is.character(mis.var) | nmis.var > length(colname)) {
    stop("Invalid mis.var object")
  }
  if (!all(err.var %in% colname)) {
    stop("Error-prone variable must be selected from the data")
  }
  if (!all(mis.var %in% colname)) {
    stop("Misspecified variable must be included in the data")
  }
  # if (!(repeated == FALSE) & !(repeated == TRUE)) {
  #   stop("Repeated indicator should only be 'TRUE' or 'FALSE'. ")
  # }
  if (!is.null(err.mat)) {
    err.mat = as.matrix(err.mat)
    if (!is.numeric(err.mat) | any(err.mat < 0)) {
      stop("Invalide err.mat object, err.mat must be a square symmetric numeric matrix")
    }
    if (nrow(err.mat) != ncol(err.mat)) {
      stop("err.mat must be a square matrix")
    }
    if (length(err.var) != nrow(err.mat)) {
      stop("SIMEXvariable and err.mat have non-conforming size")
    }
    SSigma <- err.mat
    dimnames(SSigma) <- NULL
    if (!isTRUE(all.equal(SSigma, t(SSigma)))) {
      warning("err.mat is numerically not symmetric")
    }
  }


  ###SIMEX parameters
  if (length(B) != 1) {
    stop("B must be positive integer")
  }
  if (!is.numeric(B) | B <= 0) {
    stop("B must be positive integer")
  }
  else {
    B = ceiling(B)
  }
  if (is.null(lambda)) {lambda<-seq(from=0,to=1,length.out=M)}
  if (!is.vector(lambda) | !is.numeric(lambda)) {
    stop(":Invalide lambda object")
  }
  if (any(lambda < 0)) {
    warning("Lambda should be positive values. Negative values will be ignored",
            call. = FALSE)
    lambda <- lambda[lambda >= 0]
  }

  extrapolation<-match.arg(extrapolation)

  ### Step 1: SImulation step
  temp.Results1<-Getalpha(validationdata,pimodel,qimodel,mis.var,mis.true,err.var)
  alphahat1<-temp.Results1$alphahat1
  alphahat2<-temp.Results1$alphahat2

  ### estimate the empirical covariance matrix
  if (!repeated){  
    if (missing(err.mat)|is.null(err.mat)){
      if (length(meformula)==0) {
      Models_res<-lapply(1:length(err.var),FUN = function(i){
        model<-lm(validationdata[,err.var[i]]~-1+.,data=data.frame(validationdata[,err.true]))
        return(model$residuals)
      })} else {
        if (length(meformula)!=length(err.var)) {stop(gettextf("the length of meformula is %d and should be equal to the number of error-prone covariates %d",
                                                  length(meformula), length(err.var)), domain = NA)}
      resname.me<-unlist(lapply(1:length(err.var),FUN = function(i){return(all.vars(meformula[[i]])[1])}))
      index<-match(err.var,resname.me)
      if (length(index)!=length(err.var)) {stop("incorrect measurement error model specified")}
      Models_res<-lapply(index,FUN = function(i){
          model<-lm(meformula[[i]],data=validationdata)
          return(model$residuals)
        })
      }
      Models_res<-matrix(unlist(Models_res),ncol=length(err.var))
      Sigma_ehat<-cov(Models_res)
    } else Sigma_ehat<-cov(err.mat)
  }

  ### change the Xstar in mainCovariates into X, Zstar into Z
  changeindex<-match(err.var,colnames(mainCovariates))
  colnames(mainCovariates)[changeindex]<-err.true
  colnames(mainCovariates)[which(colnames(mainCovariates)==mis.var)]<-mis.true

  main.new<-mainCovariates
  Wnames<- setdiff(colnames(main.new),c(val.covariates,unlist(repind)))
  Wmatrix<-main.new[,Wnames]
  NSam<-dim(main.new)[1]
  nbeta<-dim(mainCovariates)[2]
  if (repeated) imputeData<-data else imputeData<-mainCovariates

  betahat<-lapply(1:M,FUN=function(i){
    betab_lam<-lapply(1:B,FUN=function(x){
      X.impute<-imputeX(imputeData,validationdata,err.true,err.var,lambda[i],Sigma_ehat,nsize,repeated,repind)
      main.new[,err.true]<-X.impute
      main.new.df<-data.frame(main.new)

      phat0<-predict(alphahat1,newdata=main.new.df,type="response")
      qhat0<-predict(alphahat2,newdata=main.new.df,type="response")

      DataM <- cbind(X.impute,Wmatrix,main.new[,mis.true])
      DataM0 <- cbind(X.impute,Wmatrix,rep(0,nsize))
      DataM1 <- cbind(X.impute,Wmatrix,rep(1,nsize))

      if (is.null(scorefunction)){
        if (cppmethod) {betasolve<-multiroot(scorefun,rep(0,nbeta),Y=Response,DataM=DataM,DataM0=DataM0,DataM1=DataM1,phat=phat0,qhat=qhat0,weight=weights,offset=offset)}
         else {betasolve<-multiroot(scorefun,rep(0,nbeta),Y=Response,DataM=DataM,DataM0=DataM0,DataM1=DataM1,phat=phat0,qhat=qhat0,weight=weights,offset=offset,linkinv=linkinv,var=variance,mueta=mu.eta)}
      } else {betasolve<-multiroot(scorefun,rep(0,nbeta),Y=Response,DataM=DataM,DataM0=DataM0,DataM1=DataM1,phat=phat0,qhat=qhat0,weight=weights,offset=offset,sfun=sfun)}

      if (any(abs(betasolve$root)>8)) { return(rep(NA,nbeta*2))}

      betab_lam_boot<-matrix(nrow=nBoot,ncol=nbeta)

      for (t in 1:nBoot){
        sample.boot<-sample(NSam,replace=T)
        repn<-length(Response[sample.boot])
        Wboot<-as.matrix(main.new[sample.boot,Wnames])
        DataM <- cbind(X.impute[sample.boot,],Wboot,main.new[sample.boot,mis.true])
        DataM0 <- cbind(X.impute[sample.boot,],Wboot,rep(0,repn))
        DataM1 <- cbind(X.impute[sample.boot,],Wboot,rep(1,repn))
        
        if (is.null(scorefunction)){
          if (cppmethod) {betasolve_boot<-multiroot(scorefun, betasolve$root,Y=Response[sample.boot],DataM=DataM,DataM0=DataM0,DataM1=DataM1,phat=phat0[sample.boot],qhat=qhat0[sample.boot],weight=weights,offset=offset)}
            else {betasolve<-multiroot(scorefun,betasolve$root,Y=Response[sample.boot],DataM=DataM,DataM0=DataM0,DataM1=DataM1,phat=phat0[sample.boot],qhat=qhat0[sample.boot],weight=weights,offset=offset,linkinv=linkinv,var=variance,mueta=mu.eta)}
          } else {betasolve_boot<-multiroot(scorefun, betasolve$root,Y=Response[sample.boot],DataM=DataM,DataM0=DataM0,DataM1=DataM1,phat=phat0[sample.boot],qhat=qhat0[sample.boot],weight=weights,offset=offset,sfun=sfun)}
        if (any(abs(betasolve_boot$root)>8)) betab_lam_boot[t,]<-rep(NA,nbeta) else betab_lam_boot[t,]<-betasolve_boot$root
      }

      betab_lam_sd<-apply(betab_lam_boot,MARGIN = 2, FUN=var,na.rm=T)
      return(c(betasolve$root,betab_lam_sd))
    })


    betahat_lam_M<-matrix(unlist(betab_lam),ncol=nbeta*2,byrow=T)
    betahat_lam<-colMeans(betahat_lam_M[,1:nbeta],na.rm=T)
    Omegahat_lam<-colMeans(betahat_lam_M[,nbeta+1:nbeta],na.rm=T)
    Shat_lam<-apply(betahat_lam_M[,1:nbeta],MARGIN = 2, FUN=var,na.rm=T)
    return(c(betahat_lam,Omegahat_lam-Shat_lam))
  })

  betamatrix<-matrix(unlist(betahat),ncol=nbeta*2,byrow=T)
  if (extrapolation=="quadratic"){
    extrapomodel<-apply(betamatrix,MARGIN=2, FUN=function(x){
      lambda2<-lambda^2
      model<-lm(x~lambda+lambda2)
      newdata<-data.frame(lambda=-1,lambda2=1)
      betaresults<-predict(model,newdata)
    })
  } else if (extrapolation=="linear"){
    extrapomodel<-apply(betamatrix,MARGIN=2, FUN=function(x){
      model<-lm(x~lambda)
      newdata<-data.frame(lambda=-1)
      betaresults<-predict(model,newdata)
    })
  }

  good <- weights > 0


  coefs <- extrapomodel[1:nbeta]
  sds<-sqrt(extrapomodel[nbeta+1:nbeta])
  names(coefs)<-c(err.true,Wnames,mis.true)
  names(sds)<-c(err.true,Wnames,mis.true)
  colnames(betamatrix)<-c(names(coefs),names(sds))
  if (is.null(scorefunction)){
    eta <- drop(mainCovariates[,names(coefs)] %*% t(t(coefs)))
    mu <- linkinv(eta <- eta + offset)
    
    varmu <- variance(mu)[good]
    if (anyNA(varmu))
      stop("NAs in V(mu)")
    if (any(varmu == 0))
      stop("0s in V(mu)")
    mu.eta.val <- mu.eta(eta)
    if (any(is.na(mu.eta.val[good])))
      stop("NAs in d(mu)/d(eta)")
    
    residuals <- (Response - mu)/mu.eta(eta)
    dev<-sum(dev.resids(Response, mu, weights))
    wtdmu <- if (intercept) sum(weights * Response)/sum(weights) else linkinv(offset)
    nulldev <- sum(dev.resids(Response, wtdmu, weights))
    ## calculate df
    n.ok <- nsize - sum(weights==0)
    nulldf <- n.ok - as.integer(intercept)
    rank <- length(coefs)
    resdf  <- n.ok - rank
    ## calculate AIC
    aic.model <- aic(Response, 0, mu, weights, dev) + 2*rank
  } else {
    eta <- NA
    mu <- NA
    varmu <- NA
    residuals <- NA
    dev<- NA
    wtdmu <- NA
    nulldev <- NA
    ## calculate df
    n.ok <- nsize - sum(weights==0)
    nulldf <- n.ok - as.integer(intercept)
    rank <- length(coefs)
    resdf  <- n.ok - rank
    ## calculate AIC
    aic.model <- NA
  }
  



  if (intercept) {
    intindex<-which(names(coefs)=="(Intercept)")
    otherindex<-1:nbeta
    otherindex<-otherindex[-intindex]
    coefs <-coefs[c(intindex,otherindex)]
    sds <-sds[c(intindex,otherindex)]
    betamatrix<-betamatrix[,c(intindex,otherindex,intindex+nbeta,otherindex+nbeta)]
  }


  output<-list(coefficients=coefs,
            se=sds,
            call = call, family=family,
            formula=mainformula, pimodel=pimodel, qimodel=qimodel,
            err.var=err.var, mis.var=mis.var, err.true=err.true, mis.true=mis.true, err.mat=err.mat,
            M=M, B=B, nBoot=nBoot, extrapolation=extrapolation,
            lambda=lambda,coefmatrix=betamatrix,linear.predictors = eta, deviance = dev, aic = aic.model,
            null.deviance = nulldev, df.residual = resdf, df.null = nulldf,residuals=residuals,
            rank=rank)

  class(output)<- "augSIMEX"

  return(output)

}
