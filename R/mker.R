#' Fit the multi-kink expectile regression with longitudinal data
#'
#' @rdname mker.bea
#' @param y A vector of response
#' @param thre.x A scalar covariate with threshold
#' @param cont.z A vector of covariates with constant slopes
#' @param tau the expectile level
#' @param J/ID ... 
#' @param Cn A psotive number corresponding to different types of BIC
#' @param control A list returned by \code{fit.control}.}
#'
#' @return A list with the elements
#' \item{bet.est}{The estimated regression coefficients with intercept.}
#' \item{bet.se}{The estimated standard error of the regression coefficients.}
#' \item{psi.est}{The estimated change points.}
#' \item{psi.se}{The estiamted standard errors of threshold parameters.}
#' \item{n.psi}{The estiamted number of change points.}
#' 
# ## simulated data
# library(MASS)
# N=200;tau=0.5;J=4
# psi0=c(-2,1,3);k=length(psi0)
# PSI=matrix(rep(psi0,rep(J*N,k)),ncol=k)
# bet0=c(1,1,1,-3,4,4)
# 
# x=runif(J*N,-5,5)
# z=rnorm(J*N,1,0.5^2)
# X=matrix(rep(x,k),nrow=J*N)
# XZ=cbind(1,z,x,pmax((X-PSI),0))
# 
# rho=0.5;mu=matrix(rep(0,J));
#sigma=matrix(NA,J,J)
#for(i in 1:J)
#{for(j in 1:J)
#{sigma[i,j]=rho^abs(i-j)}}
# #homoscedastic normal errors
# L=0;err=mvrnorm(N,mu,sigma);err=as.vector(err)-efun(tau,errtype = 1)
# 
# y=XZ%*%bet0+(1+L*z)*err
# 
# fit=mker.bea(y,x,z,tau,J,Cn=(log(N)))






fit.control <-
  function(toll=1e-4,h=1,it.max=50,K.max=6,stop.if.error=TRUE,dev0=NULL,visual=FALSE,
           visualBoot=FALSE,pow=c(1,1),digits=NULL,grid=NULL,n.boot=20){
    list(toll=toll,h=h,it.max=it.max,K.max=K.max,stop.if.error=stop.if.error,
         dev0=dev0,visual=visual,n.boot=n.boot,
         visualBoot=visualBoot,pow=pow,digits=digits,grid=grid)
  }




expTLlaws <- function(X, y, tau,  max.iter, tol)
{
  b <- bLast <- rep(0, ncol(X)) # initialize
  it <- 1 # iteration index
  while (it <= max.iter){
    
    ypred <- c(X %*% b)
    w <- as.vector(tau *(y>= ypred) + (1-tau)* (y<ypred))
    b <- lsfit(X, y, w, intercept=FALSE)$coef
    if (max(abs(b - bLast)/(abs(bLast) + 0.01*tol)) < tol) break
    bLast <- b
    it <- it + 1 # increment index
  }
  if (it > max.iter) warning('maximum iterations exceeded')
  
  res=y-c(X%*%b)
  loss <- sum(w*(y-c(X%*%b))^2)
  
  list(coefficients = b, loss = loss, residuals=res, it = it)
}



pefun <- function(t, errtype){
  
  if (errtype==1){
    ## example1: normal distribution
    F <- pnorm(t, 0, 1)
    
    integrand <- function(x) {x*dnorm(x, 0, 1)}
    
  } else if (errtype==2){
    ## example2: t4 distribution
    F <- pt(t, 4)
    
    integrand <- function(x) {x*dt(x, 4)}
    
  } else if (errtype==3) {
    ## example3: mixtrue distribution of normal distributions
    prop <- 0.1  # mixture proportion
    F <- prop * pt(t, 4) + (1-prop)*pnorm(t, 0, 1)
    
    integrand <- function(x){
      x*(prop * dt(x, 4) +  (1-prop)*dnorm(x, 0, 1))
    }
  }
  
  G <- integrate(integrand, lower = -(1e+3), upper = t)$value
  gmean <-  integrate(integrand, lower = -(1e+3), upper = 1e+3)$value
  u <- G -t * F
  asy <- u/(2*u + t-gmean)
  
  return(asy)
}

efun <- function (tau, errtype){
  tau[tau > 1 | tau < 0] = NA
  zz = 0 * tau
  lower = rep(-10, length(tau))
  upper = rep(10, length(tau))
  diff = 1
  index = 1
  while (diff > 1e-10 && index < 1000) {
    root = pefun(zz, errtype) - tau
    root[is.na(root)] = 0
    lower[root < 0] = zz[root < 0]
    upper[root > 0] = zz[root > 0]
    zz = (upper + lower)/2
    diff = max(abs(root), na.rm = T)
    index = index + 1
  }
  zz[is.na(tau)] = NA
  return(zz)
}



wfun <- function(u, tau){
  abs(tau- (u <= 0))
}




brisq <-function(y, XREG,X, PSI, tau, J, opz, n.boot=20, 
                 size.boot=NULL, jt=FALSE,
                 nonParam=TRUE, random=FALSE)
{
  extract.psi<-function(lista){
    dev.values<-lista[[1]]
    psi.values<-lista[[2]]
    dev.values <- as.vector(dev.values, mode = "numeric")
    dev.ok <- min(dev.values,na.rm = TRUE)
    id.dev.ok<-which.min(dev.values)  
    if(is.list(psi.values))  psi.values<-matrix(unlist(psi.values),
                                                nrow=length(dev.values), byrow=TRUE)
    if(!is.matrix(psi.values)) psi.values<-matrix(psi.values)
    psi.ok<-psi.values[id.dev.ok,]
    r<-list(SumSquares.no.gap=dev.ok, psi=psi.ok)
    r
  }
  #-------------
  visualBoot<- opz$visualBoot
  opz.boot<-opz
  opz.boot$pow=c(1.1,1.2)  #
  opz1<-opz
  opz1$it.max <-1
  n<-length(y)
  #x <- X[,1]
  
  
  #ris <- seg.qr.fit(y,XREG,X,Z,PSI,tau=0.5,opz=seg.control())
  o0 <- try(brisq.fit(y, XREG, X, PSI, tau, opz, return.all.sol=FALSE), silent=TRUE)
  
  rangeZ <- apply(X, 2, range)
  if(!is.list(o0)) {
    o0 <- brisq.fit(y, XREG, X, PSI, tau, opz, return.all.sol=TRUE)
    o0<-extract.psi(o0)  
    if(!nonParam) {warning("using nonparametric boot");nonParam<-TRUE}
  }
  if(is.list(o0)){
    est.psi00<-est.psi0<-o0$psi
    ss00<-o0$SumSquares.no.gap
    if(!nonParam) fitted.ok<-fitted(o0$obj)  #
  } else {  
    if(!nonParam) stop("the first fit failed and I cannot extract fitted values for the boot sample")
    if(random) {
      est.psi00<-est.psi0<-apply(rangeZ,2,function(r)runif(1,r[1],r[2]))
      PSI1 <- matrix(rep(est.psi0, rep(nrow(X), length(est.psi0))), ncol = length(est.psi0))
      o0<-try(brisq.fit(y, XREG,X, PSI, tau, opz1), silent=TRUE)
      ss00<-o0$SumSquares.no.gap
    } else {
      est.psi00<-est.psi0<-apply(PSI,2,mean)
      ss00<-opz$dev0
    }
  }
  
  all.est.psi.boot<-all.selected.psi<-all.est.psi<-matrix(NA, nrow=n.boot, ncol=length(est.psi0))
  all.ss<-all.selected.ss<-rep(NA, n.boot)               
  if(is.null(size.boot)) size.boot<-n/J   
  
  #      na<- ,,apply(...,2,function(x)mean(is.na(x)))
  
  X.orig<- X
  if(visualBoot) cat(0, " ", formatC(opz$dev0, 3, format = "f"),"", "(No breakpoint(s))", "\n")
  count.random<-0
  for(k in seq(n.boot)){
    PSI <- matrix(rep(est.psi0, rep(nrow(X), length(est.psi0))), ncol = length(est.psi0))
    if(jt)  X<-apply(X.orig,2,jitter)
    if(nonParam){
      idx<-sample(n/J, size=size.boot, replace=TRUE)   
      id=idx; for(j in 1:(J-1)){id=c(id,idx+j*(n/J))}  
      o.boot<-try(brisq.fit(y[id], XREG[id,,drop=FALSE],X[id,,drop=FALSE], PSI[id,,drop=FALSE],
                            tau=tau, opz.boot), silent=TRUE)
    } else {
      yy<-fitted.ok+sample(residuals(o0$obj),size=n, replace=TRUE)#residuals(o0)
      o.boot<-try(brisq.fit(yy, XREG, X.orig, PSI, tau=tau, opz.boot), silent=TRUE)
    }
    if(is.list(o.boot)){  
      all.est.psi.boot[k,]<-est.psi.boot<-o.boot$psi
    } else {
      est.psi.boot<-apply(rangeZ,2,function(r)runif(1,r[1],r[2]))
    }
    PSI <- matrix(rep(est.psi.boot, rep(nrow(X), length(est.psi.boot))), ncol = length(est.psi.boot))
    opz$h<-max(opz$h*.9, .2)
    opz$it.max<-opz$it.max+1
    o<-try(brisq.fit(y, XREG,X.orig,  PSI, tau=tau, opz, return.all.sol=TRUE), silent=TRUE)
    if(!is.list(o) && random){  
      est.psi0<-apply(rangeZ,2,function(r)runif(1,r[1],r[2]))
      PSI1 <- matrix(rep(est.psi0, rep(nrow(X), length(est.psi0))), ncol = length(est.psi0))
      o<-try(brisq.fit(y, XREG,X,  PSI1, tau, opz1), silent=TRUE)
      count.random<-count.random+1
    }
    if(is.list(o)){
      if(!"coefficients"%in%names(o$obj)) o<-extract.psi(o) 
      all.est.psi[k,]<-o$psi    
      all.ss[k]<-o$SumSquares.no.gap
      if(o$SumSquares.no.gap<=ifelse(is.list(o0), o0$SumSquares.no.gap, 10^12)) o0<-o
      est.psi0<-o0$psi
      all.selected.psi[k,] <- est.psi0  
      all.selected.ss[k]<-o0$SumSquares.no.gap #min(c(o$SumSquares.no.gap, o0$SumSquares.no.gap))
    }
    if(visualBoot) {
      flush.console()
      spp <- if (k < 10) "" else NULL
      cat(k, spp, "", formatC(o0$SumSquares.no.gap, 3, format = "f"), "\n")
    }
  } #end n.boot
  #browser()
  
  
  all.selected.psi<-rbind(est.psi00,all.selected.psi)
  all.selected.ss<-c(ss00, all.selected.ss)
  
  SS.ok<-min(all.selected.ss)
  id.accept<- ((abs(all.ss-SS.ok)/SS.ok )<= 0.05)  
  psi.mean<-apply(all.est.psi[id.accept,,drop=FALSE], 2, mean)  
  
  #      est.psi0<-psi.mean
  #      #devi ristimare il modello con psi.mean
  #      PSI1 <- matrix(rep(est.psi0, rep(nrow(Z), length(est.psi0))), ncol = length(est.psi0))
  #      o0<-try(seg.lm.fit(y, XREG, Z, PSI1, w, offs, opz1), silent=TRUE)
  
  ris<-list(all.selected.psi=drop(all.selected.psi),all.selected.ss=all.selected.ss,
            all.psi=all.est.psi, all.ss=all.ss)
  
  if(is.null(o0$obj)){  
    PSI1 <- matrix(rep(est.psi0, rep(nrow(X), length(est.psi0))), ncol = length(est.psi0))
    o0<-try(brisq.fit(y, XREG, X, PSI1, tau, opz1), silent=TRUE)
  }
  if(!is.list(o0)) return(0)  
  
  o0$boot.restart<-ris
  return(o0)
}



brisq.fit <-function(y,XREG,X,PSI,tau,opz,return.all.sol=FALSE)
{
  
  #-----------
  psi <- PSI[1,]
  n <- length(y)
  x <- X[,1]
  c1 <- apply((X <= PSI), 2, all)
  c2 <- apply((X >= PSI), 2, all)
  if(sum(c1 + c2) != 0 || is.na(sum(c1 + c2))) stop("psi out of the range")
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  ##~~~~~~~~~~~~~
  xreg.names <- opz$xreg.names  
  #mtype <- opz$mtype
  grid <- opz$grid   
  digits <- opz$digits  
  pow<-opz$pow          
  #nomiOK<-opz$nomiOK   ##
  toll<-opz$toll
  h<-opz$h ##
  #gap<-opz$gap  ##FALSE
  stop.if.error<-opz$stop.if.error
  dev.new<-opz$dev0
  visual<-opz$visual
  #id.psi.group<-opz$id.psi.group  #
  it.max<-old.it.max<-opz$it.max
  rangeZ <- apply(X, 2, range)
  #psi<-PSI[1,]
  #names(psi)<-id.psi.group
  #H<-1
  it <- 1
  epsilon <- 10
  dev.values<-psi.values <- NULL
  id.psi.ok<-rep(TRUE, length(psi))
  #abs(epsilon) > toll
  while (it < it.max) {
    U <- pmax((X - PSI), 0)
    V <- ifelse((X > PSI), -1, 0)
    
    ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    XZ <- cbind(XREG,U,V)
    rownames(XZ) <- NULL
    colnames(XZ) <- c(xreg.names,paste("U",1:ncol(U),sep=""),paste("V",1:ncol(V),sep=""))
    ####~~~~~~~~~~~~~~~~~~~``
    
    
    #obj <- rq(y~XZ, tau = tau,  method ="br")
    XX=cbind(1,XZ)                                          
    obj=expTLlaws(X=XX, y, tau,  max.iter=100, tol=1E-4)    
    
    dev.old<-dev.new  
    dev.new <- dev.new1 <- obj$loss
    dev.values[[length(dev.values) + 1]] <- dev.new1
    if (visual) {
      flush.console()
      if (it == 1)
        cat(0, " ", formatC(dev.old, 3, format = "f"),
            "", "(No breakpoint(s))", "\n")
      spp <- if (it < 10) "" else NULL
      cat(it, spp, "", formatC(dev.new, 3, format = "f"), "",length(psi),"\n")
      #cat(paste("iter = ", it, spp," dev = ",formatC(dev.new,digits=3,format="f"), " n.psi = ",formatC(length(psi),digits=0,format="f"), sep=""), "\n")
    }
    epsilon <- (dev.new - dev.old)/(dev.old + .001)
    obj$epsilon <- epsilon  
    it <- it + 1
    obj$it <- it   
    
    
    beta.c=obj$coefficients[(2+ncol(XREG)):(1+ncol(XREG)+length(psi))]      
    gamma.c=obj$coefficients[(2+ncol(XREG)+length(psi)):(1+ncol(XREG)+2*length(psi))]
    
    
    #if (it > it.max) break
    if(it>10 && epsilon<toll) break
    #if(max(abs(gamma.c))<1e-4) break
    
    psi.values[[length(psi.values) + 1]] <- psi.old <- psi
    #       if(it>=old.it.max && h<1) H<-h
    psi <- psi.old + h*gamma.c/beta.c
    if(!is.null(digits)) psi<-round(psi, digits)  
    PSI <- matrix(rep(psi, rep(n, length(psi))), ncol = length(psi))
    #check if psi is admissible..
    a <- apply((X <= PSI), 2, all) #prima era solo <
    b <- apply((X >= PSI), 2, all) #prima era solo >
    if(stop.if.error) {
      isErr<- (sum(a + b) != 0 || is.na(sum(a + b))) || (any(diff(sort(psi))<grid))
      if(isErr) {
        if(return.all.sol) return(list(dev.values, psi.values)) else stop("(Some) estimated psi gets wrong")
      }
    } else {
      id.psi.ok<-!is.na((a+b)<=0)&(a+b)<=0
      #X <- X[,id.psi.ok,drop=FALSE]
      psi <- psi[id.psi.ok]
      PSI <- PSI[,id.psi.ok,drop=FALSE]
      X <- X[,id.psi.ok,drop=FALSE]
      #nomiOK<-nomiOK[id.psi.ok]
      #id.psi.group<-id.psi.group[id.psi.ok]
      #names(psi)<-id.psi.group
      
      ##~~~~~~~~~~~~~~~~~~~~~~~~~~
      if(any(diff(sort(psi))<grid)){
        psi <- sort(psi)
        PSI <- PSI[,order(psi)]
        X <- X[,order(psi)]
        id.psi.grid <- which(diff(psi)<grid)+1
        psi <- psi[-id.psi.grid]
        PSI <- matrix(PSI[,-id.psi.grid],nrow=n)
        X <- matrix(X[,-id.psi.grid],nrow=n)
        #nomiOK<-nomiOK[-id.psi.grid]
        #id.psi.group<-id.psi.group[-id.psi.grid]
        #names(psi)<-id.psi.group
      }
      ##~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      #if(length(psi)<=0) return(0)
      if(length(psi)==0){
        #obj <- rq(y~XREG,tau=tau,method="br")
        XX=cbind(1,XREG)                                        
        obj=expTLlaws(X=XX, y, tau,  max.iter=100, tol=1E-4)    
        obj$n.psi <- 0
        obj$psi <- NULL
        return(obj)
      }
    } #end else
    #obj$psi <- psi
  } #end while
  psi <- sort(psi)
  names(psi) <- paste("psi",1:length(psi),sep="")
  #psi<-unlist(tapply(psi, id.psi.group, sort))
  #names(psi)<-id.psi.group
  X <- matrix(rep(x,length(psi)),nrow=n)
  PSI <- matrix(rep(psi, rep(n, length(psi))), ncol = length(psi))
  U <- pmax((X - PSI), 0)
  #V <- ifelse((X > PSI), -1, 0)
  XZ <- cbind(XREG,U)
  rownames(XZ) <- NULL
  colnames(XZ) <- c(xreg.names,paste("U",1:ncol(U),sep=""))
  
  #obj <- rq(y~XZ, tau = tau,  method ="br")
  XX=cbind(1,XZ)                                          
  obj=expTLlaws(X=XX, y, tau,  max.iter=100, tol=1E-4)    
  SS.new <- obj$loss                                      
  
  #fino a qua..
  obj<-list(obj=obj,psi=psi,psi.values=psi.values,n.psi=length(psi), SumSquares.no.gap=SS.new)
  return(obj)
}




mker.bea <- function(y,thre.x,cont.z,tau,J,Cn,control=fit.control()){
  
  n <- length(y)
  dev0 <- control$dev0
  if(missing(cont.z) || is.null(cont.z)){
    p <- 2; XREG <- matrix(thre.x,nrow=n,ncol=1);xreg.names <- c("x")
  }else{
    XREG <- cbind(cont.z,thre.x); p <- ncol(XREG)+1
    if(p==3) xreg.names <- c("z","x") else xreg.names <- c(paste("z",1:ncol(cont.z),sep=""),"x")
  }
  if(is.null(dev0)){
    #obj0 <- rq(y~XREG,tau,method="br")
    XX=cbind(1,XREG)                                        
    obj0=expTLlaws(X=XX, y, tau,  max.iter=100, tol=1E-4)   
    control$dev0 <- obj0$loss                               
  }
  if(is.null(control$grid)) control$grid <- (max(thre.x)-min(thre.x))/30
  control$xreg.names <- xreg.names
  
  k <- control$K.max
  psi0 <- quantile(thre.x,seq(0,1,l=k+2)[-c(1,(k+2))],names=FALSE)
  k <- length(psi0)
  control$stop.if.error <- FALSE
  X <- matrix(rep(thre.x,k),nrow=n)
  PSI <- matrix(rep(psi0,rep(n,k)),ncol=k)
  obj <- suppressWarnings(brisq.fit(y,XREG,X,PSI,tau,control,return.all.sol=FALSE))
  if(obj$n.psi==0) return(obj)
  psi0 <- obj$psi
  k <- obj$n.psi
  bic0 <- log(obj$obj$loss) + (p + 2*k)*log(n)/2/n*Cn  ####n  n/J
  while(k>0){
    k <- k-1
    if(k==0){
      bic1 <- log(control$dev0) + p*log(n)/2/n*Cn      #####n  n/J
      obj1 <- NULL
      obj1$n.psi <- 0;obj1$psi <- NULL
    }else{
      psi0 <- quantile(thre.x,seq(0,1,l=k+2)[-c(1,(k+2))],names=FALSE)
      control$stop.if.error <- TRUE
      X <- matrix(rep(thre.x,k),nrow=n)
      PSI <- matrix(rep(psi0,rep(n,k)),ncol=k)
      obj1 <- suppressWarnings(brisq(y,XREG,X,PSI,tau,J,control)) 
      bic1 <- log(obj1$obj$loss) + (p+2*k)*log(n)/2/n*Cn          ####n  n/J
    }
    if(bic1 > bic0) break else obj <- obj1; bic0 <- bic1
  }
  psi <- obj$psi
  
  if(length(psi)==0){                                       
    XX=cbind(1,XREG)                                        
    oo=expTLlaws(X=XX, y, tau,  max.iter=100, tol=1E-4)     
    oo$bet.est=oo$coefficients;oo$psi.est=NULL              
    oo$bet.se=NULL;oo$psi.se=NULL
  }else{                                                    
    oo <- mkqr.fit(y,thre.x,cont.z,psi,k=length(psi),tau,J)  
  }
  
  oo$n.psi <- length(psi)
  return(oo)
}


mkqr.fit <- function(y,thre.x,cont.z,psi,k,tau,J,control=fit.control()){
  
  if(missing(k) || is.null(psi)){  
    if(missing(k) || is.null(k)){
      stop("either psi or k must be given")
    }else{
      psi <- quantile(thre.x,seq(0,1,l=(k+2)))[-c(1,k+2)]
    }
  }
  k <- length(psi)
  
  n <- length(y) 
  true_n=length(y)/J     
  
  dev0 <- control$dev0  
  if(missing(cont.z) || is.null(cont.z)){
    p <- 2; XREG <- matrix(thre.x,nrow=n,ncol=1);xreg.names <- c("x")
  }else{
    XREG <- cbind(cont.z,thre.x); p <- ncol(XREG)+1
    if(p==3) xreg.names <- c("z","x") else xreg.names <- c(paste("z",1:ncol(cont.z),sep=""),"x")
  }
  if(is.null(dev0)){
    #obj0 <- rq(y~XREG,tau,method="br")
    XX=cbind(1,XREG)                                        
    obj0=expTLlaws(X=XX, y, tau,  max.iter=100, tol=1E-4)   
    control$dev0 <- obj0$loss                               
  }
  control$xreg.names <- xreg.names
  control$stop.if.error <- TRUE
  X <- matrix(rep(thre.x,k),nrow=n)
  PSI <- matrix(rep(psi,rep(n,k)),ncol=k)
  obj <- suppressWarnings(brisq(y,XREG,X,PSI,tau,J,control))  
  res <- obj$obj$residuals    
  res=as.vector(res)          
  
  psi.est <- obj$psi
  bet.est <- obj$obj$coefficients
  U.est <- bet.est[(p+1):(p+k)]
  PSI <- matrix(rep(psi.est,rep(n,k)),ncol=k)
  BB <- matrix(rep(U.est,n),nrow=n,ncol=k,byrow=T)
  
  ###############
  wt <- wfun(res, tau)
  hh=true_n^(-2)
  
  hh_theta=cbind(1,XREG,(X-PSI)*(pnorm((X-PSI)/hh)),
                 -BB*pnorm((X-PSI)/hh)-BB*((X-PSI)/hh)*dnorm((X-PSI)/hh))
  h_theta=wt*res*hh_theta
  ercix=matrix(rep(diag(1,true_n),J),true_n)
  ercix=t(ercix)%*%ercix
  Sig=(4/true_n)*t(h_theta)%*%ercix%*%h_theta
  
  element_a=apply(wt*res*(-pnorm((X-PSI)/hh)-((X-PSI)/hh)*dnorm((X-PSI)/hh)),2,sum)
  element_b=apply(wt*res*((BB/hh)*dnorm((X-PSI)/hh)*(2-((X-PSI)/hh)^2)),2,sum)
  
  latofHn=matrix(0,p+2*k,p+2*k)
  latofHn[(p+1):(p+k),(p+k+1):(p+2*k)]=diag(element_a,k)
  latofHn[(p+k+1):(p+2*k),(p+k+1):(p+2*k)]=diag(element_b,k)
  latofHn[(p+k+1):(p+2*k),(p+1):(p+k)]=latofHn[(p+1):(p+k),(p+k+1):(p+2*k)]
  Hn=(2/true_n)*(t(hh_theta)%*%diag(wt)%*%hh_theta-latofHn)
  
  
  A <- solve(Hn + diag(1e-8, p+2*k)) %*% Sig %*% solve(Hn + diag(1e-8, p+2*k))
  ese <- sqrt(diag(A)/true_n)
  
  
  bet.se <- ese[1:(p+k)]
  psi.se <- ese[(p+k+1):(p+2*k)]
  names(bet.se) <- names(bet.est)
  names(psi.se) <- names(psi.est)
  
  return(list(bet.est=bet.est, bet.se=bet.se,
              psi.est=psi.est,psi.se=psi.se))
}




