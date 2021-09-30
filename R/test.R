#'
#' test the existence of change point in the multi-kink expectile regression with longitudinal data
#'
#' This function for calculating the test statistics and p-value by wild bootstrap.
#'
#' @rdname cterTest_lg
#' @param y A vector of response
#' @param x A scalar covariate with threshold
#' @param z A vector of covariates
#' @param tau the expectile level
#' @param NB resampling times
#' @param J/ID ...
#'
#' @return A list with the elements
#' \item{Tn}{The statistic based on original data.}
#' \item{Tn.NB}{The statistics by wild bootstrap.}
#' \item{p.value}{The p-value by wild bootstrap.}
#'
# ## simulated data
# library(MASS)
# N=200;tau=0.5;J=4
# t0=1.5;bet0=c(1,3,0,1)
# 
# x=runif(J*N,-5,5)
# z=rnorm(J*N,1,0.5^2)
# xz=cbind(1, x, pmax(x-t0, 0), z)
# 
# rho=0.5;mu=matrix(rep(0,J));
#for(i in 1:J)
#{for(j in 1:J)
#{sigma[i,j]=rho^abs(i-j)}}
# 
# err=mvrnorm(N,mu,sigma);err=as.vector(err)
# y=xz %*% bet0+err
# 
# fittest=cterTest_lg(y, x, z, tau = tau, NB = 200, J=4)
# fittest$p.value




cterTest_lg <- function(y, x, z, tau, NB, J){
  
  ## global variable
  max.iter = 100
  tol = 1e-5
  
  expTLlaws <- function(X, y, tau,  max.iter, tol)
  {
    # X is the model matrix
    # y is the response vector of observed proportion
    # maxIter is the maximum number of iterations
    # tol is a convergence criterion
    #X <- cbind(1, X) # add constant
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
    
    ## the loss function
    loss <- sum(w*(y-c(X%*%b))^2)
    
    # ## the variance
    # Am <- t(X) %*% diag(w) %*% X
    # if(min(eigen(Am)$values)<1e-8){Am=as.matrix(nearPD(Am)$mat)}
    # 
    # A <- solve(Am)
    # H <- w * diag(X %*% A %*% t(X))
    # B <- t(X) %*% diag(w^2*(y-ypred)^2/(1-H)) %*% X
    # Vb <- A %*% B %*% A  # variance function
    list(coefficients=b, loss=loss, iterations=it)
  }
  
  
  wfun <- function(u, tau){
    abs(tau- (u <= 0))
  }
  
  ### test statistic based on origin data
  testFun <- function(y, x, z, tau, tt, J){
    
    ## under H0
    ## by asymmetric weight least square
    X <- cbind(1, x, z)
    fit <- expTLlaws(X, y, tau, max.iter, tol)
    bet <- fit$coefficients
    res <- as.vector(y - X %*% bet)
    wt <- wfun(res, tau)
    
    n <- length(y)/J
    
    Rn <-  rep(0, length(tt))
    for (kk in 1:length(tt)){
      
      Rn[kk] <- 1/sqrt(n)*sum(
        wt*res*(x-tt[kk])*ifelse(x<=tt[kk], 1, 0)
      )
      
    }
    
    Tn <- max(abs(Rn))
    
    return(Tn)
  }
  
  ## perturbed method to calculate the p-value
  testFun.resample <- function(y, x, z, tau, tt, J){
    
    #########################
    ## permutation random errors
    n<- length(y)/J
    u <- rnorm(n, 0, 1)
    uu=rep(u,J)
    #########################
    
    X <- cbind(1, x, z)
    fit <- expTLlaws(X, y, tau, max.iter, tol)
    bet <- fit$coefficients
    res <- as.vector(y - X %*% bet)
    wt <- wfun(res, tau)
    
    ## Sn
    xz <-  cbind(1, x, z)
    Sn <- (t(xz) %*% diag(wt) %*% xz)/n #Swn
    
    
    ## under H0
    Rn <-  rep(0, length(tt))
    for (kk in 1:length(tt)){
      
      Sn.t <- J*apply(xz*(wt*(x-tt[kk])*ifelse(x <= tt[kk], 1, 0)), 2, mean) #S1n
      Rn[kk] <- 1/sqrt(n)*sum(
        uu * wt *res *((x-tt[kk])*ifelse(x<= tt[kk],1,0) -
                         xz %*% solve(Sn) %*% Sn.t)
      )
      
      
    }
    Tn <- max(abs(Rn))
    
    return(Tn)
  }
  
  
  #######################################################
  ###  calculate the p-value by wild bootstrap
  
  tt <- seq(min(x), max(x), length = 100)
  Tn <-  testFun(y, x, z, tau, tt, J)
  Tn.NB <- replicate(NB, testFun.resample(y, x, z, tau, tt, J))
  
  pv <- mean(Tn.NB > Tn,  na.rm = TRUE)
  
  return(list(Tn = Tn, Tn.NB = Tn.NB, p.value = pv))
  
  
}

