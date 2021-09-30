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
# ## The example of NGHS dataset
# NGHS_An=NGHS[c(1:3358),c(1,3,6,9)];NGHS_An=na.omit(NGHS_An)
# tau=0.3
# y=NGHS_An$SBP
# x=NGHS_An$BMI
# z=NGHS_An$AGE
# ID=NGHS_An$ID
# fittest=cterTest_lg(ID, y, x, z, tau, NB=400)
# fittest$p.value




cterTest_lg <- function(ID, y, x, z, tau, NB){
  
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
  testFun <- function(ID, y, x, z, tau, tt){
    
    ## under H0
    ## by asymmetric weight least square
    X <- cbind(1, x, z)
    fit <- expTLlaws(X, y, tau, max.iter, tol)
    bet <- fit$coefficients
    res <- as.vector(y - X %*% bet)
    wt <- wfun(res, tau)
    
    n <- length(unique(ID))      ####
    
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
  testFun.resample <- function(ID, y, x, z, tau, tt){
    
    #########################
    ## permutation random errors
    n<- length(unique(ID))       ####
    u <- rnorm(n, 0, 1)
    uu=NULL; for(i in 1:n){uu[which(ID==i)]=u[i]}
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
      
      Sn.t <- apply(xz*(wt*(x-tt[kk])*ifelse(x <= tt[kk], 1, 0)), 2, sum)/n #S1n
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
  Tn <-  testFun(ID, y, x, z, tau, tt)
  Tn.NB <- replicate(NB, testFun.resample(ID, y, x, z, tau, tt))
  
  pv <- mean(Tn.NB > Tn,  na.rm = TRUE)
  
  return(list(Tn = Tn, Tn.NB = Tn.NB, p.value = pv))
  
  
}

