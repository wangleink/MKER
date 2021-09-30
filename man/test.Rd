\name{cterTest_lg}
\alias{cterTest_lg}
\title{test the existence of change point in the multi-kink expectile regression with longitudinal data}
\usage{
  cterTest_lg(y, x, z, tau, NB, J) or cterTest_lg(ID, y, x, z, tau, NB)
}
\arguments{
  \item{y}{A vector of response.}
  
  \item{x}{A scalar covariate with threshold.}
  
  \item{z}{A vector of covariates.}
  
  \item{tau}{the expectile level.}
  
  \item{NB}{resampling times.}
 
  \item{J/ID}{...}
}
\value{
  A list with the elements
  \item{Tn}{The statistic based on original data.}
  \item{Tn.NB}{The statistics by wild bootstrap.}
  \item{p.value}{The p-value by wild bootstrap.}
}
\description{
  This function for calculating the test statistics and p-value by wild bootstrap.
}
\examples{
  
 ## simulated data
 library(MASS)
 N=200;tau=0.5;J=4
 t0=1.5;bet0=c(1,3,0,1)
 
 x=runif(J*N,-5,5)
 z=rnorm(J*N,1,0.5^2)
 xz=cbind(1, x, pmax(x-t0, 0), z)
 
 rho=0.5;mu=matrix(rep(0,J));
 s_elem=NULL
 for(j in 1:J){
   ss_elem=NULL;for(i in 1:j){ss_elem=c(ss_elem,rho^(i-1))}
   s_elem=c(ss_elem,s_elem)
 }
 S=matrix(0,J,J);S[lower.tri(S,diag=TRUE)]=s_elem
 sigma=S+t(S)-diag(1,J,J)
 
 err=mvrnorm(N,mu,sigma);err=as.vector(err)
 y=xz \%*\% bet0+err
 
 fittest=cterTest_lg(y, x, z, tau = tau, NB = 200, J=4)
 fittest$p.value



 ## The example of NGHS dataset
 NGHS_An=NGHS[c(1:3358),c(1,3,6,9)];NGHS_An=na.omit(NGHS_An)
 tau=0.3
 y=NGHS_An$SBP
 x=NGHS_An$BMI
 z=NGHS_An$AGE
 ID=NGHS_An$ID
 fittest=cterTest_lg(ID, y, x, z, tau, NB=400)
 fittest$p.value
}
\keyword{test,test_real}

