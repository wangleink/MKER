\name{mker.bea}
\alias{mker.bea}
\title{Fit the multi-kink expectile regression with longitudinal data}
\usage{
  mker.bea(y,thre.x,cont.z,tau,J,Cn) or mker.bea(ID,y,thre.x,cont.z,tau,Cn)
}
\arguments{
  \item{y}{A vector of response.}
  \item{thre.x}{A scalar covariate with threshold effect.}
  \item{cont.z}{A vector of covariates with constant slopes.}
  \item{tau}{The quantile level that belongs to (0,1).}
  \item{Cn}{A psotive number corresponding to different types of BIC.}
  \item{J/ID}{...}
  \item{control}{A list returned by \code{fit.control}.}
}
\value{
  A list with the elements
   \item{bet.est}{The estimated regression coefficients with intercept.}
   \item{bet.se}{The estimated standard error of the regression coefficients.}
   \item{psi.est}{The estimated change points.}
   \item{psi.se}{The estiamted standard errors of threshold parameters.}
   \item{n.psi}{The estiamted number of change points.}
}
\description{
  A expectile BIC combining the backward elimination algorithm to estimate the multi-kink expectile regression.
}
\examples{
  
## simulated data
library(MASS)
N=200;tau=0.5;J=4
psi0=c(-2,1,3);k=length(psi0)
PSI=matrix(rep(psi0,rep(J*N,k)),ncol=k)
bet0=c(1,1,1,-3,4,4)

x=runif(J*N,-5,5)
z=rnorm(J*N,1,0.5^2)
X=matrix(rep(x,k),nrow=J*N)
XZ=cbind(1,z,x,pmax((X-PSI),0))

rho=0.5;mu=matrix(rep(0,J));
s_elem=NULL
for(j in 1:J){
  ss_elem=NULL;for(i in 1:j){ss_elem=c(ss_elem,rho^(i-1))}
  s_elem=c(ss_elem,s_elem)
}
S=matrix(0,J,J);S[lower.tri(S,diag=TRUE)]=s_elem
sigma=S+t(S)-diag(1,J,J)
#homoscedastic normal errors
L=0;err=mvrnorm(N,mu,sigma);err=as.vector(err)-efun(tau,errtype = 1)

y=XZ\%*\%bet0+(1+L*z)*err

fit=mker.bea(y,x,z,tau,J,Cn=(log(N)))



## The example of NGHS dataset
NGHS_An=NGHS[c(1:3358),c(1,3,6,9)];NGHS_An=na.omit(NGHS_An)
tau=0.3
y=NGHS_An$SBP
x=NGHS_An$BMI
z=NGHS_An$AGE
ID=NGHS_An$ID
n=length(unique(ID))
fit=mker.bea(ID,y,x,z,tau,Cn=log(n),control=fit.control(K.max=4))

}
\keyword{mker,mker_real}

