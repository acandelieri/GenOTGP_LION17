library(MGMM)
library(DiceKriging)
library(mlegp)


# Sampling n datapoints from a d-dimensional Mixture of k Gaussians
GM.sampling <- function( n, d, k, means, covs ) {
  #TODO: stopifnot
  GM.sampling = rGMM( n=n, d=d, k=k, means=means, covs=covs )
}


# Sampling n datapoints from a 2-dimensional Swiss Roll
SR.sampling <- function( n, SR.coeff1, SR.coeff2, l1=1, l2=1, m1=0, m2=0, s1=1, s2=1 ) {
  a = runif( n, SR.coeff1*pi, SR.coeff2*pi )
  SR.sampling = cbind( a*cos(a)/l1 + rnorm(length(a),m1,s1), a*sin(a)/l2 + rnorm(length(a),m1,s1) )
}


# Learning a multi-output GP (via the 'mlegp' package)
MGP.mogp <- function( X, Y, ... ) {
  #TODO: stopifnot
  MGP.mgop = mlegp( X, Y, ... )
}


# Learning a multi-output GP (as a list of independent single-output GPS,
# via the 'DiceKriging' package)
MGP.sogps <- function( X, Y, kernels=rep("gauss",ncol(Y)), nugget.estim=T ) {
  #TODO: stopifnot + kernels
  gps = list()
  for( i in 1:ncol(Y) ) {
    gps[[i]] = km( design=data.frame(X), response=Y[,i], covtype=kernels[i],
                   nugget.estim=nugget.estim, control=list(trace=0) )
    cat("> * GP #",i,": trained!\n",sep="")
  }
  MGP.sogps = gps
}



toy2D_G2G <- function() {

  d=2
  
  mu.P=numeric(d); Sigma.P=(0.5)^2 * diag(d)
  
  mu.Q=numeric(d); Sigma.Q=0.5*diag(d)
  
  toy2D_G2G = list( mu.P=mu.P, Sigma.P=Sigma.P, mu.Q=mu.Q, Sigma.Q=Sigma.Q )
}

toy2D_G2GM <- function() {
  
  d=2 
  k=8
  r=sqrt(2)
  
  mu.P=numeric(d); Sigma.P=(0.5)^2*diag(d)
  
  mu.Q=Sigma.Q=list()
  for( i in 1:k ) {
    mu.Q[[i]] = c( r * sin((i-1)*pi/4), r * cos((i-1)*pi/4) )
    Sigma.Q[[i]] = (0.05)^2 * diag(d)
  }
  
  toy2D_G2GM = list( mu.P=mu.P, Sigma.P=Sigma.P, mu.Q=mu.Q, Sigma.Q=Sigma.Q )
}

toy2D_G2SR <- function() {
  
  d=2 
  mu.P = numeric(d); Sigma.P = (0.5)^2 * diag(d)
  
  SR.coeff1 = 1.5
  SR.coeff2 = 4.5
  l1 = l2 = 7
  m1 = m2 = 0
  s1 = s2 = 0.1
  
  toy2D_G2SR = list( mu.P=mu.P, Sigma.P=Sigma.P,
                     SR.coeff1=SR.coeff1, SR.coeff2=SR.coeff2,
                     l1=l1, l2=l2, m1=m1, m2=m2, s1=s1, s2=s2 )
}