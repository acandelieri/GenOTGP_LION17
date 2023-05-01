rm(list=ls()); graphics.off(); cat("\014")

library(transport)
library(DiceKriging)

source("use_cases.R")


# *****************************************************
# Input
# *****************************************************

# set.seed(42)

cat("> Select the use case:\n")
cat(" [1] Gaussian to Gaussian\n")
cat(" [2] Gaussian to Gaussian Mixture (8 components)\n")
cat(" [3] Gaussian to Swiss Roll\n\n")

usecase = ""
while( !usecase %in% c("1","2","3") ) {
  usecase = readline("> Your choice: ")
}

N = as.integer(readline("> Size of the TRAINING point-clouds (suggested 500): "))
M = as.integer(readline("> Size of the SOURCE point-cloud to use for generation (suggested >100): "))

trn.vis = ""
while( !(tolower(trn.vis) %in% c("n","y")) )
  trn.vis = readline("> Visualize the OT map on training data? [y/n]: ")
trn.vis = ifelse( tolower(trn.vis)=="y", T, F )
  

kernel.type = "exp"
nugget.estimation = T


# *****************************************************
# Training
# *****************************************************

if( usecase=="1" ) {
  
  cat("\n***** Gaussian to Gaussian *****\n")
  toyExample = toy2D_G2G()
  
  cat("> Sampling a SOURCE point-cloud for training...\n")
  X = GM.sampling(n=N,d=2,k=1,means=toyExample$mu.P,covs=toyExample$Sigma.P)
  cat("> Sampling a TARGET point-cloud for training...\n")
  Y = GM.sampling(n=N,d=2,k=1,means=toyExample$mu.Q,covs=toyExample$Sigma.Q)
  
  cat("\n> Sampling a NEW SOURCE point-cloud to use for generation...\n")
  X.val = GM.sampling(n=M,d=2,k=1,means=toyExample$mu.P,covs=toyExample$Sigma.P)
  
} else {
  
  if( usecase=="2" ) {
    
    cat("\n***** Gaussian to Gaussian Mixture (8 components) *****\n\n")
    toyExample = toy2D_G2GM()
    
    cat("> Sampling a SOURCE point-cloud for training...\n")
    X = GM.sampling(n=N,d=2,k=length(toyExample$mu.P),means=toyExample$mu.P,covs=toyExample$Sigma.P)
    cat("> Sampling a TARGET point-cloud for training...\n")
    Y = GM.sampling(n=N,d=2,k=length(toyExample$mu.Q),means=toyExample$mu.Q,covs=toyExample$Sigma.Q)
    
    cat("\n> Sampling a NEW SOURCE point-cloud to use for generation...\n")
    X.val = GM.sampling(n=M,d=2,k=length(toyExample$mu.P),means=toyExample$mu.P,covs=toyExample$Sigma.P)
    
  } else {
    
    cat("\n***** Gaussian to Swiss Roll *****\n\n")
    toyExample = toy2D_G2SR()
    
    cat("> Sampling a SOURCE point-cloud for training...\n")
    X = GM.sampling(n=N,d=2,k=length(toyExample$mu.P),means=toyExample$mu.P,covs=toyExample$Sigma.P)
    cat("> Sampling a TARGET point-cloud for training...\n")
    Y = SR.sampling(n=N,SR.coeff1=toyExample$SR.coeff1,SR.coeff2=toyExample$SR.coeff2,
                    l1=toyExample$l1,l2=toyExample$l2,m1=toyExample$m1,m2=toyExample$m2,
                    s1=toyExample$s1,toyExample$s2) 
    
    cat("\n> Sampling a NEW SOURCE point-cloud to use for generation...\n")
    X.val = GM.sampling(n=M,d=2,k=length(toyExample$mu.P),means=toyExample$mu.P,covs=toyExample$Sigma.P)
    
  }  
}


cat("\n> Solving the Optimal Transport problem via primal-dual...\n")
ot = transport( pp(X), pp(Y), p=2, method="primaldual" )


IJ = ot$to[order(ot$from)]

Y_ = Y[IJ,]

cat("> Generalizing the optimal transport map via GP regression...\n")
gps = MGP.sogps(X=X,Y=Y_,kernels=rep(kernel.type,ncol(Y)),nugget.estim=nugget.estimation )


if( trn.vis ) {
  plot( X[,1], X[,2], pch=19, cex=1.5, xlim=c(-2,2), ylim=c(-2,2) )
  for( i in 1:nrow(X) )
    lines( c(X[i,1], Y[i,1]), c(X[i,2],Y[i,2]), col="grey" )
  points( X[,1], X[,2], pch=19, col="pink", cex=1.5 )
  points( Y[,1], Y[,2], pch=19, col="green3", cex=1.5 )
  grid()
}


Y.mu = Y.sd = NULL
for( i in 1:length(gps)){
  pred = predict(gps[[i]],X,"UK")
  Y.mu = cbind( Y.mu, pred$mean )
  Y.sd = cbind( Y.sd, pred$sd )
}

# points( Y.mu[,1], Y.mu[,2], pch=19, col="blue", cex=1.5 )


ot = transport( pp(X), pp(Y.mu), p=2 )
IJ.pred = ot$to[order(ot$from)]
mismatch = length( which(IJ.pred != 1:nrow(X)) )

cat("\n> Training statistics:\n")
cat(" - MISMATCH: ", round(100*mismatch/length(IJ),2),"%\n", sep="" )
cat("   (RMSE: ",round(mean( sqrt(apply( (Y_ - Y.mu)^2, 1, sum )) ),8),")\n",sep="")

estCost = actCost = 0

for( i in 1:nrow(X) ) {
  estCost = estCost + sum((X[i,] - Y.mu[i])^2)
  actCost = actCost + sum((X[i,] - Y.mu[IJ.pred[i]])^2)
}
cat(" - Delta Cost [%]: ",round(100*(estCost-actCost)/actCost,2),"%\n", sep="")




# *****************************************************
# Generation
# *****************************************************

cat("\n> Using GPs to generate a TARGET point-cloud from the sampled NEW SOURCE point-cloud...\n")
Y.mu = Y.sd = NULL
for( i in 1:length(gps)){
  pred = predict(gps[[i]],X.val,"UK")
  Y.mu = cbind( Y.mu, pred$mean )
  Y.sd = cbind( Y.sd, pred$sd )
}

cat("\n> Computing the possible mismatch...\n")
ot = transport( pp(X.val), pp(Y.mu), p=2 )
IJ.pred = ot$to[order(ot$from)]
mismatch = length( which(IJ.pred != 1:nrow(X.val)) )

cat("\n> Generating statistics:\n")
cat(" - MISMATCH: ", round(100*mismatch/length(IJ.pred),2),"%\n", sep="" )
cat("   (RMSE: ",round(mean( sqrt(apply( (Y.mu - Y.mu[IJ.pred,])^2, 1, sum )) ),8),")\n",sep="")

for( i in 1:nrow(X) ) {
  estCost = estCost + sum((X.val[i,] - Y.mu[i])^2)
  actCost = actCost + sum((X.val[i,] - Y.mu[IJ.pred[i]])^2)
}
cat(" - Delta Cost [%]: ",round(100*(estCost-actCost)/actCost,2),"%\n", sep="")


plot( X.val[,1], X.val[,2], pch=21, cex=1.5, bg="pink", col="black", xlim=c(-2,2), ylim=c(-2,2) )
grid()

plot( X.val[,1], X.val[,2], pch=19, cex=1.5, xlim=c(-2,2), ylim=c(-2,2) )
for( i in 1:nrow(X.val) )
  lines( c(X.val[i,1], Y.mu[i,1]), c(X.val[i,2],Y.mu[i,2]), col="grey" )
points( X.val[,1], X.val[,2], pch=21, bg="pink", col="black", cex=1.5 )
points( Y.mu[,1], Y.mu[,2], pch=21, bg="green3",col="black", cex=1.5 )
grid()


if( usecase=="1") {
  YY = GM.sampling(n=10000,d=2,k=1,means=toyExample$mu.Q,covs=toyExample$Sigma.Q)
} else {
  if( usecase=="2") {
    YY = GM.sampling(n=10000,d=2,k=length(toyExample$mu.Q),means=toyExample$mu.Q,covs=toyExample$Sigma.Q)
  } else {
    YY = SR.sampling(n=10000,SR.coeff1=toyExample$SR.coeff1,SR.coeff2=toyExample$SR.coeff2,
                     l1=toyExample$l1,l2=toyExample$l2,m1=toyExample$m1,m2=toyExample$m2,
                     s1=toyExample$s1,toyExample$s2) 
  }
}
plot( YY[,1], YY[,2], pch=19, col=adjustcolor("grey",alpha.f=0.3),
      cex=1.5, xlim=c(-2,2), ylim=c(-2,2) )
points( Y.mu[,1], Y.mu[,2], pch=21, bg="green3", col="black", cex=1.5 )
grid()

