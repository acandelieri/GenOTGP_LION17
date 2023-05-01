ot = transport( pp(X_), pp(Ymo.m), p=2 )
IJ.mo = ot$to
mis.mo = length(which(1:nrow(X_)!=IJ.mo))

ot = transport( pp(X_), pp(Yso.m), p=2 )
IJ.so = ot$to
mis.so = length(which(1:nrow(X_)!=IJ.so))

cat("MISMATCH MO-GP:\t", round(100*mis.mo/nrow(X_),2),"%\n" )
cat("MISMATCH SO-GPs:", round(100*mis.so/nrow(X_),2), "%\n" )

estC.mo = actC.mo = estC.so = actC.so = 0

for( i in 1:nrow(X_) ) {
  estC.mo = estC.mo + sum((X_[i,] - Ymo.m[i])^2)
  actC.mo = actC.mo + sum((X_[i,] - Ymo.m[IJ.mo[i]])^2)
  estC.so = estC.so + sum((X_[i,] - Yso.m[i])^2)
  actC.so = actC.so + sum((X_[i,] - Yso.m[IJ.so[i]])^2)
}
cat("Delta Cost [%] MO-GP: ",round(100*(estC.mo-actC.mo)/actC.mo,2),"%\n")
cat("Delta Cost [%] SO-GPs:",round(100*(estC.so-actC.so)/actC.so,2),"%\n")






# Coherence of the generated transported point cloud!
NN = 5000
cat("\n\n> Generating a very large point cloud from 'Q'...\n")
if( use.case == 1 ) {
  params = toy2D_G2G()
  YY = GM.sampling( n=NN, d=2, k=1, means=params$mu.Q, covs=params$Sigma.Q )
} else {
  if( use.case == 2 ) {
    params = toy2D_G2GM()
    YY = GM.sampling( n=NN, d=2, k=length(params$mu.Q), means=params$mu.Q, covs=params$Sigma.Q )
  } else {
    # use.case == 3
    params = toy2D_G2SR()
    YY = SR.sampling( n=NN, SR.coeff1=params$SR.coeff1, SR.coeff2=params$SR.coeff2,
                     l1=params$l1, l2=params$l2, m1=params$m1, m2=params$m2,
                     s1=params$s1, s2=params$s2 )
  }
}


curr.mar = par("mar")
par(mar=c(3.1,2.1,2.1,1.1))

par(mfrow=c(2,2))

plot( YY[,1], YY[,2], pch=21, bg="grey", col="grey", cex=1.5,
      xlab="", ylab="", xlim=c(-2,2), ylim=c(-2,2),
      main="qualitative 'coherence' MO-GP (strong)")
grid()
points( Ymo.m[,1], Ymo.m[,2], pch=21, bg="green3", cex=1.5 )

plot( YY[,1], YY[,2], pch=21, bg="grey", col="grey", cex=1.5,
      xlab="", ylab="", xlim=c(-2,2), ylim=c(-2,2),
      main="qualitative 'coherence' indep. GPs (strong)")
grid()
points( Yso.m[,1], Yso.m[,2], pch=21, bg="green3", cex=1.5 )

plot( YY[,1], YY[,2], pch=21, bg="grey", col="grey", cex=1.5,
      xlab="", ylab="", xlim=c(-2,2), ylim=c(-2,2),
      main="qualitative 'coherence' MO-GPs (weak)")
grid()
for( i in 1:nrow(X_) )
  DrawEllipse( Ymo.m[i,1], Ymo.m[i,2], Ymo.s[i,1], Ymo.s[i,2],
               col=adjustcolor("green3",alpha.f=0.1), border=NA )

plot( YY[,1], YY[,2], pch=21, bg="grey", col="grey", cex=1.5,
      xlab="", ylab="", xlim=c(-2,2), ylim=c(-2,2),
      main="qualitative 'coherence' indep. GPs (weak)")
grid()
for( i in 1:nrow(X_) )
  DrawEllipse( Yso.m[i,1], Yso.m[i,2], Yso.s[i,1], Yso.s[i,2],
               col=adjustcolor("green3",alpha.f=0.1), border=NA )


par(mfrow=c(1,1))
par(mar=curr.mar)
