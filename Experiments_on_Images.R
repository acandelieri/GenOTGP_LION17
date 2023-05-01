rm(list=ls()); graphics.off(); cat("\014")

library(transport)
library(DiceKriging)
library(doParallel)


# *****************************************************
# Technical parameters
# *****************************************************

set.seed(42)

dgt = 6 # number of digits for grey-scale color
n1 = 16 # number of horizontal pixels
n2 = 16 # number of vertical pixels

source.digit = 3
target.digit = 1
N = 300 # number of training images

rescale.input = T
rescale.predictions = F

kernel.type = "gauss"
nugget.estimation = T




# *****************************************************
# Loading data
# *****************************************************

cat("> Loading data...\n")

dataset = read.delim( file="C:\\Users/Antonio Candelieri/Downloads/usps",
                      header=F, sep=" ")

# Data pre-processing

# - selecting only the first 500 examples for each digit
# - replacing "digit id" to be from 0 to 9 (instead of 1 to 10)
# - removing useless columns
# - fomratting each image as a vector
cat("> Pre-processing data...\n")
ds = NULL
for( i in 1:10 ) {
  ixs = which( dataset$V1 == i )
  tmp = dataset[ ixs[1:500], -258 ]
  tmp$V1 = rep(i-1,500)
  for( j in 2:ncol(tmp) )
    tmp[,j] = matrix(as.numeric(unlist(strsplit(tmp[,j],":",T))),ncol=2,byrow=T)[,2]
  ds = rbind( ds, tmp )
}
ds = as.matrix(ds)
rm(tmp)


if( rescale.input ) {
  cat("> Rescaling color for each image...\n")
  for( i in 1:nrow(ds) )
    ds[i,-1] = round( ( ds[i,-1] - min(ds[i,-1]) )/( max(ds[i,-1]) - min(ds[i,-1]) ), dgt )
}


# Separating between training and test set
if( N < ncol(ds)-1 ) {
  cat(" * N must be at least ",ncol(ds),". It has been modified accordingly!\n")
  N = ncol(ds)
}

trn.data = tst.data = NULL
trn.label = tst.label = NULL
for( i in 1:10 ) {
  r1 = ((i-1)*500) + 1:N
  r2 = ((i-1)*500) + (N+1):500
  trn.data = rbind( trn.data, ds[r1,-1] )
  tst.data = rbind( tst.data, ds[r2,-1] )
  trn.label = c( trn.label, ds[r1,1] )
  tst.label = c( tst.label, ds[r2,1] )
}
rm(ds)


# selecting source and target
X = trn.data[which(trn.label==source.digit),]
Y = trn.data[which(trn.label==target.digit),]

X.tst = tst.data[which(tst.label==source.digit),]
# Y.tst = tst.data[which(tst.label==target.digit),] # it is not relevant!




# generating the ground truth

cat("> Generating the ground-truth...\n")

res = transport( pp(X), pp(Y), p=2, method="primaldual" )
res = res[order(res$from),]
I = res$to


Y_ = Y[I,]


# Learning the multi-GPs model

cat("> Learning",n1,"x",n2,"=",n1*n2,"GPs...\n")

cores = max(1,detectCores()-1)
myCluster <- makeCluster( cores, type = "PSOCK") # PSOCK type on Win!
registerDoParallel(myCluster)

started = Sys.time()
cat(" * Started at",toString(started),"\n")
gps = foreach( i = 1:ncol(Y_), .packages="DiceKriging" ) %dopar% {
  gp = NULL
  # trials = 0
  # while( is.null(gp) & trials<2 ) {
    # try( expr=(
      # gp = km( design = data.frame(X), response = Y_[,i], optim.method="BFGS",
      #          covtype = kernel.type, nugget.estim = nugget.estimation, iso=T,
      #          #upper=0.5, lower=10^dgt,
      #          control=list(trace=F) )
      # ), silent=T )
  #   if( is.null(gp) )
  #     trials = trials+1
  # }
  # if( is.null(gp) )
  #   stop("IMPOSSIBLE TO LEARN A GP!")
  
  if( length(unique(Y_[,i]))>1 ) {
    gp = km( design = data.frame(X), response = Y_[,i], optim.method="BFGS",
             covtype = kernel.type, nugget.estim = nugget.estimation, iso=T,
             #upper=0.5, lower=10^dgt,
             control=list(trace=F) )
  } else {
    gp = Y_[1,i]
  }
}
finished = Sys.time()
cat(" * Finished at",toString(finished),"\n")
cat(" * Elapsed:",difftime(finished,started,units="secs"),"[secs]\n")
cat(" * with",cores,"cores\n\n")
stopCluster(myCluster)




# *****************************************************
# Testing
# *****************************************************

preds = NULL

for( gp in gps ) {
  if( class(gp)=="km" ) {
    pred = predict( gp, data.frame(X.tst), "UK" )
    p = pred$mean    
    if( rescale.predictions )
      p = (p-min(p))/(max(p)-min(p))
  } else {
    p = rep(gp,nrow(X.tst))
  }
  
  preds = cbind( preds, round(p,dgt) )
  # TODO: weak?
}

gc()


# for( i in 1:nrow(X.tst) ) {
png( paste0(source.digit,"_to_",target.digit,".png"),
     width=800, height=800 )
par(mfrow=c(20,20))
par(mar=rep(0.1,4))
for( i in 1:200 ) {

  image( matrix(X.tst[i,],n1,n2,byrow=F)[,n2:1],col=rev(grey(0:10^dgt/10^dgt)),
       xaxs="r",yaxs="r",xaxt="n",yaxt="n")
  # image( matrix(X.tst[i,],n1,n2,byrow=F)[,n2:1],col=colorRampPalette(c("white", "skyblue"))(10),
  #        xaxs="r",yaxs="r",xaxt="n",yaxt="n")
  image( matrix(preds[i,],n1,n2,byrow=F)[,n2:1],col=rev(grey(0:10^dgt/10^dgt)),
        xaxs="r",yaxs="r",xaxt="n",yaxt="n")
}
dev.off()
