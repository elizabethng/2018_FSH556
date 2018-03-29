
# setwd( "C:/Users/James.Thorson/Desktop/Project_git/2018_FSH556/Week 1 -- Likelihoods and linear models/Lab 1" )
library(TMB)

###########
# Nonlinear optimization
###########

run_animations = FALSE

if( run_animations==TRUE ){
  library("animation")  # Requires v2.12 or higher, also need to download: http://www.imagemagick.org/script/binary-releases.php

  # Nonlinear function
  RosenbrookFn = function(Params, Write=TRUE){
    Dev = (1-Params[1])^2 + 100*(Params[2]-Params[1]^2)^2
    if(Write==TRUE){
      if("Trace.txt" %in% list.files()){
        write.table(matrix(c(Dev,Params),nrow=1), file="Trace.txt", append=TRUE,col.names=FALSE,row.names=FALSE)
      }else{
        write.table(matrix(c(Dev,Params),nrow=1), file="Trace.txt", append=FALSE,col.names=FALSE,row.names=FALSE)
      }
    }
    return(Dev)
  }

  # Generate contour
  X = seq(-10,20,length.out=1e3)
  Y = seq(-10,20,length.out=1e3)
  Z = matrix(NA, ncol=length(X), nrow=length(Y))
  for(xI in 1:length(X)){
  for(yI in 1:length(Y)){
    Z[xI,yI] = RosenbrookFn(c(X[xI],Y[yI]), Write=FALSE)
  }}

  # Rosenbrook function
  png( file="Rosenbrook.png", width=6, height=6, res=200, units="in")
    par( mar=c(3,3,2,0) )
    contour(x=X, y=Y, z=Z, levels=seq(0,1000,by=2)^5+1)
    points( x=1, y=1, cex=2, col="red", lwd=2)
  dev.off()

  # Plot path -- Nelder Mead
  if("Trace.txt" %in% list.files()) file.remove("Trace.txt")
  Opt = nlminb( start=c(4,4), objective=RosenbrookFn, Write=TRUE )
    Trace = read.table("Trace.txt" )
  contour(x=X, y=Y, z=Z, levels=seq(0,1000,by=10)^3+1)
    lines(x=Trace[,2], y=Trace[,3],col="red")
    ani.options(interval=0.10, nmax=1e4)
    saveVideo(expr={for(i in 1:nrow(Trace)){contour(x=X, y=Y, z=Z, levels=seq(0,1000,by=10)^3+1); lines(x=Trace[max(1,i-10):i,2], y=Trace[max(1,i-10):i,3],col="red",lwd=3); points(x=c(1,Trace[i,2]),y=c(1,Trace[i,3]), col=c("black","red"),cex=2,lwd=3)}}, video.name="Nelder-Mead.mp4", convert = 'gm convert', clean=TRUE)         #

  # Plot path -- BFGS
  if("Trace.txt" %in% list.files()) file.remove("Trace.txt")
  Opt = optim( par=c(4,4), fn=RosenbrookFn, Write=TRUE )
    Trace = read.table("Trace.txt" )
  contour(x=X, y=Y, z=Z, levels=seq(0,1000,by=10)^3+1)
    lines(x=Trace[,2], y=Trace[,3],col="red")
    ani.options(interval=0.10, nmax=1e4)
    saveVideo(expr={for(i in 1:nrow(Trace)){contour(x=X, y=Y, z=Z, levels=seq(0,1000,by=10)^3+1); lines(x=Trace[max(1,i-10):i,2], y=Trace[max(1,i-10):i,3],col="red",lwd=3); points(x=c(1,Trace[i,2]),y=c(1,Trace[i,3]), col=c("black","red"),cex=2,lwd=3)}}, video.name="BFGS.mp4", convert = 'gm convert', clean=TRUE)         #

  # Plot path -- TMB
  # Step 1 -- make and compile template file
  compile( "Rosenbrook.cpp" )

  # Step 2 -- build inputs and object
  dyn.load( dynlib("Rosenbrook") )
  Params = list( "Params"=c(4,4) )
  Data = list( "dummy"=0 )
  Obj = MakeADFun( data=Data, parameters=Params, DLL="Rosenbrook")

  # Redefine functions
  Obj$fn_trace <- function( Params, ...){
    Dev = Obj$fn( Params, ... )
    if("Trace.txt" %in% list.files()){
      write.table(matrix(c(Dev,Params),nrow=1), file="Trace.txt", append=TRUE,col.names=FALSE,row.names=FALSE)
    }else{
      write.table(matrix(c(Dev,Params),nrow=1), file="Trace.txt", append=FALSE,col.names=FALSE,row.names=FALSE)
    }
    return( Dev )
  }

  # Step 3 -- test and optimize
  if("Trace.txt" %in% list.files()) file.remove("Trace.txt")
  Opt = nlminb( start=Obj$par, objective=Obj$fn_trace, gradient=Obj$gr )
    Trace = read.table("Trace.txt" )
  contour(x=X, y=Y, z=Z, levels=seq(0,1000,by=10)^3+1)
    lines(x=Trace[,2], y=Trace[,3],col="red")
    ani.options(interval=0.10, nmax=1e4)
    saveVideo(expr={for(i in 1:nrow(Trace)){contour(x=X, y=Y, z=Z, levels=seq(0,1000,by=10)^3+1); lines(x=Trace[max(1,i-10):i,2], y=Trace[max(1,i-10):i,3],col="red",lwd=3); points(x=c(1,Trace[i,2]),y=c(1,Trace[i,3]), col=c("black","red"),cex=2,lwd=3)}}, video.name="TMB.mp4", convert = 'gm convert', clean=TRUE)         #
}

###########
# Delta-model for canary rockfish
###########

devtools::install_github("nwfsc-assess/geostatistical_delta-GLMM") # data are in here
devtools::install_github("kaskr/TMB_contrib_R/TMBhelper") # package to ease optimization process
library( SpatialDeltaGLMM )

#
data(WCGBTS_Canary_example)
CPUE = WCGBTS_Canary_example$HAUL_WT_KG

hist(CPUE)
mean(CPUE == 0) # 92% zero obs
mean(CPUE[which(CPUE>0)]) # but a few have a lot ("dust bin" distribution)
max(CPUE)
# Distribution for data? can't use lognormal or gamma because there are zeros
# This is biomass (non-integer) sampled
# Tuesday we fit a normal distirbution (eq to least squares) but we have non-neg stuff...
# Not clear that there's a great distribution

X = cbind( "Intercept"=rep(1,length(CPUE)) ) # matrix with 1 col and n obs rows

# Step 1 -- make and compile template file
compile( "delta_model_v1.cpp" )

# Step 2 -- build inputs and object
dyn.load( dynlib("delta_model_v1") )
Params = list("log_mu"=0, 
              "log_phi"=0,
              "t_p" = 0) # must match in tmb file?
Data = list( "y_i"= CPUE, "X_ij"= X )
Obj = MakeADFun( data=Data, parameters=Params, DLL="delta_model_v1")
# simplest delta model with three paranms, start values at zero

# Step 3 -- test and optimize
Obj$fn( Obj$par )
Obj$gr( Obj$par ) # gradient, from overloaded use of Type
Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr ) # gradient declines over course of call
Opt$diagnostics = data.frame( "name"=names(Obj$par), "Est"=Opt$par, "final_gradient"=as.vector(Obj$gr(Opt$par)))
# added slot of diagnostics, final gradients are wihin our tolerance
Opt$par # estimated parameters
SD = sdreport( Obj ) # standard errors using delta method
args(Obj$report)
Obj$env$last.par # tmb keeps track of the values?? 
# use Reports to keep track of intermediate values
Opt$par # use for point estimates

### Check convergence
# Are the gradients close to zero
all(abs(Opt$diagnostics[,'final_gradient'])>0.0001) # tells about good convergiance
# Is the hessian positive definite?
(SD$pdHess==TRUE) # yay, hessian is positive definite


# Extract stuff
Report = Obj$report(par = Opt$par)
# Report = SD$par.fixed
otherreport = sdreport(Obj)
summary(SD) # to print the vlaues
names(SD)

# use dyn.unload() to unload the model
# ADREPORT macro in tmb part gives back plogis(2.54) ... does the delta method internally to calculate the sd for a derived quantity (somtimes harder in R)

# Visualize fit
png( file="Canary_histogram--with_fit.png", width=4, height=4, res=200, units="in")
  par( mar=c(3,3,2,0), mgp=c(2,0.5,0), tck=-0.02)
  hist( log(1+CPUE), freq=FALSE, col=rgb(1,0,0,0.2) )
  Sim_CPUE = (1-rbinom(1e5, size=1, prob=Report$zero_prob)) * rlnorm(1e5, meanlog=Report$linpred_i, sdlog=Report$logsd)
  hist( log(1+Sim_CPUE), freq=FALSE, add=TRUE, col=rgb(0,0,1,0.2) )
  legend( "topright", bty="n", legend=c("Observed","Predicted"), fill=c("red","blue"))
dev.off()

# So fit seems like it's fitting ok!

###########
# Crossvalidation experiment using delta-model for canary rockfish
###########

# assign every datum an integer from 1 to 10
# for first crossvalidation loop, pass the same, but now 
# we have predTF_i which pulls out the data selection for the loop
# 

# Step 0 -- make and compile template file
compile( "delta_model_v2.cpp" )
dyn.load( dynlib("delta_model_v2") )

# Step 1 -- divide into partitions
K = 10
Partition_i = sample( x=1:K, size=length(CPUE), replace=TRUE )
PredNLL_k = rep(NA, K)

# Step 2 --Loop through partitions
for(k in 1:K){
  Params = list("b_j"=rep(0,ncol(X)), "theta_z"=c(0,0))
  Data = list( "y_i"=CPUE, "X_ij"=X, predTF_i=ifelse(Partition_i==k,1,0) )
  Obj = MakeADFun( data=Data, parameters=Params, DLL="delta_model_v2")

  # Optimize
  Opt = TMBhelper::Optimize( obj=Obj, newtonsteps=1, getsd=FALSE ) # gets final gradient to be better
  SD = sdreport( Obj ) # standard errors
  final_gradient = Obj$gr( Opt$par )

  # Check convergence
  if( any(abs(final_gradient)>0.001) | SD$pdHess==FALSE ) stop("Not converged")

  # Report stuff
  Report = Obj$report()
  PredNLL_k[k] = Report$pred_jnll
}

# log-Predictive probability per datum
mean( PredNLL_k / table(Partition_i) )

# Best predictive log-likelihood with matching models

