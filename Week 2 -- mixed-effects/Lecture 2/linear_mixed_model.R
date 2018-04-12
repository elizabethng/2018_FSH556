
# setwd( "C:/Users/James.Thorson/Desktop/Project_git/2018_FSH556/Week 2 -- mixed-effects/Lecture 2" )
setwd("~/Documents/GitHub/2018_FSH556/Week 2 -- mixed-effects/Lecture 2")
Use_REML = FALSE

devtools::install_github("kaskr/TMB_contrib_R/TMBhelper")

######################
# Simulate data
######################

# Simulate predictors
group_i = rep( 1:10, each=10)
z_g = rnorm( length(unique(group_i)), mean=0, sd=1)
beta0 = 0

# Simulate response
y_i = z_g[group_i] + beta0 + rnorm( length(group_i), mean=0, sd=1)

######################
# Run in R
######################

library(lme4)
Lme = lmer( y_i ~ 1|factor(group_i), REML=Use_REML)

# partitions variance

######################
# Run in TMB
######################

library(TMB)

# Compile model
Version = "linear_mixed_model"
compile( paste0(Version,".cpp") )

# SDz is sd between groups
# SD0 is overall ??

# Build inputs
Data = list( "n_groups"=length(unique(group_i)), "g_i"=group_i-1, "y_i"=y_i)
Parameters = list( "beta0"=-10, "log_SD0"=2, "log_SDZ"=2, "z_g"=rep(0,Data$n_groups) )
Random = c("z_g") ## THIS IS HOW TMB KNOWS TO INTEGRATE AND DO LAPLACE APPX (TMB even takes this out of the Obj$par call)
if( Use_REML==TRUE ) Random = union( Random, "beta0")

# Laplace appx to neg marginal log likelihood
# first run will show the inner and outer parts (inner optimizer is newton step)
# once it's zero, it then does laplace appx (the outer part??) Nope, no outer call here actually
Obj$fn(Obj$par)
Obj$gr(Obj$par)  # restarts from previous best value. smart! outputs gradient wrt to each fixed effect

# Build object
dyn.load( dynlib("linear_mixed_model") )
Obj = MakeADFun(data=Data, parameters=Parameters, random=Random)  #

# Prove that function and gradient calls work
Obj$fn( Obj$par )
Obj$gr( Obj$par )

# Optimize
# Iterates between inner and outer until the final gradient is very small
Opt = TMBhelper::Optimize( obj=Obj, newtonsteps=1 )

# TMBhelper automatically gets the SEs for the fixed effects (inputted params), 
# which comes from the inverse hessian matrix (diagonal) using the laplace appx
# (asymptotic result)
# (inverse?) hessian matrix is an estimator for Sigma, the asymptotic variacne-covariance
# matrix for the MVN distribution of the MLEs

# e.g., predictive distribution would be something like  ~ MVN (MLEs, inv.hessian)
Opt$SD$cov.fixed
sqrt(diag(Opt$SD$cov.fixed))
Opt$SD

# Output using delta method
summary(Opt$SD, "report")


# Get reporting and SEs
Report = Obj$report()
ParHat = as.list( Opt$SD, "Estimate" )

######################
# Shrinkage estimator
######################

# jim's attempt at replicating this from first principles
Mu = mean(y_i)
  Mu_s = tapply( y_i, INDEX=group_i, FUN=mean)
Sigma = sd( Mu_s )
  Sigma_s = sd( y_i - Mu_s[group_i] )
Weights_hat = c( 1/Sigma^2, length(y_i)/length(unique(group_i))/Sigma_s^2 )
  Weights_hat = Weights_hat / sum(Weights_hat)

# Predictions
Mu_s_hat = ( Mu*Weights_hat[1] + Mu_s*Weights_hat[2] )


######################
# Compare estimates
######################

# Global mean
c( fixef(Lme), ParHat$beta0, Mu )

# Random effects
cbind( "True"=z_g, "Lme4"=ranef(Lme)[['factor(group_i)']], "TMB"=ParHat$z_g, "Shrinkage_estimator"=Mu_s-Mu )

# Variances
summary(Lme)
unlist( Report[c("SDZ","SD0")] )

# TMB gives identical answer to R, and very close to shrinkage estimator
# helps prove to ourselves that we know what we're doing