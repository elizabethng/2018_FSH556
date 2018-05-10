
# setwd( "C:/Users/James.Thorson/Desktop/Project_git/2018_FSH556/Week 7 -- spatiotemporal models/Lecture/")
setwd("~/GitHub/2018_FSH556/Week 7 -- spatiotemporal models/Lecture")

# Libraries and functions
library( TMB )
library( INLA )
library( RandomFields )
library( RANN )

# Simulate data
source( "Fn_simulate_sample_data.R")
SimList = Sim_Fn( logmean=1, Scale=0.2, SD_omega=1, SD_epsilon=1, n_per_year=100, n_years=10 )

# Make triangulated mesh
mesh = inla.mesh.create( SimList$loc_xy )
spde = inla.spde2.matern(mesh)

# Area for each location
# Make an arbitarily fine grid
# use a kmeans grid to find out how many grid cells are closest to location
# then that number is a good proxy for the area
# helps us avoid dealing with the function itself
# we pick 100 points within the domain to model (knots)
# More efficient than a Monte Carlo estimator
loc_extrapolation = expand.grid( "x"=seq(0,1,length=1e3), "y"=seq(0,1,length=1e3) )
# Find the nearest neighbors
NN_extrapolation = nn2( data=SimList$loc_xy, query=loc_extrapolation, k=1 )
# count the number of assigned locations and calc area as a proportion of total locations
a_s = table(factor(NN_extrapolation$nn.idx,levels=1:nrow(SimList$loc_xy))) / nrow(loc_extrapolation)

# Compile
Version = "spatial_index_model"
compile( paste0(Version,".cpp") )
dyn.load( dynlib(Version) )

# Make inputs
Data = list("n_s"=SimList$n_per_year,  # number obs per year
            "n_t"=SimList$n_years,     # number years
            "a_s"=a_s,                 # area of each location
            "c_i"=SimList$DF$c_i,      # observed counts
            "s_i"=SimList$DF$s_i-1,    # location for observation i (indexed from 0)
            "t_i"=SimList$DF$t_i-1,    # time interval for observation i (indexed from 0)
            "M0"=spde$param.inla$M0,
            "M1"=spde$param.inla$M1,   # Precision matrix components from INLA 
            "M2"=spde$param.inla$M2)

Params = list("beta0"=0, 
              "ln_tau_O"=log(1), 
              "ln_tau_E"=log(1), 
              "ln_kappa"=1, 
              "omega_s"=rep(0,mesh$n), 
              "epsilon_st"=matrix(0,nrow=mesh$n,ncol=Data$n_t))

Random = c("omega_s", "epsilon_st")
# annoying but there is a constant intercept through time (global intercept)
#   could have had one for each year for consistent notation
# epsilon_st is the space time effects (100 locations plus 16 buffers)
# tau controls sd of omega, which is spatial effect
# kappa is decorrelation rate for both space and time variance
# including bias.correct = TRUE to deal with retransformation bias
#   ranef in log space, so smoother will be biased when we do another transformation

# Build and run
Obj = MakeADFun( data=Data, parameters=Params, random=Random )
Opt = TMBhelper::Optimize( obj=Obj, newtonsteps=1, getsd=TRUE, bias.correct=TRUE )
Report = Obj$report()
unlist( Report[c('Range','SigmaO','SigmaE')] )

# Confirm that index is uncorrelated among years
Cov_Dt = Opt$SD$cov
Corr_Dt = cov2cor( Cov_Dt )
# Pull out the covariance
# off diagonals are all less than 0.01, since we've desinged the model to have minimal covariance
# between years

summary( Opt$SD, "report" ) # <-- handy!!!


plot(x = SimList$b_t, y = Opt$SD$value) # true vs estiamted, looks plausible, good index even though scale differes
matplot( x=1:SimList$n_years, 
         y=cbind(SimList$b_t/mean(SimList$b_t),Opt$SD$value/mean(Opt$SD$value)), 
         col=c("black","red"), type="l", lwd=2 )



