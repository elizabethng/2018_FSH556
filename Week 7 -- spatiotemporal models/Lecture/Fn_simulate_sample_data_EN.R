logmean=1
Scale=0.2
SD_omega=1
SD_epsilon=1
n_per_year=100
n_years=10

# logmean=1; Scale=0.2; SD_omega=1; SD_epsilon=1; n_per_year=100; n_years=10
Sim_Fn = function( logmean=1, Scale=0.2, SD_omega=1, SD_epsilon=1, n_per_year=100, n_years=10 ){

  # Generate n_per_year uniform random locs on 1x1 square
  # but we'll look at the same locations in each n_years
  loc_xy = cbind( "x"=runif(n_per_year), "y"=runif(n_per_year)) 

  # Set up the model distributions for GRF
  RF_omega = RMgauss(var=SD_omega^2, scale=Scale)
  RF_epsilon = RMgauss(var=SD_epsilon^2, scale=Scale)

  # Matrix for data observations?
  log_d_rt = Epsilon_rt = matrix(NA, ncol=n_years, nrow=n_per_year)
  
  # RFsimulate: simulates random field conditional on locations passed
  #  then extract the values and use column selection to unlist
  # these are the spatial (main) effects
  Omega_r = RFsimulate(model=RF_omega, x=loc_xy[,'x'], y=loc_xy[,'y'])@data[,1]
  
  # fill in the data matrix of observations
  for(t in 1:n_years){
    # simulate the temporal effects for the first year
    Epsilon_rt[,t] = RFsimulate(model=RF_epsilon, x=loc_xy[,'x'], y=loc_xy[,'y'])@data[,1]
    
    # linear predictor (spatial effects are constant in time)
    log_d_rt[,t] = logmean + Epsilon_rt[,t] + Omega_r
  }
  
  # plot the temporal trend only
  test = data.frame(x = loc_xy[, "x"],
                    y = loc_xy[, "y"])
  test = cbind(test, Epsilon_rt)
  test = melt(test, id.vars = c("x", "y"))
  ggplot(test, aes(x, y, color = value)) + geom_point() + facet_wrap(~variable)
  
  # plot true density trends
  test = data.frame(x = loc_xy[, "x"],
                    y = loc_xy[, "y"])
  test = cbind(test, log_d_rt)
  test = melt(test, id.vars = c("x", "y"))
  ggplot(test, aes(x, y, color = value)) + geom_point() + facet_wrap(~variable)
  
  # True biomass (??) in each year, given spatial and temporal variation
  b_t = colSums( exp(log_d_rt) )
  ts.plot(b_t)

  # Simulate samples for each site and year
  DF = expand.grid("s_i"=1:n_per_year, "t_i"=1:10)
  DF = cbind( DF, "log_d_rt"=log_d_rt[ as.matrix(DF[,c('s_i','t_i')]) ] ) # True log densities
  DF = cbind( DF, "c_i"=rpois(n_per_year*n_years, lambda=exp(DF[,'log_d_rt'])) ) # poisson distributed counts are observed

  # Return stuff
  Return = list("DF"=DF, "Epsilon_rt"=Epsilon_rt, "Omega_r"=Omega_r, "loc_xy"=loc_xy, "b_t"=b_t, "n_per_year"=n_per_year, "n_years"=n_years)
  return( Return )
}
