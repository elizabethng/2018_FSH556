#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_INTEGER( nt );
  DATA_VECTOR( y_t );

  // Parameters
  PARAMETER( x0 );
  PARAMETER( log_sigmaP );
  PARAMETER( log_sigmaM );
  PARAMETER( alpha );
  PARAMETER_VECTOR( x_t );
  
  // Objective funcction
  Type jnll = 0;
  
  // Probability of random coefficients
  jnll -= dnorm( x_t(0), x0, exp(log_sigmaP), true );
  for( int t=1; t<nt; t++){
    jnll -= dnorm( x_t(t), x_t(t-1) + alpha, exp(log_sigmaP), true ); // alpha is term for directed walk (drift term) and put in log_sigmaP so we can estimate on the unconstrained scale. also 1) closer to normal, so better to integrate out (if it were a random effect), 2) sample variance of normal dist (sum) is chi squared dist, and estimate might be appx ~ chi square dist m --> gives better est of CI?? 3) log maintains scale of params better (want them all to be on the same scale)
  }

  // Probability of data conditional on fixed and random effect values
  for( int t=0; t<nt; t++){
    jnll -= dnorm( y_t(t), x_t(t), exp(log_sigmaM), true );
  }
  
  // Reporting
  
  // Any benefit to REPORT outside of ADREPORT? speed...
  // also for debugging, you can get it for any point
  
  Type sigmaP = exp(log_sigmaP);
  Type sigmaM = exp(log_sigmaM);

  REPORT( sigmaP );
  REPORT( sigmaM );
  REPORT( x_t );

  ADREPORT( sigmaP );
  ADREPORT( sigmaM );
  ADREPORT( x_t );

  return jnll;
}


// Jim doesn't like ragged arrays... put in long form with a new indicator t_i which is the time interval for each obs
// 19-23 is same,but now
//  y_i(i) = x_t(t_i(i))