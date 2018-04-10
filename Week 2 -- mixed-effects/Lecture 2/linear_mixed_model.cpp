#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_INTEGER( n_groups );
  DATA_IVECTOR( g_i ); 
  DATA_VECTOR( y_i );
  
  // Parameters
  PARAMETER( beta0 );    // global intercept
  PARAMETER( log_SD0 );  // estimate unconstrained and exponetiate later
  PARAMETER( log_SDZ );  
  PARAMETER_VECTOR( z_g ); // random effects?
  
  // Objective funcction
  Type jnll = 0;
  int n_i = y_i.size();

// Code up joint neg log like (no integral here...)
  // Probability of data conditional on fixed and random effect values
  for( int i=0; i<n_i; i++){
    jnll -= dnorm( y_i(i), beta0 + z_g(g_i(i)), exp(log_SD0), true );
  }
  
  // Probability of random coefficients
    // different from linear model b/c we need to 
    // account for hyperparameter, theta2 (SDZ)
  for( int g=0; g<n_groups; g++){
    jnll -= dnorm( z_g(g), Type(0.0), exp(log_SDZ), true );
  }
  
  // Reporting
  Type SDZ = exp(log_SDZ);
  Type SD0 = exp(log_SD0);

  REPORT( SDZ );
  REPORT( SD0 );
  ADREPORT( SDZ );
  ADREPORT( SD0 );

  return jnll;
}
