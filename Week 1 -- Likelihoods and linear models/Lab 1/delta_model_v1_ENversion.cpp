
#include <TMB.hpp>

// dlnorm
// template<class Type> //need this at top of any function
// Type dlognorm(Type x, Type meanlog, Type sdlog, int give_log=0){  
	// defining a function in tmb, same syntax as dlnorm. scalar return
	//return 1/(sqrt(2*M_PI)*sd)*exp(-.5*pow((x-mean)/sd,2));
  // Type logres = dnorm( log(x), meanlog, sdlog, true) - log(x); //internal calculation for the Jacobian of transformation
  // if(give_log) return logres; else return exp(logres);
//}

// main function
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_VECTOR( y_i );  // data are a vector, match type of R object (c++ macro, code written in tmb as a way to input)
  DATA_MATRIX( X_ij );

  // Parameters
  // need different params for tweedie
  // read in restricted values here
	// PARAMETER_VECTOR( b_j );
	// PARAMETER_VECTOR( theta_z );
  PARAMETER( log_mu );  // >= 0
  PARAMETER( log_phi ); // needs to be > 0
  PARAMETER( t_p );   // in (1, 2) use logistic transformation and add 1
  

  // Objective function
  // do the parameter transformations to get unrestricted values
	// Type zero_prob = 1 / (1 + exp(-theta_z(0)));
	// Type logsd = exp(theta_z(1));
  Type mu = exp(log_mu);
  Type phi = exp(log_phi);
  Type p = 1 + 1 / (1 + exp(-t_p));
  
  Type jnll = 0;
  int n_data = y_i.size();

  // Linear predictor
	// just using a constant mu for now
	// vector<Type> linpred_i( n_data ); // elements of our vector are type Type
	// linpred_i = X_ij * b_j;

  // Probability of data conditional on fixed effect values
  for( int i=0; i<n_data; i++){
     // if(y_i(i)==0) jnll -= log( zero_prob );
     // if(y_i(i)!=0) jnll -= log( 1-zero_prob ) + dlognorm( y_i(i), linpred_i(i), logsd, true );
	 jnll -= dtweedie(y_i(i), mu, phi, p, true); // check that this is the right use of dtweedie
  }
  
  // Reporting
  // these are tmb-specific macros
  ADREPORT( mu );
  ADREPORT( phi );
  ADREPORT( p );
  // ADREPORT( zero_prob ); //way to tell it to calc sd of something. zero_prob is the logistic transform

  return jnll;
}
