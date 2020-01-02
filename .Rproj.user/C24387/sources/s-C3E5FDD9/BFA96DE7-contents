#include <Rcpp.h>
using namespace Rcpp;

//' @title Random walk Metropolis sampler using Rcpp.
//' @description Implement a random walk Metropolis sampler for generating the standard Laplace distribution
//' @param x0 initial point
//' @param sigma normal distribution variance
//' @param n sample size
//' @return x a random sample of size \code{n} 
//' @return k the number of iterations
//' @examples
//' \dontrun{
//' x0<-25
//' sigma<-0.05
//' n<-2000
//' rwMC(x0,sigma,n)
//' }
//' @export
// [[Rcpp::export]]

List rwMC(double x0, double sigma, int N){
  NumericVector x(N);
  x[0]=x0;
  int k=0;
  for(int i=1;i<N;i++){
    
    double y=as<double>(rnorm(1,x[i-1], sigma));
    double u=as<double>(runif(1));
    
    if (u<=exp(-abs(y))/exp(-abs(x[i-1]))) {
      x[i]=y;
      k=k+1;  
    }
    else x[i]=x[i-1];
  }
  
  return(List::create(Named("x")=x,Named("k")=k));
}