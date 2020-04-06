# include <Rcpp.h>
using namespace Rcpp ;

// [[Rcpp::export]] 
List estOpt(NumericVector y,
              NumericVector r,
              NumericVector trt,
              NumericVector prtx,
              NumericVector rule,
              double tau,
              double lambda){
  
  double eff = 0;
  double risk = 0;
  int ylength = y.length();
  int rlength = r.length();
  
  for(int i = 0; i < ylength; i++){
    eff += y[i] * (trt[i] * rule[i] + (1-trt[i]) * (1-rule[i])) / prtx[i];
  }
  for(int i = 0; i < rlength; i++){
    risk += r[i] * (trt[i] * rule[i] + (1-trt[i]) * (1-rule[i])) / prtx[i];
  }
  
  double output = eff / ylength - lambda * (risk / rlength - tau);
  return List::create(
    _["optValue"] = output,
    _["efficacy"] = eff / ylength, 
    _["risk"] = risk / rlength
    );
}
