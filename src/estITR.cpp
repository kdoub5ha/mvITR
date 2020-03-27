#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double estITR(List input) {
  // Extract information from data matrix
  NumericVector y = as<NumericVector>(input["y"]);
  NumericVector ae = as<NumericVector>(input["ae"]);
  NumericVector prtx = as<NumericVector>(input["prtx"]);
  NumericVector kmCens = as<NumericVector>(input["KM.cens"]);
  IntegerVector trt = as<IntegerVector>(input["trt"]);
  IntegerVector status = as<IntegerVector>(input["status"]);
  IntegerVector z = as<IntegerVector>(input["z"]);
  int n0 = as<int>(input["n0"]);
  double lambda = as<double>(input["lambda"]);
  double maxRisk = as<double>(input["maxRisk"]);
  // int aug = as<int>(input["aug"]);
  
  // Set up augmented 
  // if(aug){
  //   double first = 0;
  //   double second = 0;
  // }
  // 
  
  // Condition checks
  int check1 = 0;
  int check2 = 0;
  int check3 = 0;
  for(int i = 0; i < y.length(); i++){
    check1 += z[i];
    check2 += trt[i];
    check3 += 1-trt[i];
  }
  
  double itrOut1 = 0;
  
  if((check1 >= n0) && (check2 >= n0) && (check3 >= n0)){
    for(int i = 0; i < y.length(); i++){
      itrOut1 += (y[i] - lambda * ae[i]) * (((trt[i] * z[i] * status[i]) / (prtx[i] * kmCens[i])) + 
        (((1-trt[i]) * (1-z[i]) * status[i]) / ((1-prtx[i]) * kmCens[i]))) / y.length();
    }
  } else{
    itrOut1 = -1E10;
  }
  
  // double output(1);
  // output = itrOut1;
  
  // Return ITR
  return itrOut1 + lambda * maxRisk;
}