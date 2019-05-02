# include <Rcpp.h>
using namespace Rcpp ;

// [[Rcpp::export]] 
Rcpp::List splitContinuous(NumericVector zcut,
                           NumericMatrix zcutCat,
                           List datMatrix,
                           List parameters) {
  
  // Extract information from data matrix
  double useOtherNodes = as<double>(parameters["useOtherNodes"]);
  double isCtg = as<double>(parameters["isCtg"]);
  NumericVector y = as<NumericVector>(datMatrix["y"]);
  NumericVector x = as<NumericVector>(datMatrix["x"]);
  NumericVector prtx = as<NumericVector>(datMatrix["prtx"]);
  NumericVector status = as<NumericVector>(datMatrix["status"]);
  NumericVector kmCens = as<NumericVector>(datMatrix["kmCens"]);
  NumericVector trt = as<NumericVector>(datMatrix["trt"]);
  NumericVector trtNew = as<NumericVector>(datMatrix["trtNew"]);
  NumericVector inNode = as<NumericVector>(datMatrix["inNode"]);
  
  // Extract information from parameter list; 
  // Not Used yet bool isCtg = as<NumericVector>(parameters["isCtg"]);
  double nodeSize = as<double>(parameters["nodeSize"]);
  double trtSize = as<double>(parameters["trtSize"]);
  double maxScore = as<double>(parameters["maxScore"]);
  
  // Define several size parameters
  double totalSampSize = y.length();
  double zlength = zcutCat.ncol()*isCtg + zcut.length()*(1-isCtg);
  double ylength = y.length();
  double parentSize = 0;
  for(int i = 0; i < totalSampSize; i++){
    parentSize += inNode[i];
  }
  double divisor = totalSampSize*useOtherNodes + parentSize*(1-useOtherNodes);
  
  // =========================================================
  // Candidate splits sending trt = 1 to the left node
  // =========================================================
  
  NumericVector outputL(zlength);
  NumericVector outputR(zlength);
  NumericVector condA(zlength);
  NumericVector condB1(zlength);
  NumericVector condB2(zlength);
  
  // Determine the pass status for each candidate split 
  // condA: number of obs per daughter node 
  // condB1: enough treated in left daughter
  // condB2: enough treated in right daughter

  for(int i = 0; i < zlength; i++){
    for(int j = 0; j < ylength; j++){
      if(inNode[j] == 1){
        if(isCtg){
          condA[i] += (zcutCat(j,i));
          condB1[i] += (zcutCat(j,i))*trt[j];
          condB2[i] += (1-zcutCat(j,i))*trt[j];
        } else{
          condA[i] += (x[j] <= zcut[i]);
          condB1[i] += (x[j] <= zcut[i])*trt[j];
          condB2[i] += (x[j] > zcut[i])*trt[j];
        }
      }
    }
  }
  
  for(int i = 0; i < zlength; i++){
    if((condA[i] >= nodeSize) &&
       (parentSize - condA[i] >= nodeSize) &&
       (condB1[i] >= trtSize) &&
       (condB2[i] >= trtSize) &&
       (condA[i] - condB1[i] >= trtSize) &&
       (parentSize - condA[i] - condB2[i] >= trtSize)){
      for(int j = 0; j < ylength; j++){
        if(inNode[j] == 1){
          if(isCtg){
            outputL[i] += y[j]*(trt[j]*(zcutCat(j,i)) / (divisor*prtx[j]) + 
              (1 - trt[j])*(1-zcutCat(j,i)) / (divisor*(1 - prtx[j])));
            outputR[i] += y[j]*(trt[j]*(1-zcutCat(j,i)) / (divisor*prtx[j]) + 
              (1 - trt[j])*(zcutCat(j,i)) / (divisor*(1 - prtx[j])));
          } else{
            outputL[i] += y[j]*(trt[j]*(x[j] <= zcut[i]) / (divisor*prtx[j]) + 
              (1 - trt[j])*(x[j] > zcut[i]) / (divisor*(1 - prtx[j])));
            outputR[i] += y[j]*(trt[j]*(x[j] > zcut[i]) / (divisor*prtx[j]) + 
              (1 - trt[j])*(x[j] <= zcut[i]) / (divisor*(1 - prtx[j])));
          }
        } else{
          if(useOtherNodes){
            outputL[i] += y[j]*(trt[j]*(trtNew[j]) / (divisor*prtx[j]) + 
                    (1 - trt[j])*(1 - trtNew[j]) / (divisor*(1 - prtx[j])));
            outputR[i] += y[j]*(trt[j]*(trtNew[j]) / (divisor*prtx[j]) + 
                    (1 - trt[j])*(1 - trtNew[j]) / (divisor*(1 - prtx[j])));
          }
        }
      }
    } else{
      outputL[i] = -1E10;
      outputR[i] = -1E10;
    }
  }
  
  // Get output direction for splits
  Rcpp::StringVector outputDirection(zlength);
  NumericVector output(zlength);
  for(int i = 0; i < zlength; i++){
    if((outputL[i] > outputR[i]) && (outputL[i] > maxScore)){
      outputDirection[i] = "l";
      output[i] = outputL[i];
    } else if((outputL[i] < outputR[i]) && (outputR[i] > maxScore)){
      outputDirection[i] = "r";
      output[i] = outputR[i];
    } else{
      outputDirection[i] = "neither";
      output[i] = -1E10;
    }
  }
  
  // Compile output
  return List::create(
  _["direction"] = outputDirection, 
  _["output"] = output
  );  
}


