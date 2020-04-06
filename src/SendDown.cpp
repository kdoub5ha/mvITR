# include <Rcpp.h>
using namespace Rcpp ;

// [[Rcpp::export]] 
Rcpp::List SendDown(NumericVector cutPoint,
                    IntegerVector splitVar,
                    NumericMatrix Data, 
                    StringVector treNodes,
                    StringVector direction) {
  
  // Extract information from data matrix
  int varLength = splitVar.length();
  int n = Data.nrow();
  StringVector nodeOut(n);
  for(int i = 0; i < n; i++){
    nodeOut[i] = "0";
  }
  NumericVector leftNode(n);
  NumericVector inNodeTmp(n);
  NumericVector tmpData(n);
  IntegerVector trtOut(n);
  StringVector tmpDir(n);
  StringVector lastDir(n);

  for(int j = 0; j < varLength; j++){ // j is the index for each splitting variable
    NumericVector tmpData = Data(_, splitVar[j]-1);

    for(int i = 0; i < n; i++){
      leftNode[i] = (tmpData[i] <= cutPoint[j]);
      inNodeTmp[i] = (nodeOut[i] == treNodes[j]);
      if(inNodeTmp[i]){
        lastDir[i] = direction[j];
        if(leftNode[i]){
          nodeOut[i] += "1";
          tmpDir[i] = "1";
        } else{
          nodeOut[i] += "2";
          tmpDir[i] = "2";
        }
      }
    }
  }

  for(int i = 0; i < n; i++){
    if(tmpDir[i] == "1"){
      if(lastDir[i] == "l"){
        trtOut[i] = 1;
      } else{
        trtOut[i] = 0;
      }
    } else{
      if(lastDir[i] == "l"){
        trtOut[i] = 0;
      } else{
        trtOut[i] = 1;
      }
    }
  }
  
  // Compile output
  return List::create(
    _["node"] = nodeOut,
    _["trt.pred"] = trtOut
    );
}
