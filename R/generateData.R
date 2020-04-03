#' @title Generates example data
#' @description Generates n example data samples for constructing and rcDT model.
#' @param n numeric. Number of example samples to be generated.
#' @return Summary of a single interaction tree. Each `node` begins with "0" indicating the root node, 
#' followed by a "1" or "2" indicating the less than (or left) child node or greater than (or right) child node. 
#' Additionally, the number of observations `size`, number treated `n.1`, number on control `n.0`, and treatment effect `trt.effect`
#' summaries are provided.  The splitting information includes the column of the chosen splitting variable `var`, the variable name 'vname',
#' the direction the treatment is sent `cut.1` ("r" for right child node, and "l" for left), the chosen split value `cut.2`, 
#' and the estimated value function `score`.
#' @import randomForest
#' @export
#' @examples
#' dat <- generateData()

generateData <- function(n = 1000){
  trt <- sample(c(-1,1), n, replace = TRUE)
  X <- matrix(round(runif(10*n), 3), n, 10, dimnames = list(c(1:n), paste0("x", 1:10)))
  subgroup <- (X[,1] <= 0.7 & X[,2] <= 0.7)
  Y <- 1 - X[,3] + X[,4] + 3*subgroup*(trt == 1) + 3*(1-subgroup)*(trt == -1) + rnorm(n)
  R <- 2 + X[,4] + 2*X[,3] + (2*X[,2] + 0.1)*trt + rnorm(n, sd = 0.25)
  out <- data.frame(X, 
                    trt = ifelse(trt == -1, 0, 1), 
                    y = round(Y, 4), 
                    r = round(R, 4), 
                    prtx = 0.5)
  return(out)
}