#' @title Treatment Prediction Function
#'
#' @description Used to make treatment prediction for a single tree or random forest. If the 
#' input is a forest, then the proportion of trees voting for treatment (`trt=1`) is returned. 
#' If the input is a single tree, then the function returns the vote (treatment or control) from the tree. 
#' 
#' @param input tree or forest object from `grow.ITR` or `Build.RF.ITR`.
#' @param new.dat data for which predictions are desired
#' @param ctgs columns of categorical variables. 
#' @return A summary list of the following elements:
#' @return \item{SummaryTreat}{proportion of trees voting for treatment (trt=1). 
#' If input is a single tree then SummaryTreat is a single number. 
#' If input is a forest then SummaryTreat is a vector equal to the length of the number of trees.}
#' @return \item{trt.pred}{vector of treatment assignments {0, 1} based on the tree vote (single tree) or majority of tree votes (forest). This vector has length equal to the number of rows in `new.dat`.}
#' @return \item{n.trees}{number of tree in `input`}
#' @return \item{tree.votes}{matrix of votes for each tree for each subject in `new.dat`. Rows correspond to trees in `input` and columns correspond to subjects in `new.dat`.}
#' @return \item{data}{input data frame `new.dat`}
#' @return \item{NA.trees}{number of trees returning no votes. In a forest, this is the number of null trees.}
#' @import randomForest
#' @export
#' @examples
#' dat <- gdataM(n=1000, depth=2, beta1=3, beta2=1)
#' # Build a forest with 100 trees
#' forest <- Build.RF.ITR(dat, col.y="y", col.trt="trt", col.prtx="prtx", split.var=1:4, ntree=100)
#' # Predict treatment assignments for 1000 observations in `dat`
#' predict.ITR(forest, dat)



predict.ITR <- function(input, new.dat, split.var, ctgs = NULL){
  if(is.null(dim(input))){
    trees <- input$TREES
    n.trees <- length(trees)
  } else{
    trees <- input
    n.trees <- 1
  }
  dat <- new.dat
  n <- nrow(dat)
  out <- NULL
 
  result <- sapply(1:n.trees, function(i){
    if(is.null(dim(input))){
      tre <- trees[[i]]
    } else{
      tre <- trees
    }
    
    if(nrow(tre) > 0){
      if(!is.na(tre[1,6])){
        idx <- !is.na(tre$cut.2)
        cutPoint <- as.numeric(tre$cut.2[idx])
        splitVar <- as.integer(tre$var[idx])
        treNodes <- as.character(tre$node[idx])
        direction <- as.character(tre$cut.1[idx])
        Data <- as.matrix(dat[,split.var,drop=F])
        send <- SendDown(cutPoint, splitVar, Data, treNodes, direction)
        trt.pred <- send$trt.pred
      } else{
        trt.pred <- rep(NA, n)
      }
    } else{
      trt.pred <- rep(NA, n)
    }
    
    return(trt.pred) 
  })
  
  if(!(is.null(dim(result)) & length(result) == 1)){
    out$SummaryTreat <- apply(result, 1, FUN = mean, na.rm=T)
  } else{
    out$SummaryTreat <- result
  }
  if(is.null(dim(input))){
    out$trt.pred <- ifelse(out$SummaryTreat < 0.5, 0, 1)
  } else{
    out$trt.pred <- out$SummaryTreat
  }

  out$n.trees <- n.trees
  out$tree.votes <- result
  out$data <- new.dat
  out$NA.trees <- ifelse((!(is.null(dim(result)) & length(result) == 1)), sum(is.na(result[1,])), NA)
  return(out)
}