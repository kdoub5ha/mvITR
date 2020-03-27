#' @title Builds the random forest of interaction trees
#' 
#' @description This function constructs a random forest of ITR trees using either the IPWE or AIPWE
#' method. A forest object can be inputed into the `predict.ITR()` function along with data in order to 
#' obtain treatment predictions. A forest object can also be given to the `Variable.Importance.ITR()` 
#' function to estimate predictor importance. 
#' 
#' @param dat the data set being used to grow the random forest. Required input. 
#' @param col.y the response variable. Required input. 
#' @param col.trt the treatment indicator.  Must be binary. Required input.
#' @param col.prtx the probability of being assigned to treatment group. Required input. 
#' @param split.var vector of columns containing the desired splitting variables.  Required input. 
#' @param ctg identifies the categorical input columns.  Defaults to NA.  Not available yet. 
#' @param N0 minimum number of observations needed to call a node terminal.  Defaults to 20. 
#' @param n0 minimum number of treatment/control observations needed in a split to call a node terminal. Defaults to 5. 
#' @param test indicator that determines if testing data is also run with each tree in the forest.  Defaults to FALSE. 
#' @param max.depth controls the maximum depth of the tree. Defaults to 10. 
#' @param mtry sets the number of randomly selected splitting variables to be included. Defaults to max of length(split.var)/3 rounded down and 1.
#' @param AIPWE logical. Should AIPWE (TRUE) or IPWE (FALSE) be used. 
#' @param ntree sets the number of trees to be generated. Defaults to 500.
#' @param avoid.nul.tree controls if trees with no splits (null trees) are allowed. Defaults to FALSE.
#' @param stabilize logical. Should a numerical stabilization be used. Can be random forest ('rf') or linear model ('linear')
#' @param stabilize.type either 'rf' or 'linear' if stabilization requested.
#' @return A list of characteristics of the forest.
#' @return \item{ID.Boots.Samples}{list of bootstrap sample IDs}
#' @return \item{TREES}{list of trees}
#' @return \item{Model.Specification}{information about the input parameters of the forest}
#' @import randomForest
#' @export
#' @examples
#' dat <- gdataM(n=1000, depth=2, beta1=3, beta2=1)
#' # This builds a forest of 100 trees using the dataset called 'dat' with columns
#' # 'y', 'trt', and 'prtx' for the outcome, treatement indicator, and probability of being
#' # in treatment group, respectively.  The splitting variables are found in columns 1-4, 
#' # and we chose to avoid null trees.
#' forest<-Build.RF.ITR(dat=dat, col.y="y", col.trt="trt", col.prtx="prtx", 
#'                      split.var=1:4, ntree=100, avoid.nul.tree=TRUE)


Build.RF.ITR <- function(dat, 
                         split.var, 
                         test = NULL, 
                         col.y = "y",
                         col.r = "r",
                         col.trt = "trt", 
                         col.prtx = "prtx", 
                         ctg = NULL,
                         N0 = 20, 
                         n0 = 5,  
                         max.depth = 10,
                         ntree = 500, 
                         mtry = max(floor(length(split.var)/3), 1),
                         avoid.nul.tree = FALSE, 
                         AIPWE = FALSE, 
                         stabilize = TRUE, 
                         stabilize.type = c('linear', 'rf'), 
                         haoda.method = TRUE, 
                         haoda.ae.level = NA, 
                         verbose = FALSE, 
                         reverse.ae.scale = FALSE, 
                         use.other.nodes = TRUE, 
                         extremeRandomized = FALSE,
                         lambda = 0)
{
  require(randomForest)
  out <- as.list(NULL)
  out$ID.Boots.Samples  <- as.list(1:ntree)
  out$TREES <- as.list(1:ntree)
  out$preds.oob <- matrix(NA, nrow = nrow(dat), ncol = ntree)
  out$preds.inbag <- matrix(NA, nrow = nrow(dat), ncol = ntree)
  out$preds.cumulative.oob <- matrix(NA, nrow = nrow(dat), ncol = ntree)
  out$preds.cumulative.inbag <- matrix(NA, nrow = nrow(dat), ncol = ntree)
  out$value.oob <- rep(NA, ntree)
  out$value.inbag <- rep(NA, ntree)
  out$risk.oob <- rep(NA, ntree)
  out$risk.inbag <- rep(NA, ntree)
  out$risk.bound.values <- NULL
  
  # Initialize output if test set is included
  if(!is.null(test)) {
    out$test.preds <- matrix(NA, nrow = nrow(test), ncol = ntree)
    out$test.preds.cumulative <- matrix(NA, nrow = nrow(test), ncol = ntree)
    out$test.value <- rep(NA, ntree)
    out$test.risk <- rep(NA, ntree)
  }
  
  # set inputs
  stabilize.type <- match.arg(stabilize.type)
  if(haoda.method){
    if(is.na(haoda.ae.level)) stop("Risk allowance level must be specified numeric")
  }
  
  # set parameters for splitting criteria
  b <- 1
  while(b <= ntree){
    # TAKE BOOTSTRAP SAMPLES
    id.b <- sample(1:nrow(dat), size=nrow(dat), replace = TRUE)
    dat.b <- dat[id.b,]
    dat.test <- dat[-unique(id.b),]
    
    # Generate tree based on b-th bootstrap sample
    seed <- sample.int(100000, 1)
    set.seed(seed)
    
    tre.b <- grow.ITR(data = dat.b, test = test, min.ndsz = N0, n0 = n0, 
                      split.var = split.var, ctg = ctg, max.depth = max.depth, 
                      AIPWE = AIPWE, mtry = mtry,
                      efficacy = col.y, risk = col.r, stabilize = stabilize,
                      haoda.method = haoda.method, stabilize.type = stabilize.type,
                      haoda.ae.level = haoda.ae.level, 
                      use.other.nodes = use.other.nodes, 
                      extremeRandomized = extremeRandomized, 
                      lambda = lambda)
    
    if(avoid.nul.tree) {
      if(nrow(tre.b$tree) > 1) {
        out$ID.Boots.Samples[[b]] <- id.b
        out$TREES[[b]] <- tre.b$tree
        b <- b + 1
      }
    } else {
      out$ID.Boots.Samples[[b]] <- id.b
      out$TREES[[b]] <- tre.b$tree
      out$AEfit[[b]] <- tre.b$fit.ae
      if(nrow(tre.b$tree) > 1){
        preds <- predict.ITR(tre.b$tree, dat, split.var)$trt.pred
        inbag.idx <- seq_along(1:nrow(dat)) %in% unique(id.b)
        # tmp.risk.oob <- mean(dat[!inbag.idx,col.r] * (dat[!inbag.idx,col.trt] == preds[!inbag.idx]) / dat[!inbag.idx,col.prtx])
        # tmp.risk.inbag <- mean(dat[inbag.idx,col.r] * (dat[inbag.idx,col.trt] == preds[inbag.idx]) / dat[inbag.idx,col.prtx])
        
        out$preds.oob[,b] <- ifelse(inbag.idx, NA, preds)
        out$preds.inbag[,b] <- ifelse(inbag.idx, preds, NA)
      }
      out$preds.cumulative.oob[,b] <- ifelse(rowMeans(out$preds.oob[,1:b,drop=F], na.rm = T) > 0.5, 1, 0)
      out$preds.cumulative.inbag[,b] <- ifelse(rowMeans(out$preds.inbag[,1:b,drop=F], na.rm = T) > 0.5, 1, 0)
      
      if(b > 5){
        out$value.inbag[b] <- mean(dat[,col.y] * (out$preds.cumulative.inbag[,b] == dat[,col.trt]) / dat[,col.prtx], na.rm = T)
        out$value.oob[b] <- mean(dat[,col.y] * (out$preds.cumulative.oob[,b] == dat[,col.trt]) / dat[,col.prtx], na.rm = T)
        out$risk.inbag[b] <- mean(dat[,col.r] * (out$preds.cumulative.inbag[,b] == dat[,col.trt]) / dat[,col.prtx], na.rm = T)
        out$risk.oob[b] <- mean(dat[,col.r] * (out$preds.cumulative.oob[,b] == dat[,col.trt]) / dat[,col.prtx], na.rm = T)
      }
      
      if(!is.null(test)){
        if(nrow(tre.b$tree) > 1){
          out$test.preds[,b] <- predict.ITR(tre.b$tree, as.data.frame(test), split.var)$trt.pred
        }
        if(b > 5){
          out$test.preds.cumulative[,b] <- ifelse(rowMeans(out$test.preds[,1:b,drop=F], na.rm = T) > 0.5, 1, 0)
          out$test.value[b] <- mean(test[,col.y] * (out$test.preds.cumulative[,b] == test[,col.trt]) / test[,col.prtx], na.rm = T)
          out$test.risk[b] <- mean(test[,col.r] * (out$test.preds.cumulative[,b] == test[,col.trt]) / test[,col.prtx], na.rm = T)
        }
      }
      
      b <- b + 1
    }
    if(verbose){
      if(b %in% seq(0, ntree, 100)) print(paste0("On Tree ", b))
    }
  }
  
  Model.Specification <- as.list(NULL)
  Model.Specification$data <- dat
  Model.Specification$split.var <- split.var
  Model.Specification$ctg <- ctg
  Model.Specification$col.y <- col.y
  Model.Specification$col.r <- col.r
  Model.Specification$col.trt <- col.trt
  Model.Specification$col.prtx <- col.prtx
  Model.Specification$lambda <- lambda
  Model.Specification$haoda.ae.level <- haoda.ae.level
  out$Model.Specification <- Model.Specification
  return(out)
}