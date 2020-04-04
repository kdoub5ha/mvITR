#' @title Cross Validation for Optimal rcDT Model Selection for Given Penalty (lambda)
#' @description Performs k-fold cross validation for rcDT model to 
#' select the best subtree from the set of optimally pruned subtree generated from 
#' `prune` function. 
#' @param dat data.frame. Data used to construct rcDT model.  
#' Must contain efficacy variable (y), 
#' risk variable (r), 
#' binary treatment indicator coded as 0 / 1 (trt), 
#' propensity score (prtx),
#' candidate splitting covariates.
#' @param split.var numeric vector. Columns of spliting variables.
#' @param efficacy char. Efficacy outcome column. Defaults to 'y'.
#' @param risk char. Risk outcome column. Defaults to 'r'.
#' @param col.trt char. Treatment indicator column name. Should be of form 0/1 or -1/+1.
#' @param col.prtx char. Propensity score column name. 
#' @param risk.control logical. Should risk be controlled? Defaults to TRUE.
#' @param risk.threshold numeric. Desired level of risk control. 
#' @param lambda numeric. Penalty parameter for risk scores. Defaults to 0, i.e. no constraint.
#' @param test data.frame of testing observations. Should be formatted the same as 'data'.
#' @param N0 numeric specifying minimum number of observations required to call a node terminal. Defaults to 20.
#' @param n0 numeric specifying minimum number of treatment/control observations needed in a split to declare a node terminal. Defaults to 5. 
#' @param max.depth numeric specifying maximum depth of the tree. Defaults to 15 levels. 
#' @param mtry numeric specifying the number of randomly selected splitting variables to be included. Defaults to number of splitting variables.
#' @param stabilize logical indicating if efficacy should be modeled using residuals. Defaults to TRUE. 
#' @param stabilize.type character specifying method used for estimating residuals. Current options are 'linear' for linear model (default) and 'rf' for random forest. 
#' @param sort internal use.
#' @param use.other.nodes logical. Should global estimator of objective function be used. Defaults to TRUE. 
#' @param use.bootstrap logical. Should a bootstrap resampling be done? Defaults to FALSE.
#' @param ctg numeric vector corresponding to the categorical input columns.  Defaults to NULL.  Not available yet. 
#' @param AIPWE logical. Should AIPWE (TRUE) or IPWE (FALSE) be used. Not available yet. 
#' @param extremeRandomized logical. Experimental for randomly selecting cutpoints in a random forest model. Defaults to FALSE and users should change this at their own peril. #' @return A summary of the cross validation including optimal penalty parameter and the optimal model. 
#' @return \item{best.tree.size}{optimal rcDT model based on size}
#' @return \item{best.tree.alpha}{optimal rcDT model based on alpha parameter selection}
#' @return \item{best.alpha}{optimal lambda parameter selected from the cross validation procedure}
#' @return \item{full.tree}{unpruned tree}
#' @return \item{pruned.tree}{output from pruning of `full.tree`}
#' @return \item{data}{input data}
#' @return \item{details}{summary of model performance}
#' @return \item{subtrees}{sequence of optimally pruned subtrees of `full.tree`}
#' @return \item{in.train}{training samples from splits}
#' @return \item{in.test}{testing samples from splits}
#' @import randomForest
#' @export
#' @examples
#' 
#' # Grow large tree
#' set.seed(1)
#' dat <- generateData()
#' fit <- treeCV(dat, split.var = 1:10, nfolds = 5, lambda = 1,
#'                  risk.control = TRUE, risk.threshold = 2.75)
#' 


treeCV <- function(dat, 
                   split.var, 
                   N0 = 20, 
                   n0 = 5, 
                   efficacy = 'y',
                   risk = 'r',
                   col.trt = "trt",
                   col.prtx = "prtx",
                   lambda = 0,
                   risk.control = FALSE, 
                   risk.threshold = NA, 
                   nfolds = 10, 
                   AIPWE = FALSE, 
                   sort = TRUE, 
                   ctgs = NA, 
                   stabilize.type = c('linear', 'rf'), 
                   stabilize = TRUE, 
                   use.other.nodes = TRUE, 
                   use.bootstrap = FALSE,
                   extremeRandomized = FALSE){
  
  # Grow initial tree
  tre <- grow.ITR(data = dat, 
                  split.var = split.var, 
                  test = NULL, 
                  min.ndsz = N0, 
                  n0 = n0, 
                  efficacy = efficacy, 
                  risk = risk, 
                  col.trt = col.trt,
                  col.prtx = col.prtx,
                  lambda = lambda,
                  risk.control = risk.control, 
                  risk.threshold = risk.threshold, 
                  stabilize = stabilize,
                  stabilize.type = stabilize.type,
                  ctg = ctg, 
                  max.depth = max.depth, 
                  AIPWE = AIPWE, 
                  mtry = length(split.var),
                  use.other.nodes = use.other.nodes, 
                  extremeRandomized = extremeRandomized)
  
  input.tre <- tre$tree
  input.dat <- dat
  
  # Model residuals if requested
  stabilize.type <- match.arg(stabilize.type)
  
  # Shuffle data
  if(!use.bootstrap){
    if(sort) input.dat <- input.dat[sample(1:nrow(input.dat), size = nrow(input.dat)),]
    folds <- cut(seq(1,nrow(input.dat)), breaks = nfolds, labels = FALSE)
  } else{
    if(sort) input.dat <- input.dat[sample(1:nrow(input.dat), size = nrow(input.dat)),]
    folds <- lapply(1:nfolds, function(i) unique(sample(1:nrow(input.dat), nrow(input.dat), replace = TRUE)))
  }
  
  in.train <- in.test <- trees <- list()
  
  # divide samples into training and testing based on n.folds specified
  for(k in 1:nfolds){
    if(!use.bootstrap){ # use traditional CV
      in.train[[k]] <- input.dat[-which(folds==k,arr.ind=TRUE),]
      in.test[[k]]  <- input.dat[which(folds==k,arr.ind=TRUE),]
    } else{
      in.train[[k]] <- input.dat[folds[[k]],]
      in.test[[k]] <- input.dat[-folds[[k]],]
    }
  }
  
  trees <- lapply(1:nfolds, function(n) 
    grow.ITR(data = in.train[[n]], 
             test = in.test[[n]], 
             split.var = sp.var, 
             risk.control = risk.control, 
             risk.threshold = risk.threshold, 
             lambda = lambda,
             efficacy = efficacy,
             risk = risk, 
             col.trt = col.trt,
             col.prtx = col.prtx,
             min.ndsz = N0, 
             n0 = n0, 
             AIPWE = AIPWE,
             ctg = ctgs, 
             mtry = length(split.var),
             stabilize = stabilize, 
             stabilize.type = stabilize.type,
             max.depth = 10, 
             use.other.nodes = use.other.nodes))
  
  out <- lapply(1:length(trees), function(tt){
    if(!is.null(dim(trees[[tt]]$tree) & nrow(trees[[tt]]$tree) != 1)){
      tmp <- prune(trees[[tt]], 0, 
                   test = trees[[tt]]$test, 
                   AIPWE = FALSE, ctgs = ctgs, 
                   risk.control = risk.control, 
                   risk.threshold = risk.threshold, 
                   lambda = lambda)
      out <- as.numeric(tmp$result$V.test)
      names(out) <- tmp$result$size.tmnl
      return(out)
    }
  })
  
  outAlpha <- lapply(1:length(trees), function(tt){
    if(!is.null(dim(trees[[tt]]$tree) & nrow(trees[[tt]]$tree) != 1)){
      tmp <- prune(trees[[tt]], 0, 
                   test = trees[[tt]]$test, 
                   AIPWE = FALSE, ctgs = ctgs, 
                   risk.control = risk.control, 
                   risk.threshold = risk.threshold, 
                   lambda = lambda)
      out <- data.frame(V.test = as.numeric(tmp$result$V.test), 
               alpha = as.numeric(tmp$result$alpha))
      return(out)
    }
  })
  
  best.alpha2 <- do.call(c, lapply(outAlpha, function(tt){
    as.numeric(tt$alpha[which.max(tt$V.test)])
  }))
  best.alpha2 <- mean(best.alpha2[best.alpha2 != 9999])
  
  min.length <- min(sapply(out, length))
  max.length <- max(sapply(out, length))
  
  validSummary <- matrix(sapply(out, function(i){
    out <- rev(i)
    length(out) <- max.length
    return(out)}), ncol = nfolds)
  colnames(validSummary) <- paste0("Tree", 1:nfolds)
  rownames(validSummary) <- paste0("n.tmnl=", 1:nrow(validSummary))
  
  if(risk.control){
    validSummary2 <- matrix(sapply(out, function(i){
      out2 <- rev(i)
      length(out2) <- max.length
      return(out2)}), ncol = nfolds)
    colnames(validSummary2) <- paste0("Tree", 1:nfolds)
    rownames(validSummary2) <- paste0("n.tmnl=", 1:nrow(validSummary2))
  }
  
  m <- apply(validSummary, 1, function(i){
    ifelse(mean(is.na(i)) <= 0.1, mean(i, na.rm = TRUE), NA)
  })
  SD <- apply(validSummary, 1, function(i){
    ifelse(mean(is.na(i)) <= 0.1, sd(i, na.rm = TRUE), NA)
  })
  l <- ncol(validSummary)
  
  if(risk.control){
    m2 <- apply(validSummary2, 1, function(i){
      ifelse(mean(is.na(i)) <= 0.1, mean(i, na.rm = TRUE), NA)
    })
    SD2 <- apply(validSummary2, 1, function(i){
      ifelse(mean(is.na(i)) <= 0.1, sd(i, na.rm = TRUE), NA)
    })
    l <- ncol(validSummary2)
    
  }
  
  full.tre.prune <- prune(tre, 0, ctgs = ctgs,
                          risk.control = risk.control, 
                          risk.threshold = tre$risk.threshold, 
                          lambda = lambda)
  tmp.l <- nrow(full.tre.prune$result)
  final.length <- min(tmp.l, max.length)
  
  result <- data.frame(tail(full.tre.prune$result, n = final.length)[,1:6], 
                       m = rev(head(m, n = final.length)), 
                       SD = rev(head(SD, n = final.length)), 
                       lower = rev(head(m, n = final.length)) - rev(head(SD, n = final.length))/sqrt(l), 
                       upper = rev(head(m, n = final.length)) + rev(head(SD, n = final.length))/sqrt(l))
  if(risk.control) result$m2 <- rev(head(m2, n = final.length))
  
  result <- result[complete.cases(result),]
  
  rm(m, SD, l)
  
  tmp.idx <- which.max(result$m)
  optSize <- result$size.tmnl[tmp.idx]
  best.alpha <- result$alpha[tmp.idx]
  best.subtree <- result$subtree[tmp.idx]
  
  if(optSize == "1"){
    best.tree <- full.tre.prune$subtrees[[length(full.tre.prune$subtrees)]][1,,drop = FALSE]
    best.tree[,6:ncol(best.tree)] <- NA
  } else{
    best.tree <- full.tre.prune$subtrees[[as.numeric(best.subtree)]]
  }

  best.tree2.idx <- which.min(abs(as.numeric(full.tre.prune$result$alpha) - as.numeric(best.alpha2)))
  if(length(full.tre.prune$subtrees) < best.tree2.idx){
    best.tree.alpha <- full.tre.prune$subtrees[[length(full.tre.prune$subtrees)]][1,,drop = FALSE]
    best.tree.alpha[,6:ncol(best.tree.alpha)] <- NA
  } else{
    best.tree.alpha <- full.tre.prune$subtrees[[which.min(abs(as.numeric(full.tre.prune$result$alpha) - as.numeric(best.alpha2)))]]
  }

  setNames(list(best.tree, best.tree.alpha, best.alpha, best.alpha2, input.tre, full.tre.prune$result, 
                input.dat, result, full.tre.prune$subtrees, in.train, in.test), 
           c("best.tree.size", "best.tree.alpha", "best.alpha", "best.alpha2", "full.tree", 
             "pruned.tree", "data", "details", "subtrees", "in.train", "in.test"))
}
