#' Performs k-fold cross validation for selection optimal tuning parameter lambda. 
#' 
#' @param tre original large tree from grow.ITR() 
#' @param dat data frame used to grow tre
#' @param N0 sets the minimal number of observations allowed in terminal nodes. Defaults to 20. 
#' @param n0 sets the minimum number of observations from each treatment group to be in each child node.  Defaults to 5. 
#' @param sp.var specifies the columns of splitting variables in the input data to grow the original tree. 
#' @param nfolds number of folds in the cross validation. Defaults to 10.
#' @param param vector of numeric values to be considered as the tuning parameter. Defaults to seq(0, 0.15, 0.1) but should be modified for each specific problem. 
#' @param AIPWE logical indicating if the augmented estimator should be used. Defaults to FALSE.
#' @param sort logical indicating if data should be sorted before folds are created. Defaults to FALSE. 
#' @param stabilize logical. Should stabilization be used?
#' @param stabilize.type Type of stabilization. 'rf' for random forests (default) or 'linear' for linear model.  
#' @param ctgs columns of categorical variables. 
#' @param haoda.method logical if the conditional tree structure is being used
#' @param haoda.ae.level value used if haoda.method is TRUE.
#' @return A summary of the cross validation including optimal penalty parameter and the optimal model. 
#' @return \item{best.tree}{optimal ITR tree model}
#' @return \item{best.lambda}{optimal lambda}
#' @return \item{full.tree}{original input tree}
#' @return \item{pruned.tree}{pruning output given the optimal lambda}
#' @return \item{data}{input data}
#' @return \item{details}{summary of model performance for each lambda under consideration}
#' @import randomForest
#' @export
#' @examples
#' 
#' # Grow large tree
#' set.seed(1)
#' dat <- gdataM(n=1000, depth = 2, beta1=3, beta2=1)
#' tre <- grow.ITR(dat, split.var = 1:4)
#' 
#' # This tre should have 3 terminal nodes (correct size), but has 4 as grown.
#' 
#' cv.prune <- treeCV(tre, dat, nfolds = 5, param = seq(0, 0.15, 0.01), sp.var = 1:4)
#' 
#' # The best tree returned has the correct number of nodes
#' cv.prune$best.tree
#' #node size n.1 n.0 trt.effect var vname cut.1 cut.2  score
#' #1    0 1000 479 521  1.2936724   1    X1     r  0.28 0.8768
#' #2   01  277  86 191 -0.9439306  NA  <NA>  <NA>  <NA>     NA
#' #3   02  723 393 330  2.0819738   3    X3     r   0.1 0.9426
#' #5  021   78   9  69 -1.3600766  NA  <NA>  <NA>  <NA>     NA
#' #4  022  645 384 261  2.4226122  NA  <NA>  <NA>  <NA>     NA
#' 
#' 


treeCV <- function(tre, 
                   sp.var, 
                   dat, 
                   outcome = c("time", "ae"),
                   nfolds = 10, 
                   param = NULL, 
                   AIPWE = FALSE, 
                   N0=20, 
                   n0=5, 
                   sort = TRUE, 
                   ctgs = NA, 
                   standardize = FALSE,
                   stabilize.type = c('linear', 'rf'), 
                   stabilize = TRUE, 
                   scale.ae = FALSE,
                   haoda.method = FALSE, 
                   haoda.ae.level = NA, 
                   reverse.ae.scale = FALSE, 
                   use.other.nodes = TRUE, 
                   use.bootstrap = FALSE){

  input.tre <- tre$tree
  input.dat <- dat
  input.dat$y <- tre$y

  # Model residuals if requested
  outcome <- match.arg(outcome)
  stabilize.type <- match.arg(stabilize.type)
  
  if(is.null(param)){
    param <- seq(0, max(input.tre$score, na.rm = TRUE), length.out = 1000)
  }
  
  # Shuffle data
  if(!use.bootstrap){
    if(sort) input.dat <- input.dat[sample(1:nrow(input.dat), size = nrow(input.dat)),]
    folds <- cut(seq(1,nrow(input.dat)), breaks = nfolds, labels = FALSE)
  } else{
    if(sort) input.dat <- input.dat[sample(1:nrow(input.dat), size = nrow(input.dat)),]
    folds <- lapply(1:nfolds, function(i) unique(sample(1:nrow(input.dat), nrow(input.dat), replace = TRUE)))
  }

  in.train <- in.test <- trees <- list()

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
      grow.ITR(data = in.train[[n]], test = in.test[[n]], 
               split.var = sp.var, standardize = standardize, outcome = outcome,
               min.ndsz = N0, n0 = n0, AIPWE = AIPWE, ctg = ctgs, 
               stabilize = stabilize, stabilize.type = stabilize.type, 
               max.depth = 5, haoda.method = haoda.method, 
               haoda.ae.level = haoda.ae.level, 
               scale.ae = scale.ae, 
               reverse.ae.scale = reverse.ae.scale, 
               use.other.nodes = use.other.nodes))

  out <- lapply(1:length(trees), function(tt){
    if(!is.null(dim(trees[[tt]]$tree) & nrow(trees[[tt]]$tree) != 1)){
      tmp <- prune(trees[[tt]], 0, train = in.train[[tt]], test = in.test[[tt]], 
                   AIPWE = FALSE, ctgs = ctgs, outcome = outcome, 
                   haoda.method = haoda.method, 
                   haoda.ae.level = haoda.ae.level)
      
      if(haoda.method){
        out <- list(as.numeric(tmp$result$V.test), 
                    tmp$v.ae)
        names(out[[1]]) <- tmp$result$size.tmnl
        names(out[[2]]) <- tmp$result$size.tmnl
      } else{
        out <- as.numeric(tmp$result$V.test)
        names(out) <- tmp$result$size.tmnl
      }
      return(out)
    }
  })

  if(haoda.method) out2 <- sapply(out, "[", 2)
  out <- sapply(out, "[", 1)
  min.length <- min(sapply(out, length))
  max.length <- max(sapply(out, length))
  
  validSummary <- matrix(sapply(out, function(i){
    out <- rev(i)
    length(out) <- max.length
    return(out)}), ncol = nfolds)
  colnames(validSummary) <- paste0("Tree", 1:nfolds)
  rownames(validSummary) <- paste0("n.tmnl=", 1:nrow(validSummary))
  
  if(haoda.method){
    validSummary2 <- matrix(sapply(out2, function(i){
      out2 <- rev(i)
      length(out2) <- max.length
      return(out2)}), ncol = nfolds)
    colnames(validSummary2) <- paste0("Tree", 1:nfolds)
    rownames(validSummary2) <- paste0("n.tmnl=", 1:nrow(validSummary2))
  }
  
  m <- apply(validSummary, 1, function(i){
    ifelse(mean(is.na(i)) <= 0.25, mean(i, na.rm = TRUE), NA)
  })
  SD <- apply(validSummary, 1, function(i){
    ifelse(mean(is.na(i)) <= 0.25, sd(i, na.rm = TRUE), NA)
  })
  l <- ncol(validSummary)

  if(haoda.method){
    m2 <- apply(validSummary2, 1, function(i){
      ifelse(mean(is.na(i)) <= 0.25, mean(i, na.rm = TRUE), NA)
    })
    SD2 <- apply(validSummary2, 1, function(i){
      ifelse(mean(is.na(i)) <= 0.25, sd(i, na.rm = TRUE), NA)
    })
    l <- ncol(validSummary2)
    
  }

  full.tre.prune <- prune(tre, 0, dat, outcome, ctgs = ctgs,
                          haoda.method = haoda.method, 
                          haoda.ae.level = haoda.ae.level)
  tmp.l <- nrow(full.tre.prune$result)
  final.length <- min(tmp.l, max.length)

  result <- data.frame(tail(full.tre.prune$result, n = final.length)[,1:6], 
                       m = rev(head(m, n = final.length)), 
                       SD = rev(head(SD, n = final.length)), 
                       lower = rev(head(m, n = final.length)) - rev(head(SD, n = final.length))/sqrt(l), 
                       upper = rev(head(m, n = final.length)) + rev(head(SD, n = final.length))/sqrt(l))
  if(haoda.method) result$m2 <- rev(head(m2, n = final.length))
  result <- result[complete.cases(result),]

  rm(m, SD, l)
  
  if(haoda.method){
    ae.idx <- result$m2 <= haoda.ae.level
    if(sum(ae.idx) == 0) ae.idx <- !is.na(as.numeric(result$node.rm)) 
    tmp.idx <- which(result$m[ae.idx] == max(result$m[ae.idx]))
  } else{
    tmp.idx <- which.max(result$m)
  }
  optSize <- result$size.tmnl[tmp.idx]
  best.alpha <- result$alpha[tmp.idx]
  best.subtree <- result$subtree[tmp.idx]
  
  if(optSize == "1"){
    best.tree <- full.tre.prune$subtrees[[length(full.tre.prune$subtrees)]][1,,drop = FALSE]
    best.tree[,6:ncol(best.tree)] <- NA
  } else{
    best.tree <- full.tre.prune$subtrees[[as.numeric(best.subtree)]]
  }

  setNames(list(best.tree, best.alpha, input.tre, full.tre.prune$result, 
                input.dat, result, full.tre.prune$subtrees), 
           c("best.tree", "best.lambda", "full.tree", 
             "pruned.tree", "data", "details", "subtrees"))
}