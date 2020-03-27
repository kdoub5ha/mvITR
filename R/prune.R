#' @title Prunes a tree for a given penalty value
#' 
#' @description The `prune` function allows the user to specify a value of the penalty for a given tree. 
#' This function uses the "weakest link" criteria in order to evaluate the order in which branches are pruned
#' and gives the penalized value along with the unpenalized value. If testing data are provided for validation, 
#' then the penalized and unpenalized values from the testing data run down the tree structure are also provided.  
#' 
#' @param tre sets the tree to be pruned 
#' @param a sets the value of the splitting penalty
#' @param train the training data used to create the tree
#' @param test the testing data to be used.  Defaults to NULL.
#' @param AIPWE indicator for AIPWE estimation.
#' @param n0 minimum number of observations allowed in a treatment group. Defaults to 5. 
#' @param ctgs columns of categorical variables.
#' @return summary of pruned branches and the associated value of the tree after pruning. 
#' @return \item{result}{contains columns: `node.rm` which is the weakest link at each
#' iteration of the pruning algorithm; `size.tree` and `n.tmnl` give number of total nodes and
#' terminal nodes in each subtree; `alpha` is the weakest link scoring criteria; `V` and `V.test`
#' are the overall value of the tree for the training and tesing samples; `V.a` and `Va.test`
#' give the penalized value for the training and testing samples.}
#' @export
#' 


prune <- function(tre, 
                  a, 
                  test = NULL, 
                  AIPWE = FALSE, 
                  n0 = 5, 
                  ctgs = NULL, 
                  haoda.method = FALSE, 
                  haoda.ae.level = NA,
                  lambda = lambda){
  
  tre.in <- tre$tree
  train <- tre$data 
  if(!is.null(test) & is.null(test$y)) test$y <- tre$y.test 
  if(!is.null(test) & is.null(test$r)) test$r <- tre$test$r
  
  # Handle null tre case
  if(is.null(dim(tre.in))){
    warning("No Need to Prune Further.")
    return(NA)
  }
  
  # If there are at least three terminal nodes we will determine pruning 
  tmnl.idx <- is.na(tre.in$var)
  n.tmnl <- sum(tmnl.idx)
  subtrees <- vector("list")
  subtree <- 1
  result <- data.frame()
  if(haoda.method){
    tmp.v.ae <- vector("numeric")
    tmp.v.ae.test <- vector("numeric") 
  }
# browser()
  while(n.tmnl > 1){
    #internal keeps track of all splits which are not terminal <NA> for score value
    subtrees[[subtree]] <- tre.in
    internal <- tre.in$node[!is.na(tre.in$cut.1)]
    l <- length(internal)
    preds.tre.in <- predict.ITR(tre.in, train, split.vars)$trt.pred
    L.tre.in <- estITR(list(y = train$y, trt = train$trt, ae = train$r,
                            maxRisk = haoda.ae.level,
                            prtx = train$prtx, status = train$status, 
                            lambda = lambda, KM.cens = train$KM.cens, 
                            n0 = 0, z = preds.tre.in))
    
    #r.value is the vector of mean score values across all splits
    r.value <- 
      sapply(1:l, function(xxx){
        #branch keeps track of all splits (terminal or not)
        #branch is a single path which can be followed down a given tree
        nodes.keep <- c(tre.in$node[!tre.in$node %in% de(internal[xxx], tree = tre.in)])
        tmp <- tre.in[tre.in$node %in% nodes.keep , , drop=F]
        tmp[tmp$node == internal[xxx], 6:ncol(tmp)] <- NA
        
        if(nrow(tmp) > 1){
          trt.pred <- predict.ITR(tmp, train, tre$split.var, ctgs = ctgs)$trt.pred
          # score <- estITR(list(y = train$y, trt = train$trt, ae = train$r,
          #                      prtx = train$prtx, status = train$status, 
          #                      KM.cens = train$KM.cens, n0 = 5, z = trt.pred,
          #                      lambda = lambda, maxRisk = haoda.ae.level))
          # 
          ae.score <- mean(train$r * (train$trt == trt.pred) / train$prtx)
          y.score <- mean(train$y * (train$trt == trt.pred) / train$prtx)
          score <- estITR(list(y = train$y, trt = train$trt, ae = train$r,
                               maxRisk = haoda.ae.level,
                               prtx = train$prtx, status = train$status, 
                               lambda = lambda, KM.cens = train$KM.cens, 
                               n0 = 0, z = trt.pred))
        }else{
          scores <- sapply(0:1, function(iii)
            estITR(list(y = train$y, trt = train$trt, ae = train$r, maxRisk = haoda.ae.level,
                        prtx = train$prtx, status = train$status, lambda = lambda,
                        KM.cens = train$KM.cens, n0 = 0, z = rep(iii, nrow(train)))))
          idx.scores <- which.max(scores)
          score <- scores[idx.scores]
          ae.score <- mean(train$r * (train$trt == rep(idx.scores-1, nrow(train))) / train$prtx)
          y.score <- mean(train$y * (train$trt == rep(idx.scores-1, nrow(train))) / train$prtx)
        }
        
        if(!haoda.method){
          return(score / sum(!is.na(tmp$var)))
        } else{
          # return(c(y = score / sum(!is.na(tmp)), r = ae.score))
          return(c(y = score,
                   y.diff = L.tre.in - score,
                   y.score = y.score,
                   r.score = ae.score))
        }
      })

    if(nrow(tre.in) > 1){
      if(!haoda.method){
        alpha <- max(r.value, na.rm = TRUE)
      } else{
        alpha <- max(r.value["y",], na.rm = TRUE)
        r.value <- as.numeric(r.value["y",])
      }
    } else{
      if(!haoda.method){
        alpha <- max(r.value[-1], na.rm = TRUE)
      } else{
        # alpha <- max(r.value["y",], na.rm = TRUE)
        
        if(sum(r.value["r",] <= haoda.ae.level) > 0){
          risk.idx <- (r.value["r",] <= haoda.ae.level)
          max.order.y <- rank(r.value["y",])
          max.order.r <- rank(r.value["r",])
          tmp.alpha <- as.matrix(cbind(t(r.value), risk.idx, 
                                       max.order.y, max.order.r))[risk.idx,,drop=FALSE]
          if(is.null(dim(tmp.alpha))){
            alpha <- tmp.alpha["y"]
          } else{
            if(nrow(tmp.alpha) > 1){
              alpha <- data.frame(tmp.alpha)[which.max(as.numeric(tmp.alpha[,"y"])),"y"]
            } else{
              alpha <- data.frame(tmp.alpha)[which.max(as.numeric(tmp.alpha[,"y"])),"y"]
            }
          }
        } else{
          alpha <- r.value["y",-1][which.max(r.value["y",-1])]
        }
        r.value <- as.numeric(r.value["y",])
      }
    }
    
    nod.rm <- sample(internal[r.value == alpha], 1)
    trt.pred <- preds.tre.in
    
    V <- L.tre.in
    V.a <- V - a*sum(!is.na(tre.in$score))
    
    if(!is.null(test)){
      # Calculate value for the training set
      trt.pred <- predict.ITR(tre.in, test, tre$split.var, ctgs = ctgs)$trt.pred
      if(is.null(test$status)) test$status <- rep(1, nrow(test))
      if(is.null(test$KM.cens)) test$KM.cens <- rep(1, nrow(test))
      V.test <- estITR(list(y = .subset2(test, 'y'), 
                            trt = .subset2(test, 'trt'), 
                            ae = .subset2(test, 'r'), 
                            maxRisk = haoda.ae.level,
                            prtx = .subset2(test, 'prtx'),
                            status = .subset2(test, 'status'),
                            lambda = lambda,
                            KM.cens = .subset2(test, 'KM.cens'),
                            n0 = 0, z = trt.pred))
      # V.ae.test <- estITR(list(y = test$r, trt = test$trt, 
      #                          prtx = test$prtx, status = rep(1, nrow(test)), 
      #                          KM.cens = rep(1, nrow(test)), n0 = 0, z = trt.pred))
      # tmp.v.ae.test <- c(tmp.v.ae.test, V.ae.test)
      Va.test <- V.test - a*sum(!is.na(tre.in$score))
    }
    
    # Calculate value for testing data
    if(is.null(test)){
      result <- rbind(result, 
                      data.frame(subtree = subtree, node.rm = nod.rm, size.tree = nrow(tre.in),
                                 size.tmnl = nrow(tre.in)-l, alpha = alpha, V = V, V.a = V.a, 
                                 V.test = NA, Va.test = NA))
    }else{
      result <- rbind(result, 
                      data.frame(subtree = subtree, node.rm = nod.rm, size.tree = nrow(tre.in), 
                                 size.tmnl = nrow(tre.in)-l, alpha = alpha, V = V, V.a = V.a, 
                                 V.test = V.test, Va.test = Va.test))
    }
    
    if(length(nod.rm) > 1){
      for(k in 1:length(nod.rm)){
        tre.in <- tre.in[!tre.in$node %in% de(nod.rm[k], tre.in),]
        if(is.null(test)){
          o <- match(nod.rm[k], tre.in$node)
          if(!is.na(o)){
            tre.in[match(nod.rm[k], tre.in$node), c("var", "vname", "cut.1", "cut.2", "score")] <- NA
          }
        }else{
          o <- match(nod.rm[k], tre.in$node)
          if(!is.na(o)){
            tre.in[match(nod.rm[k], tre.in$node), c("var", "vname", "cut.1", "cut.2", "score", "score.test")] <- NA
          }
        }
        n.tmnl <- sum(is.na(tre.in$cut.1))
        subtree <- subtree + 1  
      }
    } else{
      tre.in <- tre.in[!tre.in$node %in% de(nod.rm,tre.in),]
      tre.in[tre.in$node == nod.rm, 6:(ncol(tre.in) - !is.null(test))] <- NA
      n.tmnl <- sum(is.na(tre.in$var))
      subtree <- subtree + 1
    }
  }
  
  # HANDLE THE NULL TREE WITH THE ROOT NODE ONLY
  tmp.train.v <- max(sapply(0:1, function(iii) 
    estITR(list(y = train$y, trt = train$trt, ae = train$r, maxRisk = haoda.ae.level,
                prtx = train$prtx, status = train$status, lambda = lambda,
                KM.cens = train$KM.cens, n0 = 0, z = rep(iii, nrow(train))))))
  if(!is.null(test)){
    tmp.test.v <- max(sapply(0:1, function(iii) 
      estITR(list(y = test$y, trt = test$trt, ae = test$r, maxRisk = haoda.ae.level,
                  prtx = test$prtx, status = test$status, lambda = lambda,
                  KM.cens = test$KM.cens, n0 = 0, z = rep(iii, nrow(test))))))
    result <- rbind(result, cbind(subtree=subtree, node.rm='NA',
                                  size.tree=nrow(tre.in), 
                                  size.tmnl=1, alpha=9999, 
                                  V = tmp.train.v, 
                                  V.a = tmp.train.v, 
                                  V.test = tmp.test.v, 
                                  Va.test = tmp.test.v))
  } else{
    result <- rbind(result, cbind(subtree=subtree, node.rm='NA', 
                                  size.tree=nrow(tre.in), 
                                  size.tmnl=1, alpha=9999, 
                                  V = tmp.train.v, 
                                  V.a = tmp.train.v, 
                                  V.test=NA, Va.test=NA))    
  }
  
  if(mean(train$y * train$trt / train$prtx) >  
     mean(train$y * (1-train$trt) / train$prtx)){
    tmp.trt <- rep(1,nrow(train))
    if(!is.null(test)) tmp.test.trt <- rep(1,nrow(test))
  } else{
    tmp.trt <- rep(0,nrow(train))
    if(!is.null(test)) tmp.test.trt <- rep(0,nrow(test))
  }
  # tmp.v.ae <- c(tmp.v.ae, mean(train$ae * (tmp.trt==train$trt) / train$prtx))
  # if(!is.null(test)) tmp.v.ae.test <- 
  #   c(tmp.v.ae.test, mean(test$ae * (tmp.test.trt==test$trt) / test$prtx))
  
  result <- as.data.frame(result)
  result <- result[!duplicated(result),]
  if(haoda.method){
    out <- list(result = result, subtrees = subtrees)#, v.ae = tmp.v.ae)
    # out$v.ae.test <- tmp.v.ae.test
  } else{
    out <- list(result = result, subtrees = subtrees)  
  }
  return(out)
}