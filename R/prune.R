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
                  train, 
                  outcome = c('time', 'ae'),
                  test = NULL, 
                  AIPWE = FALSE, 
                  n0 = 5, 
                  ctgs = NULL, 
                  haoda.method = FALSE, 
                  haoda.ae.level = NA){
  
  outcome <- match.arg(outcome)
  tre.in <- tre$tree
  train$y <- tre$y 
  if(!is.null(test)) test$y <- tre$y.test 
  
  # Handle null tre case
  if(is.null(dim(tre.in))){
    warning("No Need to Prune Further.")
    return(NA)
  }
  
  if(haoda.method){
    dat.haoda <- data.frame(train[,c("y", "trt", "prtx", "id", "ae")])
    if(!is.null(test)) dat.haoda.test <- data.frame(test[,c("y", "trt", "prtx", "id", "ae")])
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

  while(n.tmnl > 1){
    #internal keeps track of all splits which are not terminal <NA> for score value
    subtrees[[subtree]] <- tre.in
    internal <- tre.in$node[!is.na(tre.in$cut.1)]
    l <- length(internal)
    #r.value is the vector of mean score values across all splits
    r.value <- 
      sapply(1:l, function(xxx){
        #branch keeps track of all splits (terminal or not)
        #branch is a single path which can be followed down a given tree
        nodes.keep <- c(tre.in$node[!tre.in$node %in% de(internal[xxx], tree = tre.in)])
        tmp <- tre.in[tre.in$node %in% nodes.keep , ]
        tmp[tmp$node == internal[xxx], 6:ncol(tmp)] <- NA
        
        if(nrow(tmp) > 1){
          trt.pred <- predict.ITR(tmp, train, ctgs = ctgs)$trt.pred

          if(!haoda.method){
            score <- switch(outcome, 
                            time = s.itrtest(train, trt.pred, n0, AIPWE), 
                            ae = itrtest(train, trt.pred, n0, AIPWE)) 
          } else{
            score <- itrtest(dat.haoda, trt.pred, n0, AIPWE)
            ae.score <- mean(train$trt*trt.pred*train$ae / train$prtx + 
                               (1-train$trt)*(1-trt.pred)*train$ae / (1-train$prtx))
          }
        }else{
          if(!haoda.method){
            score <- switch(outcome, 
                            time = max(s.itrtest(train, rep(0,nrow(train)), 0, AIPWE), 
                                       s.itrtest(train, rep(1,nrow(train)), 0, AIPWE)), 
                            ae = max(itrtest(train, rep(0,nrow(train)), 0, AIPWE), 
                                     itrtest(train, rep(1,nrow(train)), 0, AIPWE)))
          } else{
            score <- max(sapply(0:1, function(xx) 
              itrtest(dat.haoda, rep(xx, nrow(dat.haoda)), -1, AIPWE)))
            ae.score <- max(mean((train$trt==0)*train$ae / train$prtx), 
                            mean((train$trt==1)*train$ae / train$prtx))
          }
        }
        if(!haoda.method){
          return(score / sum(!is.na(tmp)))
        } else{
          return(c(y = score / sum(!is.na(tmp)), r = ae.score))
        }
      })
# browser()
    if(nrow(tre.in) == 3){
      if(!haoda.method){
        alpha <- min(r.value, na.rm = TRUE)
      } else{
        alpha <- min(r.value["y",], na.rm = TRUE)
      }
    } else{
      if(!haoda.method){
        alpha <- min(r.value[-1], na.rm = TRUE)
      } else{
        if(sum(r.value[2,] <= haoda.ae.level) > 0){
          risk.idx <- (r.value["r",] <= haoda.ae.level)
          max.order <- rank(r.value["y",])
          tmp.alpha <- as.matrix(cbind(t(r.value), risk.idx, max.order))[risk.idx,]
          if(is.null(dim(tmp.alpha))){
            alpha <- tmp.alpha["y"]
          } else{
            if(nrow(tmp.alpha) > 1){
              alpha <- data.frame(tmp.alpha)[which.min(as.numeric(tmp.alpha[-1,1]))+1,1]
            } else{
              alpha <- data.frame(tmp.alpha)[which.min(as.numeric(tmp.alpha[,1])),1]
            }
          }
        } else{
          alpha <- r.value["y",-1][which.min(r.value[1,-1])]
        }
        r.value <- as.numeric(r.value["y",])
      }
    }

    nod.rm <- sample(internal[r.value == alpha], 1)
    trt.pred <- predict.ITR(tre.in, train, ctgs)$trt.pred
    
    if(!haoda.method){
      V <- switch(outcome, 
                  time = s.itrtest(train, trt.pred, n0, AIPWE),
                  ae = itrtest(train, trt.pred, n0, AIPWE))
    } else{
      V <- itrtest(dat.haoda, trt.pred, n0, AIPWE)
      V.ae <- itrtest(data.frame(y = dat.haoda$ae, dat.haoda[,2:4]), trt.pred, n0, AIPWE)
      tmp.v.ae <- c(tmp.v.ae, V.ae)
    }
    V.a <- V - a*sum(!is.na(tre.in$score))
    
    if(!is.null(test)){
      # Calculate value for the training set
      trt.pred <- predict.ITR(tre.in, test, ctgs = ctgs)$trt.pred
      if(!haoda.method){
        V.test <- switch(outcome, 
                         time = s.itrtest(dat = test, z=trt.pred, 0, AIPWE), 
                         ae = itrtest(dat = test, z=trt.pred, 0, AIPWE)) 
      } else{
        V.test <- itrtest(dat = dat.haoda.test, z=trt.pred, 0, AIPWE)
        V.ae.test <- itrtest(data.frame(y = dat.haoda.test$ae, dat.haoda.test[,2:4]), trt.pred, n0, AIPWE)
        tmp.v.ae.test <- c(tmp.v.ae.test, V.ae.test)
      }
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
# browser()
  # HANDLE THE NULL TREE WITH THE ROOT NODE ONLY
  if(!haoda.method){
    switch(outcome, ae = {
      if(!is.null(test)){
        result <- rbind(result, cbind(subtree=subtree, node.rm='NA', size.tree=nrow(tre.in), 
                                      size.tmnl=1, alpha=9999, 
                                      V=max(itrtest(train, rep(1,nrow(train)), 5, AIPWE), itrtest(train, rep(0, nrow(train)), 5, AIPWE)), 
                                      V.a=max(itrtest(train, rep(1,nrow(train)), 5, AIPWE), itrtest(train, rep(0, nrow(train)), 5, AIPWE)), 
                                      V.test=max(itrtest(test, rep(1,nrow(test)), 5, AIPWE), itrtest(test, rep(0, nrow(test)), 5, AIPWE)), 
                                      Va.test=max(itrtest(test, rep(1,nrow(test)), 5, AIPWE), itrtest(test, rep(0, nrow(test)), 5, AIPWE))))
      }else{
        result <- rbind(result, cbind(subtree=subtree, node.rm='NA', size.tree=nrow(tre.in), 
                                      size.tmnl=1, alpha=9999, 
                                      V=max(itrtest(train, rep(1,nrow(train)), 5, AIPWE), itrtest(train, rep(0, nrow(train)), 5, AIPWE)), 
                                      V.a=max(itrtest(train, rep(1,nrow(train)), 5, AIPWE), itrtest(train, rep(0, nrow(train)), 5, AIPWE)), 
                                      V.test=NA, Va.test=NA))    
      }
    }, 
    time = {
      if(!is.null(test)){
        result <- rbind(result, cbind(subtree=subtree, node.rm='NA', size.tree=nrow(tre.in), 
                                      size.tmnl=1, alpha=9999, 
                                      V=max(s.itrtest(train, rep(1,nrow(train)), 5, AIPWE), s.itrtest(train, rep(0, nrow(train)), 5, AIPWE)), 
                                      V.a=max(s.itrtest(train, rep(1,nrow(train)), 5, AIPWE), s.itrtest(train, rep(0, nrow(train)), 5, AIPWE)), 
                                      V.test=max(s.itrtest(test, rep(1,nrow(test)), 5, AIPWE), s.itrtest(test, rep(0, nrow(test)), 5, AIPWE)), 
                                      Va.test=max(s.itrtest(test, rep(1,nrow(test)), 5, AIPWE), s.itrtest(test, rep(0, nrow(test)), 5, AIPWE))))
      }else{
        result <- rbind(result, cbind(subtree=subtree, node.rm='NA', size.tree=nrow(tre.in), 
                                      size.tmnl=1, alpha=9999, 
                                      V=max(s.itrtest(train, rep(1,nrow(train)), 5, AIPWE), s.itrtest(train, rep(0, nrow(train)), 5, AIPWE)), 
                                      V.a=max(s.itrtest(train, rep(1,nrow(train)), 5, AIPWE), s.itrtest(train, rep(0, nrow(train)), 5, AIPWE)), 
                                      V.test=NA, Va.test=NA))    
      }
    })
  } else{
    if(!is.null(test)){
      result <- rbind(result, cbind(subtree=subtree, node.rm='NA', size.tree=nrow(tre.in), 
                                    size.tmnl=1, alpha=9999, 
                                    V=max(itrtest(dat.haoda, rep(1,nrow(dat.haoda)), 5, AIPWE), itrtest(dat.haoda, rep(0, nrow(dat.haoda)), 5, AIPWE)), 
                                    V.a=max(itrtest(dat.haoda, rep(1,nrow(dat.haoda)), 5, AIPWE), itrtest(dat.haoda, rep(0, nrow(dat.haoda)), 5, AIPWE)), 
                                    V.test=max(itrtest(dat.haoda.test, rep(1,nrow(dat.haoda.test)), 5, AIPWE), itrtest(dat.haoda.test, rep(0, nrow(dat.haoda.test)), 5, AIPWE)), 
                                    Va.test=max(itrtest(dat.haoda.test, rep(1,nrow(dat.haoda.test)), 5, AIPWE), itrtest(dat.haoda.test, rep(0, nrow(dat.haoda.test)), 5, AIPWE))))
    }else{
      result <- rbind(result, cbind(subtree=subtree, node.rm='NA', size.tree=nrow(tre.in), 
                                    size.tmnl=1, alpha=9999, 
                                    V=max(itrtest(dat.haoda, rep(1,nrow(dat.haoda)), 5, AIPWE), itrtest(dat.haoda, rep(0, nrow(dat.haoda)), 5, AIPWE)), 
                                    V.a=max(itrtest(dat.haoda, rep(1,nrow(dat.haoda)), 5, AIPWE), itrtest(dat.haoda, rep(0, nrow(dat.haoda)), 5, AIPWE)), 
                                    V.test=NA, Va.test=NA))    
    }
    
    if(mean(train$time * train$trt / train$prtx) >  
       mean(train$time * (1-train$trt) / train$prtx)){
      tmp.trt <- rep(1,nrow(dat.haoda))
      if(!is.null(test)) tmp.test.trt <- rep(1,nrow(dat.haoda.test))
    } else{
      tmp.trt <- rep(0,nrow(dat.haoda))
      if(!is.null(test)) tmp.test.trt <- rep(0,nrow(dat.haoda.test))
    }
    tmp.v.ae <- c(tmp.v.ae, mean(dat.haoda$ae * (tmp.trt==dat.haoda$trt) / dat.haoda$prtx))
    if(!is.null(test)) tmp.v.ae.test <- 
      c(tmp.v.ae.test, mean(dat.haoda.test$ae * (tmp.test.trt==dat.haoda.test$trt) / dat.haoda.test$prtx))
  }
  result <- as.data.frame(result)
  result <- result[!duplicated(result),]
  if(haoda.method){
    out <- list(result = result, subtrees = subtrees, v.ae = tmp.v.ae)
    out$v.ae.test <- tmp.v.ae.test
  } else{
    out <- list(result = result, subtrees = subtrees)  
  }
  return(out)
}