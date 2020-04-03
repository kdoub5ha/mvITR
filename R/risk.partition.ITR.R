#' @title Generates partition summary based on itr value. 
#' 
#' @description Determines optimal split for a given node. 
#' 
#' @param dat data.frame. Data used to identify split.  
#' @param split.var numeric vector. Columns of spliting variables.
#' @param risk.threshold numeric. Desired level of risk control. 
#' @param lambda numeric. Penalty parameter for risk scores. Defaults to 0, i.e. no constraint.
#' @param test data.frame of testing observations. Should be formatted the same as 'data'.
#' @param min.ndsz numeric specifying minimum number of observations required to call a node terminal. Defaults to 20.
#' @param n0 numeric specifying minimum number of treatment/control observations needed in a split to declare a node terminal. Defaults to 5. 
#' @param max.depth numeric specifying maximum depth of the tree. Defaults to 15 levels. 
#' @param mtry numeric specifying the number of randomly selected splitting variables to be included. Defaults to number of splitting variables.
#' @param dat.rest dataframe. Data outside current splitting node.
#' @param max.score numeric. Current score for the tree. 
#' @param name char. Name of internal node, used for ordering splits.
#' @param use.other.nodes logical. Should global estimator of objective function be used. Defaults to TRUE. 
#' @param ctg numeric vector corresponding to the categorical input columns.  Defaults to NULL.  Not available yet. 
#' @param AIPWE logical. Should AIPWE (TRUE) or IPWE (FALSE) be used. Not available yet. 
#' @param extremeRandomized logical. Experimental for randomly selecting cutpoints in a random forest model. Defaults to #' @return summary of the best split for a given data frame. 
#' @export

risk.partition.ITR <- function(dat, 
                               split.var,
                               test = NULL, 
                               risk.threshold = NA, 
                               min.ndsz = 20, 
                               n0 = 5, 
                               lambda = 0,
                               name = "0", 
                               ctg = ctg,
                               max.depth = 15, 
                               mtry = length(split.var), 
                               dat.rest = NULL, 
                               max.score = NULL, 
                               AIPWE = AIPWE,
                               use.other.nodes = TRUE, 
                               extremeRandomized = FALSE)
{   
  # inialize various variable
  call <- match.call()
  out <- match.call(expand = F)
  out$info <- NULL
  out$name.l <- NULL
  out$name.r <- NULL
  out$left <- NULL
  out$right <- NULL
  out$... <- NULL
  # label the binary tree by 1 (left) and 2 (right).
  name.l <- paste(name, 1, sep="")
  name.r <- paste(name, 2, sep="")
  # sample size
  n <- nrow(dat)
  # check whether testing data is provided
  if (!is.null(test)) {
    n.test <- nrow(test)
    score.test <- NA
  }
  # prepare for the first cut these variable were used to store final cut information
  var <- NA
  vname <- NA
  cut <- NA
  
  if(name=="0") {
    dat.comb <- dat[ ,c('y', 'trt', 'prtx', 'r', 'id')]
  } else{
    dat.comb <- rbind(dat.rest[, c('y', 'trt', 'prtx', 'r', 'id')], 
                      dat[,c('y', 'trt', 'prtx', 'r', 'id')])
  }
  
  # extract value from data
  trt <- .subset2(dat, 'trt')
  y <- .subset2(dat, 'y')
  vnames <- colnames(dat)
  # COMPUTE THE TREATMENT EFFECT IN CURRENT NODE
  trt.effect <- NA
  n.0 = length(y[trt==0])
  n.1 = n - n.0
  if (min(n.1, n.0) >0) {
    trt.effect <- mean(y[trt==1]) - mean(y[trt==0])
  }
  # CONTROL THE MAX TREE DEPTH
  # name is the current tree label.
  # only when currently depth < max.depth and n > min terminal node proceed.
  depth <- nchar(as.character(name)) 
  
  if(depth <= max.depth && n >= min.ndsz) {
    m.try <- ifelse(is.null(mtry), length(split.var), mtry)
    
    # Fulfill splitting conditions
    splitVars <- sample(split.var, size = m.try, replace = FALSE)
    # splitVars <- sample(split.var, size = m.try, replace = FALSE)
    names(splitVars) <- colnames(dat)[splitVars]
    # apply search algorithm to each potential split
    outSplit <- lapply(splitVars, function(i){
      x <- .subset2(dat, i)
      v.name <- vnames[i]
      temp <- sort(unique(x))
      is.ctg <- is.element(i, ctg)
      
      if(length(temp) > 1) {
        # handle categorial variable first, otherwise take out the final value as no cut after it.
        if(is.ctg){
          zcutCat <- power.set(temp)
          zcutCat <- do.call(cbind, lapply(zcutCat, 
                                           function(zzz) as.numeric(x %in% zzz)))
          zcut <- 1:10
        } else{
          zcut <- temp[-length(temp)]
          zcutCat <- matrix(1)
        }
        
        if(extremeRandomized){
          zcut <- sort(sample(zcut, min(length(zcut), 20)))
        }
        
        if(name == "0"){
          if(!is.ctg){
            x.tmp <- x
          } else{
            x.tmp <- as.numeric(factor(x))
          }
          
          # Set up data for splitting
          datMatrix <- list(y = .subset2(dat.comb, 'y'), 
                            ae = .subset2(dat.comb, 'r'), 
                            x = x.tmp, 
                            prtx = .subset2(dat.comb, 'prtx'), 
                            trt = .subset2(dat.comb, 'trt'), 
                            trtNew = rep(-1000, length(x)),
                            inNode = rep(1, length(x)))
        } else{
          if(!is.ctg){
            x.tmp <- c(dat.rest[,i], x)
          } else{
            x.tmp <- as.numeric(factor(c(dat.rest[,i], x)))
          }
          
          # Set up data for splitting
          datMatrix <- list(y = .subset2(dat.comb, 'y'), 
                            ae = .subset2(dat.comb, 'r'), 
                            x = x.tmp, 
                            prtx = .subset2(dat.comb, 'prtx'), 
                            trt = .subset2(dat.comb, 'trt'),
                            trtNew = c(dat.rest$trt.new, 
                                       rep(-1000, length(x))),
                            inNode = c(rep(0, length(dat.rest$y)), 
                                       rep(1, length(x))))
        }
        
        # set up split parameters
        splitParams <- list(isCtg = is.ctg, 
                            nodeSize = min.ndsz, 
                            trtSize = n0, 
                            maxRisk = risk.threshold, 
                            useOtherNodes = as.numeric(use.other.nodes), 
                            maxScore = max.score,
                            lambda = lambda)
        
        tmp <- splitConditional(zcut, zcutCat, datMatrix, splitParams)
        return(tmp)
      } else{
        return(list(output = NA, direction = NA, zcut = NA))
      }
    })
    
    # Extract best split info
    out.idx <- which.max(lapply(outSplit, function(xx) max(xx$output)))
    if(length(out.idx) != 0){
      x <- .subset2(dat, names(out.idx))
      if(extremeRandomized){
        temp <- sort(outSplit[[out.idx]][["zcut"]])
      } else{
        temp <- sort(unique(x))
      }
      tmp.idx <- which.max(outSplit[[out.idx]]$output)
      if(outSplit[[out.idx]]$direction[tmp.idx] %in% c("l", "r")){
        vname <- names(out.idx)
        var <- which(colnames(dat) == vname)
        max.score <- max(outSplit[[out.idx]]$output)
        if(!is.null(ctg) & names(out.idx) %in% colnames(dat)[ctg]){
          best.cut <- paste(unlist(power.set(temp)[tmp.idx]), collapse=",")
        } else{
          if(extremeRandomized){
            best.cut <- temp[tmp.idx]
          } else{
            best.cut <- temp[-length(temp)][tmp.idx]
          }
        }
        
        cut <- cbind(outSplit[[out.idx]]$direction[tmp.idx], 
                     as.character(best.cut))
      }
    }
  }
  # when testing data is provided, assess new treatment assignment 
  # using testing sample and the rule calculated from training sample
  # var is the covariates calcualted before where spliting adopts. 
  # best.cut is the cutting point.
  if (!is.null(test)) {
    n.test <- nrow(test)
    score.test <- NA
    if (!is.na(var)) {
      if (is.element(var,ctg)) {
        if(cut[1] == "l"){
          grp.test <- sign(is.element(test[,var], unlist(strsplit(best.cut, ","))))
        } else{
          grp.test <- sign(!is.element(test[,var], unlist(strsplit(best.cut, ","))))
        }
      } else  {
        if(cut[1]=="l"){
          grp.test <- sign(test[,var] <= best.cut)
        } else{
          grp.test <- sign(test[,var] > best.cut)
        }
      }
      
      score.test <- estITR(list(y = .subset2(test, 'y'), 
                                prtx = .subset2(test, 'prtx'), 
                                ae = .subset2(test, 'r'),
                                trt = .subset2(test, 'trt'), 
                                KM.cens = .subset2(test, 'KM.cens'), 
                                status = .subset2(test, 'status'), 
                                n0 = n0, z = grp.test, 
                                lambda = lambda, 
                                maxRisk = risk.threshold))
      
      if (!is.na(score.test)){
        out$name.l <- name.l
        out$name.r <- name.r
        if(cut[1]=="l"){
          out$left.test <- test[grp.test==1,  ]
          out$right.test <- test[grp.test==0,  ]
        } else{
          out$left.test <- test[grp.test==0,  ]
          out$right.test <- test[grp.test==1,  ]
        }
        if (is.element(var,ctg)) {
          if(cut[1] == "l"){
            out$left  <- dat[is.element(dat[,var], unlist(strsplit(best.cut, ","))),]
            out$right <- dat[!is.element(dat[,var], unlist(strsplit(best.cut, ","))), ]
          } else{
            out$left  <- dat[!is.element(dat[,var], unlist(strsplit(best.cut, ","))),]
            out$right <- dat[is.element(dat[,var], unlist(strsplit(best.cut, ","))), ]
          }
        } else {
          if(cut[1]=='l'){
            out$left  <- cbind(dat[dat[,var]<= best.cut,],new.trt=rep(1,n=sum(dat[,var]<= best.cut)))
            out$right <- cbind(dat[dat[,var]> best.cut, ],new.trt=rep(0,n=sum(dat[,var]> best.cut)))
          }else{
            out$left  <- cbind(dat[dat[,var]<= best.cut,],new.trt=rep(0,n=sum(dat[,var]<= best.cut)))
            out$right <- cbind(dat[dat[,var]> best.cut, ],new.trt=rep(1,n=sum(dat[,var]> best.cut)))
          }  
        }
      } else {
        var <- NA
        vname <- NA
        cut <- NA
        max.score <- NA
      }
      # output results from both testing and training data.
      if(!is.na(var)){  
        out$info <- data.frame(node=name, size = n, n.1=n.1, n.0=n.0, 
                               trt.effect=trt.effect,var = var, 
                               vname=vname, cut.1 = cut[1], cut.2 = cut[2], 
                               score=ifelse(max.score==-1e20, NA, max.score), 
                               score.test=score.test, size.test=n.test)
      } else{
        out$info <- data.frame(node=name, size = n, n.1=n.1, n.0=n.0, 
                               trt.effect=trt.effect, var = NA, 
                               vname=NA, cut.1 = NA, cut.2 = NA, 
                               score=NA,score.test=NA, size.test=n.test)
      }
    } else {
      out$info <- data.frame(node=name, size = n, n.1=n.1, n.0=n.0, trt.effect=trt.effect, var = NA, 
                             vname=NA, cut.1=NA, cut.2=NA, score=NA,score.test=NA, size.test=n.test)
    }
  } else {
    # if no testing data output results from training data only.
    if (!is.na(var)) {
      out$name.l <- name.l
      out$name.r <- name.r
      if (is.element(var,ctg)){
        out$left  <- dat[is.element(dat[,var], unlist(strsplit(best.cut, ","))),]
        out$right <- dat[!is.element(dat[,var], unlist(strsplit(best.cut, ","))),]
      } else {
        if(cut[1]=='l'){
          out$left  <- data.frame(dat[dat[,var]<= best.cut, colnames(dat) != "new.trt"],
                                  new.trt = rep(1,n=sum(dat[,var]<= best.cut)))
          out$right <- data.frame(dat[dat[,var]> best.cut, colnames(dat) != "new.trt"],
                                  new.trt = rep(0,n=sum(dat[,var]> best.cut)))
        }else{
          out$left  <- data.frame(dat[dat[,var]<= best.cut, colnames(dat) != "new.trt"],
                                  new.trt = rep(0,n=sum(dat[,var]<= best.cut)))
          out$right <- data.frame(dat[dat[,var]> best.cut, colnames(dat) != "new.trt"],
                                  new.trt = rep(1,n=sum(dat[,var]> best.cut)))
        }  
      }
      out$info <- data.frame(node=name, size = n, n.1=n.1, n.0=n.0, trt.effect=trt.effect, var = var, 
                             vname=vname, cut.1 = unique(cut[,1]), cut.2 = paste(cut[,2], collapse = ','),
                             score=ifelse(max.score==-1e20, NA, max.score))
    } else {
      out$info <- data.frame(node=name, size = n, n.1=n.1, n.0=n.0, trt.effect=trt.effect,var=NA, 
                             vname=NA, cut.1= NA,cut.2=NA, score=NA)
    }
  }
  out 
}